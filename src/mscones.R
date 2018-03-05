##--
##    Multi-SConES
##    Copyright (C) 2014 Mahito Sugiyama
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License along
##    with this program; if not, write to the Free Software Foundation, Inc.,
##    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
##
##    Contact: Mahito Sugiyama <mahito@ar.sanken.osaka-u.ac.jp>
##
##    Please refer the following article in your published research:
##    Sugiyama, M., Azencott, C.-A., Grimm, D., Kawahara, Y., Borgwardt, K. M.:
##    Multi-Task Feature Selection on Multiple Networks via Maximum Flows,
##    SIAM International Conference on Data Mining (SDM 2014), 2014.
##--

## input:: g: a graph, X: a data matrix (rows: objects, columns: features), Y: a matrix of response vectors (rows: objects, columns: tasks)
##         lambda, eta, mu: parameters (they should be determined by grid-search with cross-validation)
## output:: selected features for each task
mscones <- function(g, C, K, X, Y, lambda = .5, eta = .25, mu=1, method="correlation", duplicate=TRUE, quiet=FALSE) {
  if (missing(K)) K <- ncol(Y) # Number of tasks
  if (missing(C)) {
    C <- NULL
    for (i in 1:K) {
      C <- c(C, abs(make.C(X, Y[,i], method)))
    }
  }

  if (duplicate) {
    if (!quiet) cat("Making a unified graph ... ")
    n <- vcount(g)
    g <- duplicate(g, K, mu)
    if (!quiet) cat("end\n")
  } else {
    n <- vcount(g) / K
  }

  if (!quiet) cat("Run SConES ... ")
  res.total <- scones(g=g, C=C, lambda=lambda, eta=eta, method=method)
  if (!quiet) cat("end\n")

  if (!quiet) cat("Post-processing ... ")
  F.list <- as.list(rep(NA, K))
  for (i in 1:K) {
    F.list[[i]] <- res.total[n * (i - 1) + 1 <= res.total & res.total <= n * i] - (i - 1) * n
  }

  ## remove isolated vertices
  i <- 1
  for (F in F.list) {
    g.sub <- induced_subgraph(g, as.integer(F))
    g.d <- degree(g.sub)
    F.list[[i]] <- F[g.d > 0]
    names(F.list)[i] <- paste("selected features for task ", i, sep="")
    i <- i + 1
  }
  if (!quiet) cat("end\n")
  F.list
}

## cross-validation of Multi-SConES
## R package  "glmnet" is needed
cv.mscones <- function(g, X, Y, lambda.list, eta.list, mu.list, method="correlation", quiet=FALSE) {
  if (missing(lambda.list)) lambda.list <- seq(0.1, 1, 0.2)
  if (missing(mu.list)) mu.list <- 1
  K <- ncol(Y) # Number of tasks
  n <- nrow(X)
	
  res.best <- rep(Inf, K)
  lambda.best <- numeric(K)
  eta.best <- numeric(K)
  mu.best <- numeric(K)
  features.best <- as.list(rep(NA, K))

  C <- NULL
  for (i in 1:K) {
    C <- c(C, abs(make.C(X, Y[,i], method)))
  }
  if (missing(eta.list)) {
    C.s <- sort(C[1:ncol(X)], decreasing=TRUE)
    eta.min <- C.s[ncol(X) * 0.3]
    eta.max <- C.s[ncol(X) * 0.02]
    eta.list <- seq(eta.min, eta.max, by = (eta.max - eta.min) / 10)
  }

  ## decide foldid
  fold.num <- 10
  foldid.list <- as.list(rep(NA, K))
  n.each <- ceiling(n / fold.num)
  for (i in 1:K) {
    tmp <- rep(1:fold.num, n.each)[1:n]
    foldid.list[[i]] <- sample(tmp, n)
  }

  if (!quiet) cat("parameter selection start\n")
  for (mu in mu.list) {
    g.d <- duplicate(g, K, mu)
    for (lambda in lambda.list) {
      for (eta in eta.list) {
        if (!quiet) cat("mu = ", mu, " lambda = ", lambda, ", eta = ", eta, "\n", sep="")
        F.list <- mscones(g=g.d, C=C, K=K, lambda=lambda, eta=eta, mu=mu, method=method, duplicate=FALSE)
        res.score <- numeric(K)
        for (i in 1:K) {
          res.score[i] <- cv.mscones.eval(X=X[,F.list[[1]]], Y=Y[,i], foldid=foldid.list[[i]], method="rreg")
          if (!quiet) cat("------------ cvm: ", res.score[i], " -----------\n", sep="")
          if (res.score[i] < res.best[i]) {
            res.best[i] <- res.score[i]
            lambda.best[i] <- lambda
            eta.best[i] <- eta
            mu.best[i] <- mu
            features.best[[i]] <- F.list[[i]]
          }
        }
      }
    }
    if (!quiet) cat("\n")
  }

  for (i in 1:K) {
    if (!quiet) cat("For task: ", i, "\n", sep="")
    if (!quiet) cat("  best lambda = ", lambda.best[i], ", best eta = ", eta.best[i], "best mu = ", mu.best[i], "\n")
  }
  res <- list(res.best, lambda.best, eta.best, mu.best, features.best)
  names(res) <- c("cvm", "lambda", "eta", "mu", "features")
  return(res)
}
cv.mscones.eval <- function(X, Y, foldid, foldnum, method="enet") {
  if (method == "enet") {
    res <- ifelse(class(try( res <- cv.glmnet(X, Y, alpha=0.1, foldid=foldid) )) == "try-error", Inf, min(res$cvm))
  } else if (method == "rreg") {
    res <- ifelse(class(try( res <- cv.glmnet(X, Y, alpha=0, foldid=foldid) )) == "try-error", Inf, min(res$cvm))
  } else if (method == "lm") {
    res <- lm(Y ~ X)
    res <- summary(res)$adj.r.squared
  }
  return(res)
}


## SConES, which is proposed in:
## Azencott, C.-A., Grimm, D., Sugiyama, M., Kawahara, Y., Borgwardt, K. M.:
## Efficient Network-guided Multi-locus Association Mapping with Graph Cut,
## Bioinformatics, 29(13), i171—i179, 2013
## package "igraph" is needed
scones <- function(g, C, X, Y, lambda = .5, eta = .5, method="correlation") {
  if (missing(C)) C <- abs(make.C(X, Y, method))

  ## convert a given graph for max flow
  n <- vcount(g)
  g <- convert(g, C, lambda, eta)
  ## apply maxflow
  res <- graph.maxflow(g, V(g)[n + 1], V(g)[n + 2], capacity=E(g)$weight)
  if(length(res$partition1) < length(res$partition2)) {
    features <- res$partition1[-length(res$partition1)]
  } else {
    features <- res$partition2[-length(res$partition2)]
  }
  ## remove isolated vertices
  g.sub <- induced_subgraph(g, features)
  g.d <- degree(g.sub)
  return(features[g.d > 0])
}
## for cross-validation
cv.scones <- function(g, X, Y, k=Inf, lambda.list, eta.list, method="correlation", quiet=FALSE) {
  if (missing(lambda.list)) lambda.list <- seq(0.1, 1, 0.2) #lambda.list <- 0.5
  ## if (missing(eta.list)) eta.list <- seq(0.1, 0.5, 0.1)
  res.best <- Inf
  fold.num <- 10

  ## decide foldid
  n <- nrow(X)
  n.each <- ceiling(n / 10)
  foldid <- rep(1:fold.num, n.each)[1:n]
  foldid <- sample(foldid, n)

  ## make C
  C <- abs(make.C(X, Y, method))
  if (missing(eta.list)) {
    C.s <- sort(C, decreasing=TRUE)
    eta.min <- C.s[ncol(X) * 0.3]
    eta.max <- C.s[ncol(X) * 0.02]
    eta.list <- seq(eta.min, eta.max, by = (eta.max - eta.min) / 10)
  }

  if (!quiet) cat("cross-validation start\n")
  for (lambda in lambda.list) {
    res.best.each <- Inf
    for (eta in eta.list) {
      if (!quiet) cat("lambda: ", round(lambda, digits=3), ", eta: ", round(eta, digits=3), sep="")
      features <- scones(g=g, C=C, lambda=lambda, eta=eta, method=method)
      if (length(features) > k) res.cvm <- Inf
      else res.cvm <- cv.mscones.eval(X=X[,features], Y=Y, foldid=foldid, foldnum=fold.num, method="rreg")
      if (!quiet) cat(" csv: ", round(res.cvm, digits=3), "\n", sep="")

      if (res.cvm < res.best) {
        res.best <- res.cvm
        lambda.best <- lambda
        eta.best <- eta
      }
      if (res.cvm <= res.best.each) {
        res.best.each <- res.cvm
      }
    }
    if (!quiet) cat("\n")
  }

  if (!quiet) cat("best lambda = ", lambda.best, ", best eta = ", eta.best, "\n")
  res <- scones(g=g, X=X, Y=Y, lambda=lambda.best, eta=eta.best, method=method)
  if (!quiet) cat("end\n")
  res <- list(res.best, lambda.best, eta.best, res)
  names(res) <- c("Best cvm", "Best lambda", "Best eta", "Selected features")
  return(res)
}


## duplicate a given graph
## (used in mscones)
duplicate <- function(g, K, mu=1) {
  n <- vcount(g)
  E <- get.edges(g, 1:ecount(g))
  E.vec <- rep(as.vector(t(E)), K)
  E.vec <- E.vec + rep(0:(K - 1), each=ecount(g) * 2) * vcount(g)
  E <- matrix(E.vec, ncol=2, byrow=TRUE)

  idx <- which(degree(g) > 0)

  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      E <- rbind(E, cbind(idx + (n * (i - 1)), idx + (n * (j - 1))))
    }
  }
  g.d <- graph(n = n * K, edges=t(E), directed=FALSE)
  E(g.d)$weight <- c(rep(E(g)$weight, K), rep(mu, sum(degree(g) > 0) * K * (K - 1) / 2))

  return(g.d)
}


## convert a given graph for max flow
## (used in scones)
convert <- function(g, c, lambda, eta) {
  n <- vcount(g)
  s <- n + 1
  t <- n + 2
  E(g)$weight <- E(g)$weight * lambda
  g <- add.vertices(g, 2)
  s.idx <- c > eta
  t.idx <- c < eta
  s.E <- matrix(c(rep(s, n), 1:n), nrow=2, byrow=TRUE)[,s.idx]
  t.E <- matrix(c(rep(t, n), 1:n), nrow=2, byrow=TRUE)[,t.idx]
  g <- g + edges(s.E, weight = (c - eta)[s.idx])
  g <- g + edges(t.E, weight = (eta - c)[t.idx])
  return(g)
}

## make a vector of correlations from a given data matrix X and Y
## (used in scones)
make.C <- function(X, Y, method="correlation", quiet=FALSE) {
  if (missing(X) || missing(Y)) stop("Input C or a design matrix X and Y!!")
  method <- match.arg(method, choices=c("correlation", "inner"))
  if (method == "correlation") {
    if (!quiet) cat("Determine C by correlation coefficient ... ")
    C <- as.vector(cor(X, Y))
    C[is.na(C)] <- 0
    if (!quiet) cat("end\n")		
  } else if (method == "inner") {
    if (!quiet) cat("Determine C by inner product ... ")
    Y <- (Y - mean(Y)) / sd(Y)
    C <- Y %*% X
    if (!quiet) cat("end\n")		
  }
  return(C)
}


##------------------------------------------------------------------------------------#
## Simulated data (From Li and Li; 2008)
## For example:
## d <- generate.data(100, 1)
generate.data <- function(size=200, model=1, cor=0.7, g.sd=sqrt(0.51), seed) {
  beta <- switch(model,
                 c(5, rep(5/sqrt(10), 10), -5, rep(-1 * 5 / sqrt(10), 10),
                   3, rep(3/sqrt(10), 10), -3, rep(-1 * 3 / sqrt(10), 10), numeric(2156)),
                 c(5, rep(-1 * 5/sqrt(10), 3), rep(5/sqrt(10), 7), -5, rep(5 / sqrt(10), 3), rep(-1 * 5/sqrt(10), 7),
                   3, rep(-1 * 3/sqrt(10), 3), rep(3/sqrt(10), 7), -3, rep(3/sqrt(10), 3), rep(-1 * 3 / sqrt(10), 7), numeric(2156)),
                 c(5, rep(5/10, 10), -5, rep(-1 * 5 / 10, 10),
                   3, rep(3/10, 10), -3, rep(-1 * 3 / 10, 10), numeric(2156)),
                 c(5, rep(-1 * 5/10, 3), rep(5/10, 7), -5, rep(5 / 10, 3), rep(-1 * 5/10, 7),
                   3, rep(-1 * 3/10, 3), rep(3/10, 7), -3, rep(3/10, 3), rep(-1 * 3 / 10, 7), numeric(2156))
                 )
  tr.idx <- seq(1, 2200, 11)

  X <- matrix(numeric(size * 2200), nrow=size)
  if (!missing(seed)) set.seed(seed)
  X[,tr.idx] <- matrix(rnorm(size * 200), nrow=size)

  for (i in 1:200) {
    gene.idx <- (tr.idx[i] + 1):(tr.idx[i] + 10)
    for (j in 1:size) {
      if (!missing(seed)) set.seed(seed + j)
      X[j,gene.idx] <- rnorm(10, X[j,tr.idx[i]] * cor, g.sd)
    }
  }
  sigma2 <- sum(beta^2) / 4
  set.seed(runif(1))
  Y <- as.vector(X %*% beta) + rnorm(size, 0, sqrt(sigma2))
  d <- list(X, Y, 1:44, beta)
  names(d) <- c("x", "y", "betaIndex", "beta")
  return(d)
}

## Simulated graph (From Li and Li; 2008)
generate.graph <- function() {
  tr.idx <- seq(1, 2200, 11)
  e1 <- rep(tr.idx, each=10)
  tmp <- rep(1:10, 200)
  e2 <- e1 + tmp
  g <- graph(n=2200, edges=as.vector(rbind(e1,e2)), directed=FALSE)
  E(g)$weight <- 1
  return(g)
}
