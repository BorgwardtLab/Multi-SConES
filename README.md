# Multi-SConES
A multi-task version of SConES, which achieves multi-task feature selection coupled with multiple network regularizers using a maximum-flow algorithm.

Please see the following paper for detailed information:
* M. Sugiyama, C.-A. Azencott, D. Grimm, Y. Kawahara, K. M. Borgwardt:
**Multi-Task Feature Selection on Multiple Networks via Maximum Flows**,
*Proceedings of the SIAM International Conference on Data Mining (SDM 2014)*, 199-207, 2014
[[PDF]](http://epubs.siam.org/doi/pdf/10.1137/1.9781611973440.23)

## Usage
To load files, type in R (without the `>`, which signifies the prompt):
```
> source("make.R")
> make()
```

To run Multi-SConES, type in R:
```
> mscones(g = g, X = X, Y = Y, lambda = lambda, eta = eta, mu = mu)
```

* Two R packages `igraph` and `glmnet` need to be installed
* `g` is a graph (in igraph format)  
* `X` is a data matrix (rows: objects, columns: features, each feature corresponds to each vertex in `g`)
* `Y` is a matrix of response vectors (rows: objects, columns: tasks)  
* `lambda`, `eta`, `mu` are parameters (they should be determined by grid-search with cross-validation)  
* output: selected features for each task  

## Example
```
> source("make.R")
> make()
> d1 <- generate.data(200, 1, seed = 1)
> d2 <- generate.data(200, 2, seed = 1)
> X <- d1$x; Y <- cbind(d1$y, d2$y)
# simulate two tasks d1$y and d2$y, and d1$x and d2$x are the same
# features from 1 to 44 are causal
> g <- generate.graph()
> res <- mscones(g = g, X = X, Y = Y)
> res
$`selected features for task 1`
+ 44/4402 vertices, from 53a17ed:
 [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
[26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44

$`selected features for task 2`
+ 44/4402 vertices, from 53a17ed:
 [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
[26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44
```


## Contact

* Author: [Mahito Sugiyama](http://mahito.info/index_e.html)   
* Affiliation: ISIR, Osaka University  
* E-Mail:  mahito@ar.sanken.osaka-u.ac.jp  
