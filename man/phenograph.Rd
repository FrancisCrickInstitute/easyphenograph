\name{phenograph}
\alias{phenograph}
\title{Using phenograph}
\usage{
phenograph(data, k = 30, method = c("louvain", "leiden"),
  resolution = 1, knn = c("exact", "hnsw"), n_threads = 1,
  seed = 42, verbose = TRUE)
}
\description{
Drop-in replacement for cytofkit2::Rphenograph() using only CRAN packages
that ship pre-built binaries for macOS / Windows: RANN, Matrix, igraph.
No Rcpp / compiler toolchain required.

Algorithm (identical to cytofkit2::Rphenograph):
  1. exact kNN (k neighbours, self excluded) with RANN::nn2 (kd-tree)
  2. for every kNN edge i -> j: u = |N(i) ∩ N(j)|, w = u / (2k - u)   (Jaccard)
     edges with u = 0 are dropped; each directed edge carries w / 2, so a
     mutual pair ends up with total weight w and a one-way pair with w / 2
  3. Louvain (or Leiden) community detection on the undirected weighted graph

Step 2 replaces the Rcpp jaccard_coeff() with sparse algebra:
  A = n x n kNN indicator matrix;  (A %*% t(A))[i, j] = |N(i) ∩ N(j)|
  multiplied element-wise by A so only kNN edges are kept, done in row
  chunks to bound memory.
}
\examples{
install.packages(c("RANN", "Matrix", "igraph", "RcppHNSW"))   # all binaries
source("phenograph_lite.R")
cl <- Rphenograph(expr_matrix, k = 30)                          # same call as cytofkit2
cl <- phenograph(expr_matrix, k = 30, knn = "hnsw", n_threads = 4, seed = 42)
clusters <- igraph::membership(cl)
# or
clusters <- cl$membership
}
