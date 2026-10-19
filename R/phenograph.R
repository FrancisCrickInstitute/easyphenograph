find_neighbors <- function(data, k, knn = c("exact", "hnsw"), n_threads = 1) {
  ## k + 1 so the point itself can be removed. cytofkit2 drops column 1
  ## blindly; with duplicated rows (or approximate search) column 1 is not
  ## always the point itself, so remove self explicitly and otherwise drop
  ## the furthest neighbour.
  knn <- match.arg(knn)
  nn <- if (knn == "exact") {
    RANN::nn2(data, data, k = k + 1, treetype = "kd", searchtype = "standard")$nn.idx
  } else {
    if (!requireNamespace("RcppHNSW", quietly = TRUE))
      stop("knn = 'hnsw' needs install.packages('RcppHNSW')")
    RcppHNSW::hnsw_knn(data, k = k + 1, distance = "euclidean",
                       M = 16, ef = 100, n_threads = n_threads)$idx
  }
  self <- nn == seq_len(nrow(nn))
  no_self <- rowSums(self) == 0
  self[no_self, k + 1] <- TRUE
  matrix(t(nn)[!t(self)], ncol = k, byrow = TRUE)
}

jaccard_edges <- function(nn, chunk_size = 20000) {
  n <- nrow(nn); k <- ncol(nn)
  A <- Matrix::sparseMatrix(i = rep(seq_len(n), times = k), j = as.vector(nn),
                            x = 1, dims = c(n, n))
  At <- Matrix::t(A)                    # column c of At = neighbour set of cell c
  starts <- seq(1, n, by = chunk_size)
  out <- lapply(starts, function(s) {
    rows <- s:min(s + chunk_size - 1, n)
    Atb <- At[, rows, drop = FALSE]     # column slice: cheap for CSC storage
    Ab <- Matrix::t(Atb)
    shared <- Matrix::crossprod(Atb, At) * Ab   # |N(i) ∩ N(j)|, kNN edges only
    tr <- Matrix::summary(as(shared, "TsparseMatrix"))
    tr <- tr[tr$x > 0, ]
    cbind(from = rows[tr$i], to = tr$j, weight = tr$x / (2 * k - tr$x) / 2)
  })
  out <- do.call(rbind, out)
  out[order(out[, "from"], out[, "to"]), , drop = FALSE]
}

phenograph <- function(data, k = 30, method = c("louvain", "leiden"),
                       resolution = 1, knn = c("exact", "hnsw"), n_threads = 1,
                       seed = 42, verbose = TRUE) {
  method <- match.arg(method)
  knn <- match.arg(knn)
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data))
    stop("Wrong input data, should be a numeric data frame or matrix!")
  if (k < 1) stop("k must be a positive integer!")
  if (k > nrow(data) - 2) stop("k must be smaller than the total number of points!")
  if (!is.null(seed)) set.seed(seed)       # igraph's Louvain/Leiden use R's RNG
  msg <- function(...) if (verbose) message(...)

  msg("phenograph: ", nrow(data), " cells x ", ncol(data), " markers, k = ", k)
  t1 <- system.time(nn <- find_neighbors(data, k, knn = knn, n_threads = n_threads))
  msg("  kNN search (", knn, ") .... ", round(t1[3], 2), " s")
  t2 <- system.time(links <- jaccard_edges(nn))
  msg("  Jaccard weights ....... ", round(t2[3], 2), " s")

  ## All n vertices are created in row order, so membership() always lines up
  ## with the input rows (cytofkit2 builds the graph from the edge list, which
  ## can silently drop or reorder cells that have no qualifying edge).
  t3 <- system.time({
    g <- igraph::make_empty_graph(n = nrow(data), directed = FALSE)
    g <- igraph::add_edges(g, t(links[, 1:2]), weight = links[, 3])
    g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum"))
  })
  msg("  graph build ........... ", round(t3[3], 2), " s")

  t4 <- system.time({
    community <- if (method == "louvain") {
      igraph::cluster_louvain(g, resolution = resolution)
    } else {
      ## objective_function must be "modularity": igraph's default (CPM with
      ## resolution 1) returns ~one cluster per cell on Jaccard weights.
      igraph::cluster_leiden(g, objective_function = "modularity",
                             resolution = resolution, n_iterations = -1)
    }
  })
  msg("  ", method, " clustering ..... ", round(t4[3], 2), " s")
  msg("  modularity = ", round(igraph::modularity(g, igraph::membership(community)), 4),
      "; clusters = ", length(unique(igraph::membership(community))))
  attr(community, "graph") <- g
  community
}

## Backwards-compatible alias with the cytofkit2 signature
Rphenograph <- function(data, k = 30) phenograph(data, k = k)
