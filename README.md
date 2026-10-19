# phenograph.R

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


## Claude Science Trace

I'll start by pulling the original source to confirm exactly what `Rphenograph` computes. Then I'll write a replacement that uses only CRAN binary packages and check it against the original.

The original algorithm uses exact kNN via RANN, Rcpp-computed Jaccard edge weights, and Louvain clustering on an igraph graph, with the Jaccard step being the only compiled code—one I can replace with sparse matrix algebra from Matrix. Next I'll check what's available locally and look for existing alternatives.

The environment is ready. Before writing the code, I'll check whether the course repo already has cytometry data I could test on.

No cytometry data in the repo, so I'll test using a simulated cytometry-like dataset and the iris example from the docs. I've rebuilt the implementation using only RANN, Matrix, and igraph, rewriting the Jaccard step as a chunked sparse product.

Now the validation. I'll compile the original Rcpp `jaccard_coeff` as a reference, then compare edge weights exactly and clustering agreement across random seeds.

Since the sandbox blocks compiled shared libraries, I'll port the C++ loop line-for-line into plain R as the reference implementation instead—slower, but exact in semantics.

The edge lists match exactly (594,414 edges, zero weight difference). Since Louvain in igraph is randomized, I'll now compare clusterings across 10 seeds, benchmark scaling to 200k cells, and test the Leiden option.

Perfect agreement (ARI=1) across all seeds means this simulation is too easy to differentiate implementations, so I'll add a harder test with overlapping populations. I also noticed superlinear scaling (89s→358s from 100k to 200k cells), likely from O(nnz) row-subsetting of the sparse matrix in the chunk loop—I'll profile to confirm.

kNN search is the main bottleneck (74 of ~88s at 100k cells), since RANN's kd-tree struggles in 25 dimensions just like cytofkit2's approach. I'll add an optional approximate search via RcppHNSW, validate its recall against exact kNN, and switch the Jaccard chunking to column slices.

HNSW gives a 19× speedup (3.8s vs 73s) with 99.1% recall, but in 0.03% of cases it doesn't return the point itself first, so self-removal needs to check identity rather than just drop column 1. I'll update the implementation accordingly.

I left a meaningless guard line in that edit. I'll remove it and pass the new `knn` and `n_threads` arguments through `phenograph()`.

The edges are still identical after the refactor. The 10-seed comparison on the harder data repeats the exact kNN search 20 times, so it takes a while. I'll wait for it.

All four methods agree perfectly across 10 seeds, each finding 14 of 15 simulated populations. Next I'll check scaling behavior with HNSW and the duplicate-row edge case from the iris example.

I've rewritten `Rphenograph` in about 100 lines of plain R that needs no compiler. On my test data it gives the same result as the original.

**Why it fails to install.** The clustering only uses about 20 lines of C++. But that C++ forces `cytofkit2` to be built from GitHub source, and that needs a compiler toolchain. The classic failure is "Could not find tools necessary to compile a package", and on Windows the fix requires Rtools, which needs admin rights (see this GitHub issue). The package also pulls in reticulate, VGAM, flowCore, shiny and more. The original `cytofkit` was removed from Bioconductor in version 3.9.

**What the replacement does.** It follows the same three steps as the original:
1. **Nearest neighbours:** `RANN::nn2`, the same exact search as the original, with k+1 neighbours and the cell itself removed.
2. **Jaccard weights:** the C++ loop is replaced by sparse matrix algebra from `Matrix`. If A is the cells × neighbours indicator matrix, `crossprod(A) * A` gives the number of shared neighbours for every neighbour pair. The weights are `u/(2k−u)/2`, exactly as in the C++, computed in chunks of rows to limit memory.
3. **Clustering:** `igraph::cluster_louvain`, with an optional `cluster_leiden`.

The only dependencies are `RANN`, `Matrix` and `igraph`, all CRAN packages with ready-made binaries for macOS and Windows. `RcppHNSW` is optional.

**Validation.** The sandbox won't load compiled C++, so I couldn't run `cytofkit2` itself. I compared against a line-by-line R port of its `jaccard_coeff.cpp` plus its exact pipeline.

| Check | Result |
|---|---|
| Jaccard edges, 20k simulated cells | 594,414 edges, identical pairs, maximum weight difference 0 |
| Clusters vs original, 10 seeds, 20k cells / 12 populations | ARI 1.000 every seed; both match the true labels (ARI 1.000) |
| Clusters vs original, harder overlapping data (30k / 15 populations) | ARI 1.000 with the original for exact search, fast search and Leiden; all find 14 clusters, ARI 0.99 vs truth |
| Fast search vs exact, 100k cells | 99.1% of neighbours match; 3.8 s vs 73 s |

This was tested on simulated cytometry-like data, not real CyTOF files. Running it once on one of your course datasets next to an existing `cytofkit2` install would be a worthwhile final check.

**Speed on this machine (fast search, 8 threads).** 100k cells take 16 s, 200k take 44 s and 500k take 128 s, with R using at most about 4.5 GB. The exact search is the bottleneck, taking 74 of the 88 s at 100k cells; `cytofkit2` uses the same search, so it is no faster.

**Small differences from the original:**
- **Row order:** the graph always contains every cell, in input order, so `membership()` lines up with the data rows. The original builds the graph from its edge list, so a cell without edges could drop out or end up out of order.
- **Duplicate rows:** the original simply drops the first neighbour column, assuming it is the cell itself. With duplicated rows that isn't always true; in `iris` one row keeps itself as a neighbour. The new version checks which neighbour is the cell itself.
- **Leiden pitfall:** the code sets `objective_function = "modularity"`. igraph's default for Leiden puts nearly every cell in its own cluster on these weights.
- **Seeds:** igraph's Louvain is random; use `seed=` to get the same clusters every run.

**Usage**
```r
install.packages(c("RANN", "Matrix", "igraph", "RcppHNSW"))   # all binaries
source("phenograph_lite.R")
cl <- Rphenograph(expr_matrix, k = 30)                          # same call as cytofkit2
cl <- phenograph(expr_matrix, k = 30, knn = "hnsw", n_threads = 4, seed = 42)
clusters <- igraph::membership(cl)
```

Bioconductor's `bluster` (`NNGraphParam(type = "jaccard", cluster.fun = "louvain")`) is another way to avoid compiling. It builds the shared-neighbour graph differently, though, so its clusters won't exactly match PhenoGraph's, and I didn't test it.

Files:
- [phenograph_lite.R]({{artifact:fa24fbaa-943e-431a-8cbf-834bf103c5c4}}) — drop-in replacement for `Rphenograph`
- [validation_hard.csv]({{artifact:a82b5eb3-22b7-468e-a671-bcd5782fe2bf}}) — 10-seed comparison on the overlapping data
- [validation_seeds.csv]({{artifact:06f79a21-5d35-47f5-b069-6be8166af3c9}}) — 10-seed comparison on the easier data
- [scaling_hnsw.csv]({{artifact:c305f5e4-2f00-4ff7-b194-b3fbf8ff200e}}) — run times and memory



