#' @title Compute Pair to Pair Distances
#' @description Compute the pairwise distances between points in matrix `X`.
#' @param X Numeric matrix.
#' @param distance String specifying the metric to compute the distances.
#' @return Numeric matrix of pairwise distances with self-distances set to 
#' `Inf`.
#' @keywords internal
.compute_pair_to_pair_dists <- function(X, distance="euclidean") {
  if(distance == "sqeuclidean"){
    dists <- as.matrix(dist(X, method = "euclidean"))
    dists <- dists^2
  }else{
    dists <- as.matrix(dist(X, method = distance))
  }
  dists[dists < 1e-12] <- 1e-12
  diag(dists) <- Inf
  return(dists)
}


#' @title Get Internal Objects
#' @description Builds the minimum spanning tree of the mutual reachability
#' graph of one cluster (Prim's algorithm on the dense distance matrix) and 
#' returns its internal nodes and the density sparseness (largest internal 
#' edge weight). Ties in mutual reachability distance are broken by the raw
#' distance.
#' @param dists Numeric matrix of within-cluster distances, `Inf` on the 
#' diagonal.
#' @param core_dists Numeric vector of core distances.
#' @return A list with the indices of internal nodes and `dsc`.
#' @importFrom Matrix Diagonal
#' @keywords internal
.get_internal_objects <- function(dists, core_dists) {
  # Rows of the mutual reachability matrix are built on the fly, so only
  # `dists` is held in memory. ponytail: O(n^2) R-level Prim, n iterations of
  # vector ops; move to C++ if single clusters get much beyond ~50k cells.
  # Mutual reachability distances are heavily tied (most MST edges equal a
  # core distance), so the MST is not unique and the internal nodes depend on
  # tie-breaking. Ties are broken by the raw distance, which makes the MST
  # unique for continuous data and independent of the order of the samples.
  n <- nrow(dists)
  mrd <- function(j) pmax(dists[, j], core_dists[j], core_dists)
  in_tree <- logical(n)
  in_tree[1] <- TRUE
  best <- mrd(1)
  best_raw <- dists[, 1]
  best[1] <- Inf
  parent <- rep(1L, n)
  from <- to <- integer(n - 1)
  w <- numeric(n - 1)
  for (k in seq_len(n - 1)) {
    cand <- which(best == min(best))
    j <- cand[which.min(best_raw[cand])]
    from[k] <- parent[j]
    to[k] <- j
    w[k] <- best[j]
    in_tree[j] <- TRUE
    best[j] <- Inf
    dj <- mrd(j)
    rj <- dists[, j]
    upd <- !in_tree & (dj < best | (dj == best & rj < best_raw))
    best[upd] <- dj[upd]
    best_raw[upd] <- rj[upd]
    parent[upd] <- j
  }
  internal <- tabulate(c(from, to), n) > 1
  # As in the reference implementation: tiny clusters have no (or a single)
  # internal node, so fall back to all nodes / the whole MST.
  dsc <- if (sum(internal) >= 2) max(w[internal[from] & internal[to]]) 
         else max(w)
  internal_node_inds <- if (any(internal)) which(internal) else seq_len(n)
  list(internal_node_inds = internal_node_inds, dsc = dsc)
}

#' @title Compute Cluster Core Distance
#' @description Computes the core distance for each point in a cluster.
#' @param dists Numeric matrix of distances.
#' @param d Integer, the dimensionality.
#' @return Numeric vector of core distances for each point.
#' @keywords internal
.compute_cluster_core_distance <- function(dists, d) {
  n <- nrow(dists)
  core_dists <- rowSums(dists^(-d)) / (n - 1)
  core_dists <- core_dists^(-1/d)
  return(as.numeric(core_dists))
}

#' @title Check Duplicated Samples
#' @description Checks for duplicated samples in matrix `X`.
#' @param X Numeric matrix of samples.
#' @param threshold Numeric, the distance threshold to consider samples as 
#' duplicates.
#' @keywords internal
#' @return None
.check_duplicated_samples <- function(X, threshold = 1e-9) {
  if (nrow(X) <= 1) {
    return()
  }
  dists <- as.matrix(dist(X))
  if (any(dists < threshold)) {
    stop("Duplicated samples have been found in X.")
  }
}

#' @title Convert Singleton Clusters to Noise
#' @description Converts clusters containing a single instance to noise.
#' @param labels Integer vector of cluster IDs.
#' @param noise_id Integer, the ID for noise.
#' @keywords internal
#' @return Integer vector with singleton clusters converted to noise.
.convert_singleton_clusters_to_noise <- function(labels, noise_id) {
  cluster_sizes <- table(labels)
  singleton_clusters <- names(cluster_sizes[cluster_sizes == 1])
  if (length(singleton_clusters) == 0) {
    return(labels)
  }
  return(ifelse(labels %in% singleton_clusters, noise_id, labels))
}

#' @title Compute Cross Distances
#' @description Computes the cross distances between two sets of indices in matrix `X`.
#' @param X Numeric matrix of samples.
#' @param inds_a Integer vector of indices for the first set.
#' @param inds_b Integer vector of indices for the second set.
#' @param distance String specifying the distance metric.
#' @return Numeric matrix of cross distances.
#' @keywords internal
.compute_cross_dists <- function(X, inds_a, inds_b, distance) {
  if (length(inds_a) == 0 || length(inds_b) == 0) {
    return(matrix(numeric(0), nrow = length(inds_a), ncol = length(inds_b)))
  }
  
  if (distance == "euclidean") {
    a <- X[inds_a, , drop = FALSE]
    b <- X[inds_b, , drop = FALSE]
    aa <- rowSums(a^2)
    bb <- rowSums(b^2)
    ab <- tcrossprod(a, b)
    dists <- sqrt(pmax(outer(aa, bb, "+") - 2 * ab, 0))
  } else if (distance == "sqeuclidean") {
    a <- X[inds_a, , drop = FALSE]
    b <- X[inds_b, , drop = FALSE]
    aa <- rowSums(a^2)
    bb <- rowSums(b^2)
    ab <- tcrossprod(a, b)
    dists <- pmax(outer(aa, bb, "+") - 2 * ab, 0)
  } else {
    # General method for other distances
    dists <- as.matrix(dist(rbind(X[inds_a, , drop = FALSE], 
                            X[inds_b, , drop = FALSE]), 
                         method = distance))
    dists <- dists[seq_along(inds_a), length(inds_a) + seq_along(inds_b), 
                   drop = FALSE]
  }
  dists[dists < 1e-12] <- 1e-12
  return(dists)
}

#' @title Calculate DBCV Metric 
#' @description Compute the DBCV (Density-Based Clustering Validation) metric.
#' 
#' @param X Numeric matrix of samples.
#' @param labels Integer vector of cluster IDs.
#' @param distance String specifying the distance metric. `"sqeuclidean"`, or 
#' possible `method` in `stats::dist()`. By default `"euclidean"`.
#' @param noise_id Integer, the cluster ID in `y` for noise (default `-1`). 
#' @param check_duplicates Logical flag to check for duplicate samples.
#' @param use_igraph_mst Ignored; kept for backward compatibility. The MST is
#' now computed with a dense Prim's algorithm.
#' @param BPPARAM BiocParallel params for multithreading (default none)
#' @param ... Ignored
#' 
#' @importFrom BiocParallel SerialParam bplapply
#' 
#' @return A list:
#' \item{vcs}{Numeric vector of validity index for each cluster, named by
#' cluster ID. Noise and singleton clusters are not included.}
#' \item{dbcv}{Numeric value representing the overall DBCV metric.}
#'
#' @importFrom utils combn
#' 
#' @references Davoud Moulavi, et al. 2014; 10.1137/1.9781611973440.96.
#' 
#' @details This implementation will not fully reproduce the results of other 
#' existing implementations (e.g. \url{https://github.com/FelSiq/DBCV}).
#' Mutual reachability distances are heavily tied, so the Minimum Spanning
#' Tree, and hence the internal nodes, are not unique. Here ties are broken by
#' the raw distance, which makes the result independent of sample order.
#' @export
#' @examples
#' data(noisy_moon)
#' data <- noisy_moon
#' dbcv(data[, c("x", "y")], data$kmeans_label)
#' dbcv(data[, c("x", "y")], data$hdbscan_label)
dbcv <- function(X, labels, distance = "euclidean", noise_id = -1, 
                 check_duplicates = FALSE,
                 use_igraph_mst = TRUE, BPPARAM = BiocParallel::SerialParam(), 
                 ...) {
  X <- as.matrix(X)
  labels <- as.integer(labels)
  n <- nrow(X)
  if (n != length(labels)) {
    stop("Mismatch in number of samples and cluster labels.")
  }
  
  labels <- .convert_singleton_clusters_to_noise(labels, noise_id = noise_id)
  non_noise_inds <- labels != noise_id
  X <- X[non_noise_inds, , drop = FALSE]
  labels <- labels[non_noise_inds]
  
  if (length(labels) == 0) {
    return(list(vcs = numeric(0), dbcv = 0.0))
  }
  
  if (check_duplicates) {
    .check_duplicated_samples(X)
  }
  
  cluster_ids <- sort(unique(labels))
  n_clusters <- length(cluster_ids)
  cls_inds <- lapply(cluster_ids, function(id) which(labels == id))
  
  # Initialize storage
  dscs <- numeric(n_clusters)
  min_dspcs <- rep(Inf, n_clusters)
  internal_objects_per_cls <- vector("list", n_clusters)
  internal_core_dists_per_cls <- vector("list", n_clusters)
  
  # Process each cluster
  density_sparseness_results <- BiocParallel::bplapply(
    seq_along(cls_inds),
    function(i) {
      inds <- cls_inds[[i]]
      X_cls <- X[inds, , drop = FALSE]
      dists_intra <- .compute_pair_to_pair_dists(X_cls, distance = distance)
      core_dists <- .compute_cluster_core_distance(dists_intra, d = ncol(X))
      internal_objects <- .get_internal_objects(dists_intra, core_dists)
      dsc <- internal_objects$dsc
      internal_node_inds_global <- inds[internal_objects$internal_node_inds]
      internal_core_dists <- core_dists[internal_objects$internal_node_inds]
      list(
        dsc = dsc,
        internal_node_inds_global = internal_node_inds_global,
        internal_core_dists = internal_core_dists
      )
    },
    BPPARAM = BPPARAM
  )
  
  # Extract results
  for (i in seq_along(density_sparseness_results)) {
    dscs[i] <- density_sparseness_results[[i]]$dsc
    internal_objects_per_cls[[i]] <- density_sparseness_results[[i]]$internal_node_inds_global
    internal_core_dists_per_cls[[i]] <- density_sparseness_results[[i]]$internal_core_dists
  }
  
  # Process cluster pairs
  if (n_clusters > 1) {
    pairs <- utils::combn(n_clusters, 2, simplify = FALSE)
    density_separation_results <- BiocParallel::bplapply(
      pairs,
      function(pair) {
        i <- pair[1]
        j <- pair[2]
        inds_i <- internal_objects_per_cls[[i]]
        inds_j <- internal_objects_per_cls[[j]]
        if (length(inds_i) == 0 || length(inds_j) == 0) {
          return(list(cls_i = i, cls_j = j, dspc_ij = Inf))
        }
        dists_inter <- .compute_cross_dists(
          X, inds_i, inds_j, distance
        )
        core_i <- internal_core_dists_per_cls[[i]]
        core_j <- internal_core_dists_per_cls[[j]]
        mat_core_i <- matrix(core_i, nrow = length(core_i), ncol = length(core_j), byrow = FALSE)
        mat_core_j <- matrix(core_j, nrow = length(core_i), ncol = length(core_j), byrow = TRUE)
        sep <- pmax(dists_inter, mat_core_i, mat_core_j)
        dspc_ij <- min(sep)
        list(cls_i = i, cls_j = j, dspc_ij = dspc_ij)
      },
      BPPARAM = BPPARAM
    )
    
    for (result in density_separation_results) {
      i <- result$cls_i
      j <- result$cls_j
      min_dspcs[i] <- min(min_dspcs[i], result$dspc_ij)
      min_dspcs[j] <- min(min_dspcs[j], result$dspc_ij)
    }
  }
  
  # Compute validity scores
  vcs <- (min_dspcs - dscs) / pmax(min_dspcs, dscs)
  vcs[is.nan(vcs)] <- 0.0
  names(vcs) <- cluster_ids
  # Noise (incl. singleton clusters) counts in the denominator (Moulavi et al.)
  dbcv <- sum(vcs * table(labels)) / n
  return(list(vcs = vcs, dbcv = dbcv))
}