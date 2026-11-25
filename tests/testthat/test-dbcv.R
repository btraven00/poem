set.seed(42)
d1 <- mockData()
x <- as.matrix(d1[, seq_len(2)])
y <- as.integer(d1$class)

test_that("dbcv counts noise in the denominator", {
  yn <- y
  yn[1:20] <- -1L
  with_noise <- dbcv(x, yn)$dbcv
  without <- dbcv(x[-(1:20), ], y[-(1:20)])$dbcv
  expect_equal(with_noise, without * (length(y) - 20) / length(y))
})

test_that("class-level dbcv gives NA for singleton classes", {
  lab <- factor(as.character(d1$class), levels = c(levels(d1$class), "Z"))
  lab[1] <- "Z"
  m <- getEmbeddingMetrics(x, lab, metrics = "dbcv", level = "class")
  expect_equal(nrow(m), nlevels(lab))
  expect_true(is.na(m$dbcv[m$class == "Z"]))
  ref <- getEmbeddingMetrics(x[-1, ], droplevels(lab[-1]), metrics = "dbcv",
                             level = "class")
  expect_equal(m$dbcv[m$class != "Z"], ref$dbcv)
})

test_that("non-euclidean dbcv uses the between-cluster block", {
  # in 1-D, manhattan and euclidean distances coincide
  x1 <- x[, 1, drop = FALSE]
  expect_equal(dbcv(x1, y, distance = "manhattan"),
               dbcv(x1, y, distance = "euclidean"))
})

test_that("dbcv handles clusters too small to have internal MST nodes", {
  y2 <- y
  y2[1:2] <- 99L
  expect_no_warning(r <- dbcv(x, y2))
  expect_true(all(is.finite(r$vcs)))
  y3 <- y
  y3[1:3] <- 99L
  expect_true(all(is.finite(dbcv(x, y3)$vcs)))
})

test_that("Prim MST matches igraph's when the MST is unique", {
  # Mutual reachability distances are full of ties (most MST edges equal a
  # core distance), so the MST is only unique without core distances.
  set.seed(7)
  for (n in c(2, 3, 4, 50)) {
    d <- as.matrix(dist(matrix(rnorm(n * 3), n)))
    diag(d) <- Inf
    got <- poem:::.get_internal_objects(d, numeric(n))
    g <- igraph::mst(igraph::graph_from_adjacency_matrix(
      d, mode = "undirected", weighted = TRUE, diag = FALSE))
    int <- igraph::degree(g) > 1
    w <- igraph::E(g)$weight
    ends <- igraph::ends(g, igraph::E(g), names = FALSE)
    dsc <- if (sum(int) >= 2) max(w[int[ends[, 1]] & int[ends[, 2]]]) else max(w)
    expect_equal(got$dsc, dsc)
    expect_equal(got$internal_node_inds,
                 if (any(int)) unname(which(int)) else seq_len(n))
  }
})
