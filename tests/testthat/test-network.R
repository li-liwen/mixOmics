context("network")

test_that("network works for rcc", {
  
  ## network representation for objects of class 'rcc'
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  ## create a tmp file
  tmp.file <- tempfile("network", fileext = ".jpeg")
  network.rcc.res <- network(nutri.res, comp = 1:3, cutoff = 0.6, save = "jpeg", 
                             name.save = tmp.file)
  expect_equal(names(network.rcc.res), c("gR", "M", "cutoff"))
  .expect_numerically_close(sum(network.rcc.res$M), 10.8786, digits = 3)
  unlink(tmp.file)
})

test_that("network works for spls", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  ## create a tmp file
  tmp.file <- tempfile("network", fileext = ".jpeg")
  network.spls.res <- network(toxicity.spls, comp = 1:3, cutoff = 0.8,
          color.node = c("mistyrose", "lightcyan"),
          shape.node = c("rectangle", "circle"),
          color.edge = color.spectral(100),
          lty.edge = "solid", lwd.edge =  1,
          show.edge.labels = FALSE, interactive = FALSE, save = "jpeg", 
          name.save = tmp.file)
  
  expect_equal(names(network.spls.res), c("gR", "M", "cutoff"))
  .expect_numerically_close(sum(network.spls.res$M), 45.6061, digits = 3)
  unlink(tmp.file)
})

test_that("network works when `shape.node` == 'none'", {
  data(breast.TCGA)
  X <- list(mRNA = breast.TCGA$data.train$mrna[,1:100],
            proteomics = breast.TCGA$data.train$protein[,1:100])
  Y <- breast.TCGA$data.train$subtype
  
  model.pls <- pls(X$mRNA, X$proteomics)
  model.block <- block.splsda(X, Y)
  
  tmp.file <- tempfile("network", fileext = ".jpeg")
  
  net <- network(model.pls, cutoff = 0.5,
                 shape.node = c("none", "none"), 
                 save = "jpeg", 
                 name.save = tmp.file)
  expect_equal(names(net), c("gR", "M", "cutoff"))
  .expect_numerically_close(sum(net$M), 11.17439, digits = 3)
  unlink(tmp.file)
  
  net <- network(model.block, cutoff = 0.5,
                 shape.node = c("none", "none"), 
                 save = "jpeg", 
                 name.save = tmp.file)
  expect_equal(names(net), c("gR", "M_mRNA_proteomics", "cutoff"))
  .expect_numerically_close(sum(net$M), 14.36216, digits = 3)
  unlink(tmp.file)
})

test_that("catches invalid `graph.scale` values", {
  data(breast.TCGA)
  X <- list(mRNA = breast.TCGA$data.train$mrna[,1:100],
            proteomics = breast.TCGA$data.train$protein[,1:100])
  Y <- breast.TCGA$data.train$subtype
  
  model <- block.splsda(X, Y)
  
  expect_error(network(model, cutoff = 0.5,
                       graph.scale = -1),
               "graph.scale")
  expect_error(network(model, cutoff = 0.5,
                       graph.scale = 2),
               "graph.scale")
})

test_that("catches invalid `size.node` values", {
  data(breast.TCGA)
  X <- list(mRNA = breast.TCGA$data.train$mrna[,1:100],
            proteomics = breast.TCGA$data.train$protein[,1:100])
  Y <- breast.TCGA$data.train$subtype
  
  model <- block.splsda(X, Y)
  
  expect_error(network(model, cutoff = 0.5,
                       size.node = -1),
               "size.node")
  expect_error(network(model, cutoff = 0.5,
                       size.node = 2),
               "size.node")
})

unlink(list.files(pattern = "*.pdf"))


test_that("network plot.graph parameter does not affect numerical output", {
  data("nutrimouse")
  X <- nutrimouse$gene
  Y <- nutrimouse$lipid
  
  pls.obj <- pls(X, Y)
  
  network.obj.F <- network(pls.obj, plot.graph = F)
  network.obj.T <- network(pls.obj, plot.graph = T)
  
  expect_equal(network.obj.F$M, network.obj.T$M)
})

test_that("color.node.individual accepts valid named color vector", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  tmp.file <- tempfile("network", fileext = ".jpeg")
  
  # Test with named color vector
  custom_colors <- c("gene1" = "red", "gene2" = "blue")
  net <- network(nutri.res, comp = 1:3, cutoff = 0.6, 
                 color.node.individual = custom_colors,
                 save = "jpeg", name.save = tmp.file, plot.graph = FALSE)
  
  expect_equal(names(net), c("gR", "M", "cutoff"))
  unlink(tmp.file)
})

test_that("color.node.individual errors when not a named vector", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(network(nutri.res, comp = 1:3, cutoff = 0.6,
                       color.node.individual = c("red", "blue")),
               "'color.node.individual' must be a named vector")
})

test_that("color.node.individual errors on invalid colors", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(network(nutri.res, comp = 1:3, cutoff = 0.6,
                       color.node.individual = c(gene1 = "notacolor")),
               "'color.node.individual' contains invalid color values")
})

test_that("color.node.individual empty vector treated as NULL", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  tmp.file <- tempfile("network", fileext = ".jpeg")
  
  # Should not error
  net <- network(nutri.res, comp = 1:3, cutoff = 0.6,
                 color.node.individual = c(),
                 save = "jpeg", name.save = tmp.file, plot.graph = FALSE)
  
  expect_equal(names(net), c("gR", "M", "cutoff"))
  unlink(tmp.file)
})

test_that("color.node.individual applies colors to matching nodes", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  # Get first two gene names
  gene_names <- rownames(nutrimouse$gene)[1:2]
  custom_colors <- setNames(c("hotpink", "cyan"), gene_names)
  
  tmp.file <- tempfile("network", fileext = ".jpeg")
  net <- network(nutri.res, comp = 1:3, cutoff = 0.6,
                 color.node.individual = custom_colors,
                 save = "jpeg", name.save = tmp.file, plot.graph = FALSE)
  
  # Verify colors were applied
  node_indices <- which(V(net$gR)$name %in% gene_names)
  applied_colors <- V(net$gR)$color[node_indices]
  expect_equal(sort(applied_colors), sort(c("hotpink", "cyan")))
  
  unlink(tmp.file)
})

test_that("color.node.individual non-matching names are silently ignored", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  # Use a name that doesn't exist
  custom_colors <- c("nonexistent_gene_xyz" = "red")
  
  tmp.file <- tempfile("network", fileext = ".jpeg")
  # Should not error and should return normal result
  net <- network(nutri.res, comp = 1:3, cutoff = 0.6,
                 color.node.individual = custom_colors,
                 save = "jpeg", name.save = tmp.file, plot.graph = FALSE)
  
  expect_equal(names(net), c("gR", "M", "cutoff"))
  unlink(tmp.file)
})
