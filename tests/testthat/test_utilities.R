txdb <- TxDb.Drerio.UCSC.danRer10.refGene
org <- org.Dr.eg.db
extPath <- system.file('extdata', package='annoLinker')
## load peaks
peaks <- rtracklayer::import(file.path(extPath, 'peaks.bed'))
## load interactions
interactions <- rtracklayer::import(file.path(extPath, 'interaction.bedpe'))
## load annotation data
annoData <- suppressMessages(genes(txdb))

test_that("validate_inputs_graph warks not correct", {
  expect_error(validate_inputs_graph(peaks=letters),
               "'peaks' must be a GRanges object")
  expect_error(validate_inputs_graph(peaks=peaks, annoData=letters),
               "'annoData' must be an annoGR or GRanges object")
  expect_error(validate_inputs_graph(peaks=peaks, annoData=GRanges()),
               "'annoData' is empty.")
  expect_error(validate_inputs_graph(peaks=peaks,
                                     annoData=GRanges('chr1', IRanges(1, 2))),
               "names of 'annoData' is missing!")
  expect_error(
    validate_inputs_graph(peaks=peaks,
                          annoData=GRanges('1',
                                           IRanges(1, 2,
                                                   names = 'a'))),
    "Please check the seqlevels style of your 'peaks' and 'annoData'")
  expect_error(validate_inputs_graph(peaks=peaks, annoData=annoData,
                                     interactions = letters),
               "'interactions' must be an GInteractions object")
  expect_error(validate_inputs_graph(
    peaks=peaks, annoData=annoData,
    interactions = GInteractions(GRanges('1', IRanges(1, 2)),
                                GRanges('1', IRanges(3, 4)))),
    "Please check the seqlevels style of your 'peaks' and 'interactions'")
  expect_error(validate_inputs_graph(peaks=peaks, annoData=annoData,
                                     interactions = interactions,
                                     bindingRegion = 1),
               "'bindingRegion' must have exactly 2 elements")
  expect_error(validate_inputs_graph(peaks=peaks, annoData=annoData,
                                     interactions = interactions,
                                     bindingRegion = c(-5, 5),
                                     interactionDistanceRange=1),
               "'interactionDistanceRange' must have exactly 2 elements")
  expect_error(validate_inputs_graph(peaks, letters),
               "'annoData' must be an annoGR or GRanges object")
  expect_error(validate_inputs_graph(peaks, letters),
               "'annoData' must be an annoGR or GRanges object")
})

test_that("filterInteractions works not correct", {
  testInteraction <- GInteractions(
    GRanges('1', IRanges(seq(1, 5), width = 1)),
    GRanges('1', IRanges(seq(15, 11), width = 1))
  )
  fi <- filterInteractions(interactions=testInteraction,
                           interactionDistanceRange=c(7, 12))
  expect_true(all(pairdist(fi)>=7 & pairdist(fi)<=12))
})

test_that("get_annotation_regions works not correct", {
  anno <- GRanges('1', IRanges(c(3, 7, 9), width = c(1, 2, 3)),
                  strand = c('+', '*', '-'))
  body <- get_annotation_regions(
    annoData=anno,
    bindingType='body',
    bindingRegion=c(-2, 1))
  expect_equal(start(body), start(anno)-c(2, 2, 1))
  expect_equal(end(body), end(anno)+c(1, 1, 2))
  promoter <- get_annotation_regions(
    annoData=anno,
    bindingType='startSite',
    bindingRegion=c(-2, 1))
  expect_equal(start(promoter)[c(1, 2)], start(anno)[c(1, 2)]-2)
  expect_equal(end(promoter)[c(3)], end(anno)[c(3)]+2)
  expect_equal(end(promoter)[c(1, 2)], start(anno)[c(1, 2)])
  terminator <- get_annotation_regions(
    annoData=anno,
    bindingType='endSite',
    bindingRegion=c(-2, 1))
  expect_equal(start(terminator)[c(1, 2)], end(anno)[c(1, 2)]-2)
  expect_equal(start(terminator)[c(3)], end(anno)[c(3)]-2)
  expect_equal(end(terminator), end(anno))
})

test_that("build_interaction_graph or detect_clusters works not correct", {
  # Create three groups with different structures
  set.seed(123)

  # Group 1: dense connections
  g1 <- sample_gnp(6, 0.7)

  # Group 2: moderate connections
  g2 <- sample_gnp(6, 0.5)

  # Group 3: sparse connections
  g3 <- sample_gnp(6, 0.3)

  # Relabel vertices so they don't overlap
  g2 <- set_vertex_attr(g2, "name", value = paste0("g2_", 1:vcount(g2)))
  g3 <- set_vertex_attr(g3, "name", value = paste0("g3_", 1:vcount(g3)))
  g1 <- set_vertex_attr(g1, "name", value = paste0("g1_", 1:vcount(g1)))

  # Combine the graphs
  g <- disjoint_union(g1, g2, g3)

  # Add a weak connection between groups
  g <- add_edges(g, c("g1_1", "g2_1", "g2_6", "g3_1"))

  anchors <- as.list(as.data.frame(as_edgelist(g)))
  methods <- c("components", "louvain", "walktrap", "infomap")
  res <- lapply(methods, build_interaction_graph, anchors=anchors, weight=NULL)
  if(interactive()){
    ## help to understand the process
    par(mfrow=c(2,2), mar=c(1,1,2,1))
    plot_ <- function(v, main){
      ids <- v$clusters$cluster_id
      names(ids) <- v$cluster$anchor_id
      ids <- ids[V(g)$name]
      plot(g, vertex.color=ids, main=main)
    }
    plot_(res[[1]], main="components()")
    plot_(res[[2]], main="Louvain")
    plot_(res[[3]], main="Walktrap")
    plot_(res[[4]], main="Infomap")
  }
  # components can not split the
  # | Method               | Type             | Detects                     | Can split connected graph? | Tends to find                        |
  # | -------------------- | ---------------- | --------------------------- | -------------------------- | ------------------------------------ |
  # | `components()`       | Structural       | Disconnected subgraphs only | No                         | Whole connected piece                |
  # | `cluster_infomap()`  | Flow-based       | Modules via info theory     | Yes                        | Large, information-dense communities |
  # | `cluster_louvain()`  | Modularity-based | Density-based communities   | Yes                        | Balanced, modular clusters           |
  # | `cluster_walktrap()` | Random-walk      | Local communities           | Yes                        | Many small, local clusters           |
  expect_equal(res[[1]]$clusters['g2_1', "cluster_id"],
               res[[2]]$clusters['g1_1', "cluster_id"])
  expect_equal(res[[2]]$clusters['g2_1', "cluster_id"],
               res[[2]]$clusters['g2_2', "cluster_id"])
  expect_equal(res[[2]]$clusters['g1_1', "cluster_id"],
               res[[2]]$clusters['g1_2', "cluster_id"])
  expect_equal(res[[3]]$clusters['g2_1', "cluster_id"],
               res[[3]]$clusters['g1_1', "cluster_id"])
  expect_equal(res[[4]]$clusters,
               res[[2]]$clusters)
})

test_that("add_cluster_id works not correct", {
  hits <- Hits(from = c(1L, 2L, 3L), to=c(1L, 2L, 4L),
               nLnode = 5L, nRnode = 5L)
  cluster_df <- data.frame(anchor_id=seq(1L, 5L),
                           cluster_id=seq(1, 5),
                           cluster_size=1)
  res <- add_cluster_id(hits, cluster_df)
  expect_equal(res$cluster_id, c(1, 2, 4))
})

test_that("annoLinker works not correct", {
  addEvidence <- interactive()
  plan(multisession, workers=2)
  anno1 <- annoLinker(peaks, annoData, interactions,
                      addEvidence = addEvidence,
                      verbose=TRUE, parallel=TRUE)
  anno2 <- annoLinker(peaks, annoData, interactions,
                      addEvidence = addEvidence,
                      verbose=FALSE, parallel=TRUE)
  plan(sequential)
  anno3 <- annoLinker(peaks, annoData, interactions,
                      addEvidence = addEvidence,
                      verbose=TRUE, parallel=FALSE)
  anno4 <- annoLinker(peaks, annoData, interactions,
                      addEvidence = addEvidence,
                      verbose=FALSE, parallel=FALSE)
  expect_identical(anno_peaks(anno1), anno_peaks(anno2))
  expect_identical(anno_peaks(anno1), anno_peaks(anno3))
  expect_identical(anno_peaks(anno1), anno_peaks(anno4))

  for(n in c(1, 5, 7, 9)){
    plotEvidence(anno1, event=n, output='trackPlot')
  }
  p <- lapply(c(1, 5), plotEvidence, anno = anno1, output='htmlWidget')
  expect_s3_class(p[[1]], 'htmlwidget')
  expect_s3_class(p[[2]], 'visNetwork')
})


