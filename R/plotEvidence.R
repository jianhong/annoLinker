#' Plot interaction network for visualization
#'
#' @param g Interaction graph output by \link{annoLinker}
#' @param evidence The evidence to highlighted
#' @param cluster_method Character, clustering method: "components" (connected
#'   components), "louvain", "walktrap", or "infomap"
#' @param ... parameters for cluster. see \link[igraph]{cluster_louvain},
#'   \link[igraph]{cluster_walktrap}, and \link[igraph]{cluster_infomap}.
#' @param output Output of the plot.
#' @param txdb,org The TxDb and OrgDb object used for annotation plot.
#' @export
#' @importFrom visNetwork toVisNetworkData visOptions visNetwork
#' @importFrom igraph induced_subgraph as_edgelist
#' @examples
#' anno <- readRDS(system.file('extdata', 'sample_res.rds', package='annoLinker'))
#'
#'
plotEvidence <- function(
  g, evidence,
  cluster_method = c("components", "louvain", "walktrap", "infomap"),
  ...,
  output = c("graph", "htmlWidget", "trackPlot"),
  txdb, org
) {
  stopifnot(is.character(evidence))
  stopifnot(length(evidence) == 1)
  stopifnot(is(g, "igraph"))
  cluster_method <- match.arg(cluster_method)
  output <- match.arg(output)
  if (output == "track") {
    if (!missing(txdb) && !missing(org)) {
      stopifnot(is(txdb, "TxDb"))
      stopifnot(is(org, "OrgDb"))
    }
  }

  # Detect clusters
  clusters <- detect_clusters(g, cluster_method, ...)

  # extract the cluster id
  evi_names <- unique(strsplit(evidence, ";|\\|"))
  evi_names <- gsub(pattern = "\\s+", replacement = "", evi_names[[1]])
  cluster_id <- clusters$cluster_id[clusters$anchor_id %in% evi_names]
  if (!all(cluster_id == cluster_id[1])) {
    stop("Not all cluster id is same. Check the cluster_method and parameters.")
  }
  # extract all the components in the same cluster
  css <- clusters$anchor_id[clusters$cluster_id == cluster_id[1]]
  sg <- induced_subgraph(g, css)
  if (output == "track") {
    return(plotTrack(as_edgelist(sg, names = TRUE), evi_names, txdb, org))
  }
  data <- toVisNetworkData(sg)
  data$nodes$title <- data$nodes$id
  data$nodes$color <- ifelse(data$nodes$id %in% evi_names, "tomato", "lightgray")
  data$nodes$size <- ifelse(data$nodes$id %in% evi_names, 20, 15)
  graph <- visNetwork(data$nodes, data$edges)
  switch(output,
    "htmlWidget" = {
      visOptions(graph, highlightNearest = TRUE, nodesIdSelection = TRUE)
    },
    "graph" = graph
  )
}

#' @importFrom trackViewer geneTrack trackList viewTracks setTrackStyleParam gi2track
#' @importFrom InteractionSet GInteractions regions
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom IRanges subsetByOverlaps
#' @importFrom AnnotationDbi select
plotTrack <- function(edges, evi_names, txdb, org) {
  evi_names <- matrix(evi_names, ncol = 2, byrow = TRUE)
  evi <- GInteractions(GRanges(evi_names[, 1]), GRanges(evi_names[, 2]))
  evi$score <- 2
  gi <- GInteractions(GRanges(edges[, 1]), GRanges(edges[, 2]))
  gi$score <- 1
  gi <- unique(c(evi, gi))
  gi <- gi[order(gi$score)]
  track <- gi2track(gi)
  setTrackStyleParam(track, "tracktype", "link")
  setTrackStyleParam(
    track, "color",
    c("white", "gray", "red")
  )
  range <- range(regions(gi))
  if (!missing(txdb) && !missing(org)) {
    genes <- suppressMessages(genes(txdb))
    genes <- subsetByOverlaps(genes, range)
    symbols <- suppressMessages(select(org, names(genes), "SYMBOL", "ENTREZID"))
    anno <- geneTrack(symbols$ENTREZID, txdb, symbols$SYMBOL,
      type = "gene",
      asList = FALSE
    )
    highlight <- queryHits(findOverlaps(anno$dat, regions(evi)))
    anno$dat$color <- "black"
    anno$dat$color[highlight] <- "red"
    tl <- trackList(genes = anno, links = track)
  } else {
    tl <- trackList(links = track)
  }
  viewTracks(tl, gr = range, autoOptimizeStyle = TRUE)
}
