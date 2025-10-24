#' Plot interaction network for visualization
#'
#' @param anno An object of annoLinkerResult output by \link{annoLinker}
#' @param event Number to indicate the event to be plot
#' @param output Output of the plot.
#' @param colors Colors settting for the plot.
#' @param txdb,org The TxDb and OrgDb object used for annotation plot.
#' @export
#' @importFrom visNetwork toVisNetworkData visOptions visNetwork
#' @importFrom igraph induced_subgraph as_edgelist shortest_paths add_vertices add_edges
#' @examples
#' anno <- readRDS(system.file('extdata', 'sample_res.rds', package='annoLinker'))
#' library(org.Dr.eg.db)
#' library(TxDb.Drerio.UCSC.danRer10.refGene)
#' n <- 1 #length(anno$annotated_peaks$evidences)
#' plotEvidence(anno, event=n,
#'  output='htmlWidget')
#' plotEvidence(anno, event=n,
#'  output='trackPlot')
plotEvidence <- function(
    anno, event,
    output = c("graph", "htmlWidget", "trackPlot"),
    colors = c(peak='darkgreen', feature='brown',
               node='tomato', background='lightgray'),
    txdb, org
) {
  stopifnot(is(anno, 'annoLinkerResult'))
  stopifnot(is.numeric(event))
  stopifnot(length(event) == 1)
  if(event>length(anno)){
    stop("event parameter is greater than available data.")
  }
  if(event<0){
    stop("event number could not smaller than 1")
  }
  stopifnot(all(c('peak', 'feature', 'node', 'background') %in% names(colors)))
  event <- round(event)
  output <- match.arg(output)
  if (output == "trackPlot") {
    if (!missing(txdb) && !missing(org)) {
      stopifnot(is(txdb, "TxDb"))
      stopifnot(is(org, "OrgDb"))
    }
  }

  # prepare the data
  clusters <- anno_clusters(anno)
  g <- anno_graph(anno)
  peakRegion <- anno_event(anno, event)
  fetureRegion <- anno_feature(anno, event)
  evidence <- anno_evidence(anno, event)
  peakbin <- anno_peakbin(anno, event)
  featurebin <- anno_featurebin(anno, event)
  if(evidence==''){
    ## find evidence
    evidence <- shortest_paths(
      graph = g,
      from = peakbin,
      to   = featurebin,
      mode = "all",
      output = "epath"
    )$epath[[1]]
    evi_names <- as_edgelist(g)[evidence, , drop=FALSE]
  }else{
    # extract the cluster id
    evi_names <- strsplit(evidence, ";|\\|")
    evi_names <- gsub(pattern = "\\s+", replacement = "", evi_names[[1]])
  }
  if(length(evi_names)==0){
    stop("No indirect evidence is available. Maybe they are connect directly.")
  }

  cluster_id <- clusters$cluster_id[clusters$anchor_id %in% evi_names]
  if (!all(cluster_id == cluster_id[1])) {
    stop("Not all cluster id is same. Check the cluster_method and parameters.")
  }
  # extract all the components in the same cluster
  css <- clusters$anchor_id[clusters$cluster_id == cluster_id[1]]
  sg <- induced_subgraph(g, css)
  if (output == "trackPlot") {
    vp <- plotTrack(as_edgelist(sg, names = TRUE), evi_names, txdb, org,
                    peakRegion, fetureRegion, colors)
    return(vp)
  }
  ## add peakRegion and fetureRegion to graph
  A <- as.character(peakRegion)
  B <- as.character(fetureRegion)
  sg <- add_vertices(sg, 2, name=c(A, B))
  sg <- add_edges(sg, c(A, peakbin, B, featurebin))
  data <- toVisNetworkData(sg)
  data$nodes$title <- data$nodes$id
  data$nodes$label[data$nodes$id %in% A] <- 'peak'
  data$nodes$label[data$nodes$id %in% B] <- 'feature'

  data$nodes$color <- ifelse(data$nodes$id %in% evi_names,
                             colors['node'],
                             ifelse(data$nodes$id %in% A,
                                    colors['peak'],
                                    ifelse(data$nodes$id %in% B,
                                           colors['feature'],
                                           colors['background'])))
  data$nodes$size <- ifelse(data$nodes$id %in% evi_names,
                            20,
                            ifelse(data$nodes$id %in% c(A, B),
                                   25, 15))
  graph <- visNetwork(data$nodes, data$edges)
  switch(output,
    "htmlWidget" = {
      visOptions(graph, highlightNearest = TRUE, nodesIdSelection = TRUE)
    },
    "graph" = graph
  )
}

#' @importFrom trackViewer geneTrack trackList viewTracks setTrackStyleParam gi2track addGuideLine setTrackYaxisParam addArrowMark
#' @importFrom InteractionSet GInteractions regions
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom IRanges subsetByOverlaps
#' @importFrom AnnotationDbi select
plotTrack <- function(edges, evi_names, txdb, org, peakRegion, fetureRegion, colors) {
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
    c("white", colors['background'], colors['node'])
  )
  setTrackYaxisParam(track, 'draw', FALSE)
  range <- range(regions(gi))
  if (!missing(txdb) && !missing(org)) {
    genes <- suppressMessages(genes(txdb))
    genes <- subsetByOverlaps(genes, range)
    symbols <- suppressMessages(select(org, names(genes), "SYMBOL", "ENTREZID"))
    symbols <- symbols[!is.na(symbols$SYMBOL), ]
    anno <- geneTrack(symbols$ENTREZID, txdb, symbols$SYMBOL,
      type = "gene",
      asList = FALSE
    )
    highlight <- queryHits(findOverlaps(anno$dat, regions(evi)))
    anno$dat$color <- colors['background']
    anno$dat$color[highlight] <- colors['node']
    tl <- trackList(genes = anno, links = track)
  } else {
    tl <- trackList(links = track)
  }
  vp <- viewTracks(tl, gr = range, autoOptimizeStyle = TRUE)
  if(!missing(peakRegion)){
    addGuideLine(guideLine=c(start(peakRegion), end(peakRegion)),
                 col=colors['peak'], vp=vp)
    addGuideLine(guideLine=c(start(fetureRegion), end(fetureRegion)),
                 col=colors['feature'], vp=vp)
    addArrowMark(pos = list(x=(start(peakRegion) + end(peakRegion))/2,
                            y=length(tl)),
                 label = 'peak', col = colors['peak'], vp = vp)
    addArrowMark(pos = list(x=(start(fetureRegion) + end(fetureRegion))/2,
                            y=length(tl)),
                 label = 'feature', col = colors['feature'], vp = vp)
  }
  return(invisible())
}
