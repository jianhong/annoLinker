#' Annotate Peaks with DNA Interaction Networks Using Graph Clustering
#'
#' @description
#' Fast annotation of genomic peaks using DNA interaction data by building
#' interaction networks with igraph. Peaks overlapping any node in a connected
#' subgraph are annotated with all genes in that subgraph.
#'
#' @param peaks GRanges object containing peak regions
#' @param annoData annoGR or GRanges object with gene annotations
#' @param interactions GInteractions, or Pairs object with
#'   interaction data (e.g., Hi-C, ChIA-PET)
#' @param bindingType Character, one of "startSite", or "endSite"
#' @param bindingRegion Numeric vector of length 2 defining promoter window
#'   (e.g., c(-5000, 5000))
#' @param cluster_method Character, clustering method: "components" (connected
#'   components), "louvain", "walktrap", or "infomap"
#' @param extend_anchors Integer, bp to extend interaction anchors for
#'   overlap detection (default: 0)
#' @param interactionDistanceRange Numeric vector of length 2 defining the minimal
#'   and maximal distance of interactions. This is used to make sure the annotations
#'   are not supper far away.
#' @param parallel Logical, use future_lapply to do parallel computing or not.
#' @param verborse Logical, print the message or not
#' @param ... Parameters for cluster. see \link[igraph]{cluster_louvain},
#' \link[igraph]{cluster_walktrap}, and \link[igraph]{cluster_infomap}.
#' @return GRanges object with peaks annotated by gene clusters, or empty
#'   GRanges if no annotations found
#'
#' @export
#' @importFrom igraph vcount ecount V<-
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom InteractionSet anchorIds regions
#' @examples
#' ## read the peaks and interactions
#' library(rtracklayer)
#' extPath <- system.file('extdata', package='annoLinker')
#' peaks <- rtracklayer::import(file.path(extPath, 'peaks.bed'))
#' interactions <- rtracklayer::import(file.path(extPath, 'interaction.bedpe'))
#' library(TxDb.Drerio.UCSC.danRer10.refGene)
#' annoData <- genes(TxDb.Drerio.UCSC.danRer10.refGene)
#' anno <- annoLinker(peaks, annoData, interactions, verborse=TRUE)
annoLinker <- function(
  peaks,
  annoData,
  interactions,
  bindingType = c("startSite", "endSite"),
  bindingRegion = c(-5000, 5000),
  cluster_method = c("components", "louvain", "walktrap", "infomap"),
  extend_anchors = 0,
  interactionDistanceRange = c(10000, 10000000),
  parallel = FALSE,
  verborse = FALSE,
  ...
) {
  # Validate inputs
  bindingType <- match.arg(bindingType)
  cluster_method <- match.arg(cluster_method)
  validate_inputs_graph(peaks, annoData, interactions, bindingRegion, interactionDistanceRange)
  if (verborse) {
    message("Step 1/5: Building interaction network graph...")
  }
  # Filter interactions
  interactions <- filterInteractions(interactions, interactionDistanceRange)
  # Extract interaction anchor regions
  anchors <- anchorIds(interactions)
  # Build interaction network using igraph
  interaction_graph <- build_interaction_graph(
    anchors,
    cluster_method,
    ...
  )

  if (is.null(interaction_graph$graph)) {
    message("Failed to build interaction graph")
    return(GRanges())
  }

  if (verborse) {
    message(sprintf(
      "Graph has %d nodes and %d edges in %d clusters, on average %f nodes in each cluster",
      vcount(interaction_graph$graph),
      ecount(interaction_graph$graph),
      length(unique(interaction_graph$clusters$cluster_id)),
      vcount(interaction_graph$graph) /
        length(unique(interaction_graph$clusters$cluster_id))
    ))
  }
  if (verborse) {
    message("Step 2/5: Finding interaction anchors overlapping genes and peaks...")
  }
  interRegion <- regions(interactions)
  # Find peaks that overlap with any anchor region
  if (inherits(annoData, "annoGR")) {
    annoData <- as(annoData, "GRanges")
  }
  annoRegion <- get_annotation_regions(annoData, bindingType, bindingRegion)
  annoOL <- findOverlaps(annoRegion, interRegion, maxgap = extend_anchors)
  if (length(annoOL) == 0) {
    message("No annotation found overlapping interaction anchors")
    return(GRanges())
  }
  peakOL <- findOverlaps(peaks, interRegion, maxgap = extend_anchors)
  if (length(peakOL) == 0) {
    message("No peak found overlapping interaction anchors")
    return(GRanges())
  }

  if (verborse) {
    message("Step 3/5: Finding genes and peaks in the same cluster...")
  }
  # convert interRegion hit id to cluster id
  annoOL <- add_cluster_id(annoOL, interaction_graph$clusters)
  peakOL <- add_cluster_id(peakOL, interaction_graph$clusters)
  # merge genes hits and peaks hits by the cluster id
  peak_ol_anno <- merge(peakOL, annoOL,
    by = "cluster_id", all = FALSE,
    suffixes = c(".peak", ".annotation")
  )

  if (nrow(peak_ol_anno) == 0) {
    message("No peaks and annotation found in same interaction cluster")
    return(GRanges())
  }

  if (verborse) {
    message("Step 4/5: Finding evidence of the annotation...")
  }
  evidences <- find_shortest_path(peak_ol_anno, interaction_graph, parallel)

  if (verborse) {
    message("Step 5/5: Annotating peaks with gene clusters...")
  }

  annotated_peaks <- annotate_peaks_with_clusters(
    peaks,
    peak_ol_anno,
    evidences,
    annoData,
    interRegion
  )

  if (verborse) {
    message(sprintf(
      "Annotated %d peaks with gene clusters",
      length(annotated_peaks)
    ))
  }
  # rename interaction graph from region id to region
  g <- interaction_graph$graph
  V(g)$name <- as.character(interRegion[as.numeric(V(g)$name)])

  return(
    list(
      annotated_peaks = annotated_peaks,
      interaction_graph = g
    )
  )
}
