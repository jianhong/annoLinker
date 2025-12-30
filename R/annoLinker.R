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
#' @param bindingType Character, one of "startSite", "body", or "endSite"
#' @param bindingRegion Numeric vector of length 2 defining promoter window
#'   (e.g., c(-5000, 5000))
#' @param cluster_method Character, clustering method: "components" (connected
#'   components), "louvain", "walktrap", or "infomap"
#' @param maxgap Integer, bp to extend interaction anchors for
#'   overlap detection (default: 0)
#' @param interactionDistanceRange Numeric vector of length 2 defining the
#'   minimal and maximal distance of interactions. This is used to make sure
#'   the annotations are not supper far away.
#' @param addEvidence Logical, add evidence to the metadata or not.
#' @param parallel Logical, use future_lapply to do parallel computing or not.
#' @param verbose Logical, print the message or not
#' @param ... Parameters for cluster. see \link[igraph]{cluster_louvain},
#' \link[igraph]{cluster_walktrap}, and \link[igraph]{cluster_infomap}.
#' @return An \link{annoLinkerResult} object
#' or NULL if no annotations found
#'
#' @export
#' @importFrom igraph vcount ecount V<-
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom InteractionSet anchorIds regions
#' @importFrom methods new
#' @examples
#' ## read the peaks and interactions
#' library(rtracklayer)
#' extPath <- system.file("extdata", package = "annoLinker")
#' peaks <- rtracklayer::import(file.path(extPath, "peaks.bed"))
#' interactions <- rtracklayer::import(file.path(extPath, "interaction.bedpe"))
#' library(TxDb.Drerio.UCSC.danRer10.refGene)
#' annoData <- genes(TxDb.Drerio.UCSC.danRer10.refGene)
#' anno <- annoLinker(peaks, annoData, interactions, verbose = TRUE)
annoLinker <- function(peaks, annoData, interactions,
    bindingType = c("startSite", "body", "endSite"),
    bindingRegion = c(-5000, 5000),
    cluster_method = c("components", "louvain", "walktrap", "infomap"),
    maxgap = 0, interactionDistanceRange = c(10000, 10000000),
    addEvidence = FALSE, parallel = FALSE, verbose = FALSE, ...) {
    bindingType <- match.arg(bindingType)
    cluster_method <- match.arg(cluster_method)
    validate_inputs_graph(peaks, annoData, interactions,
        bindingRegion, interactionDistanceRange)
    totalSteps <- ifelse(addEvidence, 5, 4)
    verboseMessage(verbose, '1', totalSteps=totalSteps)
    # Filter interactions
    interactions <- filterInteractions(interactions, interactionDistanceRange)
    # Build interaction network using igraph
    interaction_graph <- graph_by_interaction(interactions, cluster_method, ...)
    if (is.null(interaction_graph$graph)) {
        return(msgNULL("Failed to build interaction graph")) }
    verboseMessage(verbose, '1.1', g = interaction_graph)
    verboseMessage(verbose, '2', totalSteps=totalSteps)
    interRegion <- regions(interactions)
    # Find peaks that overlap with any anchor region
    if (inherits(annoData, "annoGR"))  annoData <- as(annoData, "GRanges")
    annoRegion <- get_annotation_regions(annoData, bindingType, bindingRegion)
    annoOL <- findOverlaps(annoRegion, interRegion, maxgap = maxgap)
    if (length(annoOL) == 0) {
        return(msgNULL("No annotation found overlapping interaction anchors"))}
    peakOL <- findOverlaps(peaks, interRegion, maxgap = maxgap)
    if (length(peakOL) == 0) {
        return(msgNULL("No peak found overlapping interaction anchors"))}
    verboseMessage(verbose, '3', totalSteps=totalSteps)
    # merge genes hits and peaks hits by the cluster id
    peak_ol_anno <- get_peak_ol_anno(annoOL, peakOL, interaction_graph)
    if (nrow(peak_ol_anno) == 0) {
        return(msgNULL("No peaks & annotation found in same cluster"))}
    evidences <- NULL
    if (addEvidence) {
        verboseMessage(verbose, '4', totalSteps=totalSteps)
        evidences <- find_shortest_path(
            peak_ol_anno, interaction_graph, parallel, verbose) }
    verboseMessage(verbose, '5', totalSteps=totalSteps)
    annotated_peaks <- annotate_peaks_with_clusters(
        peaks, peak_ol_anno, evidences, annoData, interRegion)
    verboseMessage(verbose, '6', num=length(annotated_peaks))
    # rename interaction graph from region id to region
    g <- fix_g_name(interaction_graph, interRegion)
    interaction_graph <- fix_insteraction_name(interaction_graph, interRegion)
    return(new("annoLinkerResult", annotated_peaks = annotated_peaks,
            graph = g, clusters = interaction_graph$clusters))
}
