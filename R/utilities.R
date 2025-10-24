# ============================================================================
# VALIDATION
# ============================================================================

#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom S4Vectors elementNROWS
#' @importFrom methods as is
validate_inputs_graph <- function(peaks, annoData, interactions,
                                  bindingRegion, interactionDistanceRange) {
  if (!inherits(peaks, "GRanges")) {
    stop("'peaks' must be a GRanges object", call. = FALSE)
  }

  if (!inherits(annoData, c("GRanges"))) {
    stop("'annoData' must be an annoGR or GRanges object", call. = FALSE)
  }

  if (length(bindingRegion) != 2) {
    stop("'bindingRegion' must have exactly 2 elements", call. = FALSE)
  }

  if (length(interactionDistanceRange) != 2) {
    stop("'interactionDistanceRange' must have exactly 2 elements", call. = FALSE)
  }

  if (is.null(annoData$gene_id) && is.null(annoData$symbol) &&
    is.null(annoData$gene_name)) {
    warning("annoData lacks gene identifiers (gene_id, symbol, gene_name)")
  }

  if (!inherits(interactions, c("GRanges", "GInteractions", "Pairs"))) {
    stop("interactions must be an GInteractions object", call. = FALSE)
  }

  if (is(interactions, "GRanges")) {
    if (length(interactions$blocks) == 0) {
      stop("If interactions is GRanges, the blocks metadata is required", call. = FALSE)
    }
    if (any(elementNROWS(interactions$blocks) != 2)) {
      stop("If interactions is GRanges, the length of each metadata blocks must be 2.", call. = FALSE)
    }
  }

  if (length(intersect(seqlevelsStyle(peaks), seqlevelsStyle(annoData))) < 1) {
    stop("Please check the seqlevels style of your peaks and annoData", call. = FALSE)
  }

  if (length(intersect(seqlevelsStyle(peaks), seqlevelsStyle(first(interactions)))) < 1) {
    stop("Please check the seqlevels style of your peaks and interactions", call. = FALSE)
  }
}

# ============================================================================
# FILTER INTERACTIONS
# ============================================================================
#' @importFrom S4Vectors first second
#' @importFrom IRanges IRangesList shift
#' @importFrom BiocGenerics start
#' @importFrom GenomicRanges GRanges
#' @importFrom Seqinfo seqnames
#' @importFrom InteractionSet GInteractions pairdist
filterInteractions <- function(interactions, interactionDistanceRange) {
  if (!is(interactions, "GInteractions")) {
    if (is(interactions, "Pairs")) {
      ## is Pairs
      interactions <- GInteractions(
        first(interactions),
        second(interactions),
        score=mcols(interactions)$score
      )
    }
    if (is(interactions, "GRanges")) {
      HiC_FIRST <- lapply(interactions$blocks, `[`, 1)
      HiC_SECOND <- lapply(
        interactions$blocks, `[`,
        2
      )
      HiC_FIRST <- unlist(IRangesList(HiC_FIRST))
      HiC_SECOND <- unlist(IRangesList(HiC_SECOND))
      HiC_FIRST <- shift(HiC_FIRST, start(interactions) -
        1)
      HiC_SECOND <- shift(HiC_SECOND, start(interactions) -
        1)
      interactions <-
        GInteractions(
          GRanges(seqnames(interactions), HiC_FIRST),
          GRanges(seqnames(interactions), HiC_SECOND),
          score = mcols(interactions)$score
        )
    }
  }
  d <- pairdist(interactions)
  gi_filt <- interactions[!is.na(d) &
    d >= interactionDistanceRange[1] &
    d <= interactionDistanceRange[2]]
  ## clean up
  GInteractions(first(gi_filt), second(gi_filt))
}

# ============================================================================
# FIND ANNOTATION REGIONS
# ============================================================================
#' @importFrom GenomicRanges promoters terminators
get_annotation_regions <- function(annoData, bindingType, bindingRegion) {
  switch(bindingType,
    nearestBiDirectionalPromoters = {
      suppressWarnings(promoters(annoData,
        upstream = abs(bindingRegion[1]),
        downstream = bindingRegion[2]
      ))
    },
    startSite = {
      suppressWarnings(promoters(annoData,
        upstream = abs(bindingRegion[1]),
        downstream = bindingRegion[2]
      ))
    },
    endSite = {
      suppressWarnings(terminators(annoData,
        upstream = abs(bindingRegion[1]),
        downstream = bindingRegion[2]
      ))
    },
    stop("Unsupported binding type: ", bindingType, call. = FALSE)
  )
}


# ============================================================================
# BUILD INTERACTION GRAPH
# ============================================================================
#' @importFrom igraph graph_from_data_frame simplify
build_interaction_graph <- function(anchors, weight, cluster_method, ...) {
  # Create edge list: genes connected through same interaction
  edge_list <- do.call(cbind, anchors)

  if (nrow(edge_list) == 0) {
    return(list(graph = NULL, clusters = NULL))
  }

  if(!missing(weight)){
    edge_list <- cbind(edge_list, weight=weight)
  }

  # Build graph
  g <- graph_from_data_frame(edge_list, directed = FALSE)

  # Simplify: remove self-loops and multiple edges
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  # Detect clusters
  clusters <- detect_clusters(g, cluster_method, ...)

  return(list(
    graph = g,
    clusters = clusters
  ))
}

# ============================================================================
# DETECT CLUSTERS
# ============================================================================
#' @importFrom igraph components cluster_louvain cluster_walktrap cluster_infomap membership V
detect_clusters <- function(g, cluster_method, ...) {
  clusters_obj <- switch(cluster_method,
    components = {
      components(g)
    },
    louvain = {
      cluster_louvain(g, ...)
    },
    walktrap = {
      cluster_walktrap(g, ...)
    },
    infomap = {
      cluster_infomap(g, ...)
    },
    stop("Unknown cluster method: ",
      cluster_method,
      call. = FALSE
    )
  )

  # Extract membership
  if (cluster_method == "components") {
    membership <- clusters_obj$membership
  } else {
    membership <- membership(clusters_obj)
  }

  # Create cluster data frame
  cluster_df <- data.frame(
    anchor_id = V(g)$name,
    cluster_id = membership,
    stringsAsFactors = FALSE
  )

  # Add cluster size
  cluster_sizes <- table(membership)
  cluster_df$cluster_size <- cluster_sizes[as.character(cluster_df$cluster_id)]

  return(cluster_df)
}

# ============================================================================
# ADD CLUSTER ID TO HITS
# ============================================================================

add_cluster_id <- function(hits, cluster_df) {
  stopifnot(is(hits, "Hits"))
  stopifnot(is(cluster_df, "data.frame"))
  hits <- as.data.frame(hits)
  hits$cluster_id <- cluster_df$cluster_id[match(
    hits$subjectHits,
    cluster_df$anchor_id
  )]
  return(hits)
}

# ============================================================================
# SHORTEST PATH
# ============================================================================
#' @importFrom future.apply future_lapply
#' @importFrom igraph shortest_paths as_edgelist
#' @importFrom progressr with_progress progressor
find_shortest_path <- function(peak_ol_anno, interaction_graph,
                               parallel=FALSE, verbose=FALSE) {
  peak_ol_anno_subset <- unique(peak_ol_anno[, c(
    "subjectHits.peak",
    "subjectHits.annotation"
  )])
  ## do not try shortest_paths in batch, unpredictable

  # choose the apply function based on `parallel`
  apply_fun <- if (parallel) future_lapply else lapply

  # define the core worker function
  get_shortest_path_edges <- function(i, p = NULL) {
    sp_edges <- as_edgelist(interaction_graph$graph)[
      shortest_paths(
        interaction_graph$graph,
        from = as.character(peak_ol_anno_subset$subjectHits.peak[[i]]),
        to   = as.character(peak_ol_anno_subset$subjectHits.annotation[[i]]),
        mode = "all",
        output = "epath"
      )$epath[[1]],
      ,
      drop = FALSE
    ]

    if (!is.null(p)) p()  # update progress if verbose

    sp_edges
  }

  # wrap everything depending on verbosity
  n_pairs <- nrow(peak_ol_anno_subset)
  args <- list(X=seq_len(n_pairs), FUN=get_shortest_path_edges)
  if(parallel){
    args$future.seed <- TRUE
  }
  sp <- if(verbose){
    on.exit(message('Please do not forget to run future::plan(future::sequential) to release the resources'))
    with_progress({
      args$p <- progressor(steps = n_pairs)
      do.call(apply_fun, args)
    })
  }else{
    do.call(apply_fun, args)
  }

  names(sp) <- apply(peak_ol_anno_subset, 1, paste, collapse = ",")
  sp
}

# ============================================================================
# ANNOTATE PEAKS WITH CLUSTERS
# ============================================================================
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom BiocGenerics start end strand
annotate_peaks_with_clusters <- function(
  peaks,
  peak_ol_anno,
  evidences,
  annoData,
  interRegion
) {
  annoted_peaks <- peaks[peak_ol_anno$queryHits.peak]
  annotations <- annoData[peak_ol_anno$queryHits.annotation]
  annoted_peaks$feature_name <- names(annotations)
  annoted_peaks$feature_start <- start(annotations)
  annoted_peaks$feature_end <- end(annotations)
  annoted_peaks$feature_strand <- strand(annotations)
  mcols(annoted_peaks) <- cbind(mcols(annoted_peaks), mcols(annotations))
  annoted_peaks$peak_bin <-
    as.character(interRegion[peak_ol_anno$subjectHits.peak])
  annoted_peaks$feature_bin <-
    as.character(interRegion[peak_ol_anno$subjectHits.annotation])
  keep <- !duplicated(annoted_peaks) |
    !duplicated(mcols(annoted_peaks))
  annoted_peaks <- annoted_peaks[keep]
  peak_ol_anno <- peak_ol_anno[keep, , drop=FALSE]
  if(length(evidences)>0){
    annoted_peaks$evidences <- vapply(evidences, function(.ele) {
      if (nrow(.ele)) {
        paste(paste(interRegion[as.numeric(.ele[, 1])],
                    interRegion[as.numeric(.ele[, 2])],
                    sep = " | "
        ), collapse = "; ")
      } else {
        ""
      }
    }, FUN.VALUE = character(1L))[
      apply(
        peak_ol_anno[, c(
          "subjectHits.peak",
          "subjectHits.annotation"
        )], 1,
        paste,
        collapse = ","
      )
    ]
  }

  return(annoted_peaks)
}
