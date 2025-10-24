#' Class \code{"annoLinkerResult"}
#'
#' An object of class \code{"annoLinkerResult"} represent the annotated peaks,
#'  which is a GRanges object with peaks annotated by gene clusters, and
#' interaction graph, which is an igraph graph.
#'
#' @name annoLinkerResult-class
#' @aliases annoLinkerResult
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("annoLinkerResult", annotated_peaks, graph, clusters)}.
#' @keywords classes
#' @export
#' @examples
#' library(igraph)
#' library(GenomicRanges)
#' new("annoLinkerResult",
#'  annotated_peaks=GRanges(),
#'  graph=make_empty_graph(),
#'  clusters=data.frame())
setClass("annoLinkerResult",
         representation = representation(
           annotated_peaks="GRanges",
           graph="ANY",
           clusters = "data.frame"
         ),
         validity = function(object){
           re <- NULL
           mc <- mcols(object@annotated_peaks)
           if(nrow(mc)){
             if(!all(c('feature_name',
                       'feature_start',
                       'feature_end',
                       'feature_strand',
                       'peak_bin',
                       'feature_bin') %in% colnames(mc))){
               re <- c(re, "One of more of columns feature_name, feature_start, feature_end, feature_strand peak_bin, feature_bin, is missing for the annotated_peaks")
             }
           }
           if(nrow(object@clusters)){
             if(!all(c('anchor_id', 'cluster_id') %in% colnames(object@clusters))){
               re <- c(re, "'anchor_id', 'cluster_id' is missing for the cluster slot.")
             }
           }
           if(!inherits(object@graph, "igraph")){
             re <- c(re, "graph must be an object of 'igraph'")
           }
           re
         })

#' @name coerce
#' @rdname annoLinkerResult-class
#' @aliases coerce,annoLinkerResult,GRanges-method
#' @exportMethod coerce
#' @importFrom methods coerce
setAs(from="annoLinkerResult", to="GRanges", function(from){
  from@annotated_peaks
})

#' @rdname annoLinkerResult-class
#' @param x,object An annoLinkerResult object.
#' @param row.names,optional,... parameters used by \link[base]{as.data.frame}
#' @exportMethod as.data.frame
#' @aliases as.data.frame,annoLinkerResult-method
setMethod("as.data.frame", signature(x="annoLinkerResult"), function(x, row.names = NULL, optional = FALSE, ...){
  as.data.frame(x@annotated_peaks, ...)
})

#' @export
#' @rdname annoLinkerResult-class
setGeneric('anno_peaks', function(x) standardGeneric('anno_peaks'))
#' @rdname annoLinkerResult-class
#' @exportMethod anno_peaks
#' @aliases anno_peaks,annoLinkerResult-method
setMethod("anno_peaks", signature(x="annoLinkerResult"), function(x){
  x@annotated_peaks
})

#' @export
#' @rdname annoLinkerResult-class
setGeneric('anno_graph', function(x) standardGeneric('anno_graph'))
#' @rdname annoLinkerResult-class
#' @exportMethod anno_graph
#' @aliases anno_graph,annoLinkerResult-method
setMethod("anno_graph", signature(x="annoLinkerResult"), function(x){
  x@graph
})

#' @export
#' @rdname annoLinkerResult-class
setGeneric('anno_clusters', function(x) standardGeneric('anno_clusters'))
#' @rdname annoLinkerResult-class
#' @exportMethod anno_clusters
#' @aliases anno_clusters,annoLinkerResult-method
setMethod("anno_clusters", signature(x="annoLinkerResult"), function(x){
  x@clusters
})

#' @export
#' @rdname annoLinkerResult-class
setGeneric('anno_evidence', function(x, i) standardGeneric('anno_evidence'))
#' @rdname annoLinkerResult-class
#' @exportMethod anno_evidence
#' @param i Numeric, index value.
#' @aliases anno_evidence,annoLinkerResult,numeric-method
setMethod("anno_evidence", signature(x="annoLinkerResult"), function(x, i){
  if(length(x@annotated_peaks$evidences)){
    x@annotated_peaks$evidences[i]
  }else{
    ''
  }
})

#' @export
#' @rdname annoLinkerResult-class
setGeneric('anno_event', function(x, i) standardGeneric('anno_event'))
#' @rdname annoLinkerResult-class
#' @exportMethod anno_event
#' @aliases anno_event,annoLinkerResult,numeric-method
setMethod("anno_event", signature(x="annoLinkerResult"), function(x, i){
  peaks <- x@annotated_peaks
  if(length(peaks)){
    mcols(peaks) <- NULL
    peaks[i]
  }else{
    GRanges()
  }
})

#' @export
#' @rdname annoLinkerResult-class
setGeneric('anno_feature', function(x, i) standardGeneric('anno_feature'))
#' @rdname annoLinkerResult-class
#' @exportMethod anno_feature
#' @aliases anno_feature,annoLinkerResult,numeric-method
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
setMethod("anno_feature", signature(x="annoLinkerResult"), function(x, i){
  peaks <- x@annotated_peaks
  if(length(peaks)){
    GRanges(as.character(seqnames(peaks[i])),
            IRanges(peaks$feature_start[i],
                    peaks$feature_end[i]),
            strand=peaks$feature_strand[i])
  }else{
    GRanges()
  }
})

#' @export
#' @rdname annoLinkerResult-class
setGeneric('anno_peakbin', function(x, i) standardGeneric('anno_peakbin'))
#' @rdname annoLinkerResult-class
#' @exportMethod anno_peakbin
#' @aliases anno_peakbin,annoLinkerResult,numeric-method
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
setMethod("anno_peakbin", signature(x="annoLinkerResult"), function(x, i){
  peaks <- x@annotated_peaks
  if(length(peaks)){
    peaks$peak_bin[i]
  }else{
    ''
  }
})

#' @export
#' @rdname annoLinkerResult-class
setGeneric('anno_featurebin', function(x, i) standardGeneric('anno_featurebin'))
#' @rdname annoLinkerResult-class
#' @exportMethod anno_featurebin
#' @aliases anno_featurebin,annoLinkerResult,numeric-method
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
setMethod("anno_featurebin", signature(x="annoLinkerResult"), function(x, i){
  peaks <- x@annotated_peaks
  if(length(peaks)){
    peaks$feature_bin[i]
  }else{
    ''
  }
})

#' @rdname annoLinkerResult-class
#' @exportMethod length
#' @aliases length,annoLinkerResult-method
setMethod("length", signature(x="annoLinkerResult"), function(x){
  length(x@annotated_peaks)
})

#' @rdname annoLinkerResult-class
#' @exportMethod show
#' @importFrom methods show
#' @aliases show,annoLinkerResult-method
setMethod("show", signature(object="annoLinkerResult"), function(object){
  cat("An object of annoLinkerResult with annotated_peaks and interaction graph.\n")
  show(object@annotated_peaks)
  show(object@graph)
})
