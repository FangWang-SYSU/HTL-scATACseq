

#' @description
#' This function will creat a peak matrix form peaks list
#' @param peaks description
#' @param object description
#' @param savePath The path to the result output
#' @param file.neme a name of result output
#' @export
AddPeaksSelf <- function(peaks,object,savePath = NULL,file.neme = 'seurObj_addPeaks'){
    frags = object@assays$ATAC@fragments
    counts = FeatureMatrix(fragments = frags, features = peaks)
    chrom_assay <- CreateChromatinAssay(counts = counts,fragments = frags)
    object@assays[["peaks"]] <- chrom_assay
    object@assays$peaks@annotation = object@assays$ATAC@annotation
    if(!is.null(savePath)){
        saveRDS(object,file = file.path(savePath,paste0(file.neme,'.RDS')))
    }
    return(object)
}