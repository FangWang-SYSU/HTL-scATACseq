#' @description 
#' @param fragPath
#' @param outPath
#' @param Species
CreateArchProject <- function(fragPath =NULL,
                              outPath,
                              Species ='mm10'){
    library(ArchR)
    library(ggplot2)
    library(ggsci)
    library(Seurat)
    addArchRThreads(threads = 8) 
    if (Species == 'mm10') {
       addArchRGenome("mm10")
    }
    sampleN = list.files(fragPath)
    inputFiles = paste0(fragPath,sampleN,'/outs/fragments.tsv.gz')
    names(inputFiles) = sampleN
    ## 创建Arrow文件
    ArrowFiles <- createArrowFiles(inputFiles = inputFiles,sampleNames = names(inputFiles),
                                   minTSS = 2, #Dont set this too high because you can always increase later
                                   minFrags = 1000, 
                                   addTileMat = TRUE,addGeneScoreMat = TRUE,verbose = F)
    ## 创建ArchRProject
    projHeme1 <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "ArchRProject",copyArrows = TRUE)
    saveArchRProject(ArchRProj = projHeme1, outputDirectory = file.path(outPath,"Save-ProjHeme1"), load = FALSE)
    return(projHeme1)
}

#' @description 
#' @param object
#' @param projHeme
ChangeArchRUmap <- function(object,
                            projHeme){
    UMAP.Info = as.data.frame(object@reductions[["umap"]]@cell.embeddings)
    UMAP.Info$CellBarcode = object$ArchRCellBarcode
    rownames(UMAP.Info) = UMAP.Info$CellBarcode
    UMAP= UMAP[projHeme2$cellNames,]
    UMAP = UMAP[,1:2]
    colnames(UMAP)= c('IterativeLSI#UMAP_Dimension_1','IterativeLSI#UMAP_Dimension_2')
    projHeme@embeddings@listData[["UMAP"]]@listData[["df"]] = UMAP
    return(projHeme)
}


#' @description 
#' @param projHeme
#' @param groupBy
AddMatrix_ArchR <- function(projHeme =NULL,
                            groupBy = 'Clusters',pathToMacs2=findMacs2()){
            projHeme <- addGroupCoverages(ArchRProj = projHeme, groupBy = groupBy)
            
            projHeme <- addReproduciblePeakSet(ArchRProj = projHeme, groupBy = groupBy, pathToMacs2 = pathToMacs2)
            projHeme <- addPeakMatrix(projHeme)
            # Add Motif
            projHeme <- addMotifAnnotations(ArchRProj = projHeme, motifSet = "cisbp", name = "Motif")
            projHeme <- addBgdPeaks(projHeme)
            projHeme <- addDeviationsMatrix(ArchRProj = projHeme, peakAnnotation = "Motif",force = TRUE)
            return(projHeme)
}
