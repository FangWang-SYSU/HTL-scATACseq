
library(cowplot)
# source package
library(EnsDb.Hsapiens.v86)
library(ggsci)
library(Signac)
library(Seurat)
library(ggplot2)
library(ggsci)
library(ggthemes)
library(VennDiagram)
library(cowplot)
library(ggpubr)
library(ArchR)
library(IRanges)
library(GenomicRanges)
library(viridis)


#' @description
#' read .bed file
#' @param file the file path 
#' @param header description
readBED <-  function(file,header = F){
    mat = read.delim(file = file,header = header)
    gr = GRanges(seqnames = mat$V1,ranges = IRanges(start = mat$V2,end = mat$V3) )
    index = dim(mat)[2]
    mcols(gr) = mat[,4:index]
    return(gr)
}




#' @describeIn 
#' @param peaks a GRanges of peak list
#' @param Ref.gr a GRanges of known peak list 
#' @param Ref.Name
VennPlot.cCRE.Overlaps.self <- function(peaks,Ref.gr,PeakName = 'Myself',Ref.Name = 'Encode',plot.title = 'Overlaps'){
    peaks =    GenomicRanges::reduce(peaks)
    Ref.gr =  GenomicRanges::reduce(Ref.gr)
    inter_Peaks = subsetByOverlaps(peaks,Ref.gr)
    
    name1 = paste0(Ref.Name,'(',length(Ref.gr),')')
    name2 = paste0(PeakName,'(',length(peaks),')')
    
    plot.aa = draw.pairwise.venn(area1 = length(Ref.gr),    
                                 area2 = length(peaks),
                                 cross.area = length(inter_Peaks),
                                 category = c(name1,name2),
                                 lwd = 1,
                                 cat.cex = 1,  
                                 cat.pos = c(-50,130),
                                 cat.dist = 0,
                                 cat.just = list(c(1,1), c(0, 0)),
                                 col= c("#440154ff", '#21908dff'),
                                 fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                                 cat.col = c("#440154ff", '#21908dff'),
                                 cat.default.pos = "outer",
                                 cex = 0.8, 
                                 fontfamily = "sans",
                                 filename = NULL,
                                 ext.line.lty = "dashed",
                                 output=F,
                                 compression = "lzw" )
    
    plot.aa = ggplotify::as.ggplot(cowplot::as_grob(plot.aa)) +
        labs(title = plot.title) +
        theme( plot.title = element_text(hjust = 0.5),plot.margin = unit(c(1.5,0,1.5,0),"cm"))+
        xlim(-0.4, 1.4)
    grid.newpage()
    dev.off()
    return(plot.aa)
}







#' 
#' #' @description
#' #' A short description...
#' #' @param name description
#' #' @param
#' 
#' PlotPeaksCompare.self <- function(Ref.PeakList,PeaksGroup.list,ncol = 3){
#'     PrintPlot = list()
#'     
#'     for(a in names(Ref.PeakList)){
#'         
#'         Ref_gr = Ref.PeakList[[a]]
#'         
#'         for(b in names(PeaksGroup.list)){
#'             nAme = paste0(a,'_',b)
#'             #
#'             Peaks.List_N = PeaksGroup.list[[b]]
#'             
#'             plot.list = list()
#'             for (c in names(Peaks.List_N)) {
#'                 peaks = Peaks.List_N[[c]]
#'                 plot.list[[c]] = VennPlot.cCRE.Overlaps.self(peaks = peaks,Ref.gr = Ref_gr,Ref.Name = a, plot.title = c)
#'             }
#'             
#'             num_plot = length(names(plot.list))
#'             
#'             if(num_plot >= 10){
#'                 ncol = 4
#'             }else{
#'                 ncol = 3
#'             }
#'             
#'            
#'             plot_grid_args <- c(plot.list[c(1:num_plot)],ncol = ncol)
#'             PrintPlot[[nAme]] = do.call(plot_grid, plot_grid_args)
#' 
#'         }
#'         
#'         
#'     }
#'     return(PrintPlot)
#'     
#' }


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


