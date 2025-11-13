suppressMessages(library(EnsDb.Mmusculus.v79))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(patchwork))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(data.table))
suppressMessages(library(cowplot))
suppressMessages(library(purrr))
suppressMessages(library(ggthemes))
suppressMessages(library(ArchR))

cb_pattern = c("#4dbbd5ff" ,"#E64b35ff", "#00a087ff" ,"#a65628", "#FF95A8FF","#BCBD22FF", "#fdbf6f", "#3c5488ff",
               "#f39b7fff", "#b09c85ff", "#7876b1ff", "#377eb8",  "#4daf4a","#97A1A7FF" ,
               "#984ea3" , "#ff7f00",  "#f781bf", "#b2df8a", "#5050FFFF", "#82581FFF" , "#E5614CFF",
               "#F0E685FF", "#D595A7FF", "#CDDEB7FF","#612A79FF" ,"#AE1F63FF","#1B1919FF",
               "#99CC00FF","#CC9900FF" ,"#9467BDFF", "#EFD500FF" , "#ADE2D0FF",pal_igv()(50))


#'  @description
#' CreateSeuratObject.scATAC
#' @param dataPath cellrange的结果的文件夹里边含有样本结果文件夹和相应.mro文件
#' @param sampleName description
#' @param min.cells  default min.cells = 10
#' @param min.features default min.features = 200
#' @param Peaks default Peaks = 'macs2'; peaks count 数据产生的方式，1.通过macs2方法进行call peaks 获得，2.10x 产生的h5文件中自带；3.同时导入两张数据

CreateSeuratObject.scATAC <- function(dataPath,Species = 'hg38',
                                      sampleName,
                                      min.cells = 10,
                                      min.features = 200,
                                      Peaks = 'cellranger'){
    
    #' 读取每个样本detected cellular barcodes的单细胞数据（raw count数据）
    #' # peak注释
    if(Species=='hg38'){
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
        seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
        genome(annotations) <- 'GRCh38'
    }else if(Species=='mm10'){
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
        seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
        genome(annotations) <- "mm10"
        
    }
    
    #读取重要的三个文件(peak_bc_matrix、singlecell.csv、fragments.tsv.gz)
    
    if (file.exists(file.path(dataPath,"filtered_peak_bc_matrix.h5"))) {
        peak_path = file.path(dataPath,"filtered_peak_bc_matrix.h5")
        meta_path=file.path(dataPath,"singlecell.csv")
        fragment_path = file.path(dataPath,"fragments.tsv.gz");
        
        #读取peak-count data
        counts <- Read10X_h5(filename = peak_path)
        #建立样本名字meta
        metadata <- read.csv(file = meta_path,header = TRUE, row.names = 1)
        barcodes = rownames(subset(metadata,is__cell_barcode=='1'))
        
    }else{
        # 文件地址
        peak_path <-  file.path(dataPath,"matrix.mtx")
        barcodes <- file.path(dataPath,"barcodes.tsv")
        peaknames_Path <- file.path(dataPath,'peaks.bed')
        fragment_path = file.path(dataPath,"fragments.tsv.gz");
        
        # 读取文件
        counts <- Matrix::readMM(peak_path)
        barcodes <- readLines(barcodes)
        peaknames  <- read.table(peaknames_Path, sep="\t")
        peaknames <- paste(peaknames$V1, peaknames$V2, peaknames$V3, sep="-")
        colnames(counts) <- barcodes
        rownames(counts) <- peaknames
        
        metadata = data.frame(row.names = barcodes,Barcode = barcodes)
        head(metadata)
    }
    
    # peak calling
    if(Peaks == 'macs2'){
        frags <- CreateFragmentObject(path = fragment_path, cells = barcodes,verbose = F)
        peaks <- CallPeaks(frags,verbose = F)
        
        # Quantify fragments in each peak
        counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)
        chrom_assay <- CreateChromatinAssay(counts = counts, 
                                            fragments = fragment_path,
                                            annotation = annotations,
                                            sep = c(":", "-"),
                                            min.cells = min.cells,
                                            min.features = min.features,
                                            verbose = F)
        #统一使用seurat对象
        temp.seurat <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", meta.data =  metadata[colnames(chrom_assay),])
        temp.seurat$sampleID = sampleName
        
    }else if (Peaks == 'cellranger'){
        # chrom_assay
        chrom_assay <- CreateChromatinAssay(counts = counts,
                                            fragments = fragment_path,
                                            annotation = annotations,
                                            sep = c(":", "-"),
                                            min.cells = min.cells,
                                            min.features = min.features,
                                            verbose = F)
        #统一使用seurat对象
        temp.seurat <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", meta.data = metadata[colnames(chrom_assay),])
        temp.seurat$sampleID = sampleName
    }else if(Peaks == 'both'){
        #读取peak-count data
        counts <- Read10X_h5(filename = peak_path)
        #建立样本名字meta
        metadata <- read.csv(file = meta_path,header = TRUE, row.names = 1)
        # chrom_assay
        chrom_assay <- CreateChromatinAssay(counts = counts,
                                            sep = c(":", "-"),
                                            fragments = fragment_path,
                                            min.cells = min.cells,
                                            annotation = annotations,
                                            min.features = min.features,
                                            verbose = F)
        
        #统一使用seurat对象
        temp.seurat <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", meta.data =  metadata[colnames(chrom_assay),])
        temp.seurat$sampleID = sampleName
        peaks = CallPeaks(temp.seurat)
        peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
        peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
        # quantify counts in each peak
        macs2_counts <- FeatureMatrix(fragments = Fragments(temp.seurat), features = peaks, cells = colnames(temp.seurat))
        temp.seurat[["Peaks"]] <- CreateChromatinAssay(counts = macs2_counts, fragments = fragment_path, annotation = annotations )
    }else{
        warning('Please choose a peak count method')
    }
    return(temp.seurat)
}

#' @description 
#' CalculateQC.index
#' @param obj seurat对象
#'
CalculateQCindexSelf <- function(obj){
    obj <- NucleosomeSignal(object = obj,verbose = F) #fragment ratio 147-294: <147  ---  mononucleosome:nucleosome-free
    obj <- TSSEnrichment(object = obj,fast = FALSE,verbose = F)
    obj$log10_nFeature_ATAC <- log10(obj$nFeature_ATAC)
    obj$log10_nCount_ATAC <- log10(obj$nCount_ATAC)
    obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
    obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments
    if(any(names(obj@assays)%in%'peaks')){
        obj$log10_nFeature_peaks <- log10(obj$nFeature_peaks)
        obj$log10_nCount_peaks <- log10(obj$nCount_peaks)
    }
    return(obj)
}



#' @description
#' QC.index
#' @param obj seurat对象
#' @param feature 需要输出的特征
#'
returnQCindex <- function(obj,feature){
    fea = colnames(obj@meta.data)
    fea = intersect(fea,feature)
    dat = obj@meta.data[,fea]
    dat2 = as.data.frame(t(apply(dat, 2, quantile)))
    dat2[,'Mean'] = colMeans(dat)
    return(round(dat2,2))
}


# define lable fun
give.n = function(x) {
    la.df = data.frame(y = mean(x)+ max(x) / 10,label = round(mean(x), 2));
    return(la.df)
}

#' @description VlnPlotSelf
#' A short description...
#' @param name description
VlnPlotSelf = function(obj,
                        fontsize = 15,
                        ncol = 5,
                        linesize = 0.35,
                        group.by = NULL,
                        feature = NULL){
    
    fea1 = colnames(obj@meta.data)
    feature = intersect(fea1,feature)
    
    meta = as.data.table(obj@meta.data)
    df = meta[, feature, with = FALSE]
    
    # plot num
    num_plot = length(feature)
    
    # assign colors
    cols = pal_igv()(num_plot)
    col.ls = setNames(cols,feature)
    column.M = colMeans(df)
    column.ls = setNames(column.M,feature)
    if(is.null(group.by)){
        df[,'orig.ident'] = meta[,'orig.ident']
    }else{
        df[,'orig.ident'] = meta[,group.by, with = FALSE]
    }
    
    # plot
    gp.ls = df[,feature, with = FALSE] %>% imap( ~ {
        ggplot(data = df, aes(x = orig.ident, y = .x)) +
            geom_violin(trim = FALSE, fill = col.ls[.y]) +
            ggtitle(label = .y) + ylab(label = .y) +
            theme_classic() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #strip.background = element_blank(),
                panel.border = element_blank()
            ) +
            theme(
                axis.text.x = element_text(size = fontsize,angle = 60,hjust = 1),
                axis.title.x = element_blank(),
                axis.ticks.length = unit(.05, "cm"),
                plot.title = element_text(size = fontsize + 4, hjust = 0.5),
                legend.position = 'none'
            ) + 
            stat_summary(fun = mean, geom = "point", col = "black") +  # Add points to plot
            stat_summary(fun.data = give.n,size = fontsize-5,position = 'identity',
                         geom = "text",
                         col = "black")
        
    })
    plot_grid_args <- c(gp.ls[c(1:num_plot)],ncol = ncol)
    do.call(plot_grid, plot_grid_args)
    
}


#' @description
#' A short description...
#' @param obj description
#' @param features description

QC.of.each.Cluster= function(obj,features= c('nFeature_ATAC','nCount_ATAC','TSS.enrichment')){
    dat = obj@meta.data[,c('seurat_clusters',features)]
    aa = aggregate(dat[,2:ncol(dat)], by=list(Cluster=dat$seurat_clusters),mean)
    aa$Cluster = paste0('C',aa$Cluster)
    aa[,2:ncol(dat)] = round(aa[,2:ncol(dat)],2)
    All  = round(colMeans(dat[,2:ncol(dat)]),2)
    cat('Mean of all cells:\n')    
    for(i in features){
        cat(i,' = ',All[i],'\n')
    }
    cat('--------------------------------------------------------------------\n')
    for(i in features){
        cat('Cluster :', aa[which(aa[i]>All[i]),1],'> mean', i, '\n')
    }
    return(t(aa))
}

#' @description
#' A short description...
#' @param obj description
#' @param features description
#' @param ncol description
VlnPlotSelf2 = function(obj,
                         features = NULL,
                         ncol = 1){
    give.m = function(x) {
        la.df = data.frame(y = median(x) + max(x)/8,label = paste0(round(median(x), 2)));
        return(la.df)
    }
    P = lapply(features,FUN = function(x){
        VlnPlot(obj,features =x,cols = cb_pattern,pt.size = 0,raster=FALSE,) +
            stat_summary(fun = median, geom = "point", col = "black") +
            stat_summary(fun.data = give.m,size = 4,position = 'identity',
                         geom = "text",
                         col = "black")+NoLegend()
    })
    plot_grid_args <- c(P[c(1:length(features))],ncol = ncol)
    do.call(plot_grid, plot_grid_args)
}



