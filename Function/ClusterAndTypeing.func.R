suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(cowplot))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(patchwork))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(ggthemes))
suppressMessages(library(ArchR))
suppressMessages(library(pheatmap))
 




cb_pattern = c("#4dbbd5ff" ,"#E64b35ff", "#00a087ff" ,"#a65628", "#FF95A8FF","#BCBD22FF", "#fdbf6f", "#3c5488ff",
               "#f39b7fff", "#b09c85ff", "#7876b1ff", "#377eb8",  "#4daf4a","#97A1A7FF" ,
               "#984ea3" , "#ff7f00",  "#f781bf", "#b2df8a", "#5050FFFF", "#82581FFF" , "#E5614CFF",
               "#F0E685FF", "#D595A7FF", "#CDDEB7FF","#612A79FF" ,"#AE1F63FF", 
               "#99CC00FF","#CC9900FF" ,"#9467BDFF", "#EFD500FF" , "#ADE2D0FF",pal_igv()(50))


# define lable fun
give.n = function(x) {
    la.df = data.frame(y = mean(x)+ max(x) / 10,label = round(mean(x), 1));
    return(la.df)
}

#' #' @description
#' #' A short description...
#' #' @param name description
#' VlnPlot.self = function(obj,
#'                         fontsize = 15,
#'                         ncol = 5,
#'                         linesize = 0.35,
#'                         group.by = NULL,
#'                         feature = c('nCount_ATAC', 'nFeature_ATAC','log10_nCount',"log10_nFrags",
#'                                     'TSS.enrichment', 'nucleosome_signal',
#'                                     'pct_reads_in_peaks', 'peak_region_fragments',   'blacklist_ratio')){
#'     
#'     fea1 = colnames(obj@meta.data)
#'     feature = intersect(fea1,feature)
#'     
#'     meta = as.data.table(obj@meta.data)
#'     df = meta[, feature, with = FALSE]
#'     
#'     # plot num
#'     num_plot = length(feature)
#'     
#'     # assign colors
#'     cols = pal_igv()(num_plot)
#'     col.ls = setNames(cols,feature)
#'     column.M = colMeans(df)
#'     column.ls = setNames(column.M,feature)
#'     if(is.null(group.by)){
#'         df[,'orig.ident'] = meta[,'orig.ident']
#'     }else{
#'         df[,'orig.ident'] = meta[,group.by, with = FALSE]
#'     }
#'     
#'     # plot
#'     gp.ls = df[,feature, with = FALSE] %>% imap( ~ {
#'         ggplot(data = df, aes(x = orig.ident, y = .x)) +
#'             geom_violin(trim = FALSE, fill = col.ls[.y]) +
#'             ggtitle(label = .y) + ylab(label = .y) +
#'             theme_classic() +
#'             theme(
#'                 panel.grid.major = element_blank(),
#'                 panel.grid.minor = element_blank(),
#'                 #strip.background = element_blank(),
#'                 panel.border = element_blank()
#'             ) +
#'             theme(
#'                 axis.text.x = element_text(size = fontsize,angle = 60,hjust = 1),
#'                 axis.title.x = element_blank(),
#'                 axis.ticks.length = unit(.05, "cm"),
#'                 plot.title = element_text(size = fontsize + 4, hjust = 0.5),
#'                 legend.position = 'none'
#'             ) + 
#'             stat_summary(fun = mean, geom = "point", col = "black") +  # Add points to plot
#'             stat_summary(fun.data = give.n,size = fontsize-5,position = 'identity',
#'                          geom = "text",
#'                          col = "black")
#'         
#'     })
#'     plot_grid_args <- c(gp.ls[c(1:num_plot)],ncol = ncol)
#'     do.call(plot_grid, plot_grid_args)
#'     
#' }



#' @description getGeneActivity
#' @param obj description
#' 
getGeneActivity <- function(obj){
    gene.activities <- GeneActivity(obj)
    # add the gene activity matrix to the Seurat object as a new assay and normalize it
    obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
    obj <- NormalizeData(
        object = obj,
        assay = 'RNA',
        normalization.method = 'LogNormalize',
        scale.factor = median(obj$nCount_RNA))
    return(obj)
}


 
#'@description
#' give.mat.RNA
#' @param obj description
#' @param marker description
give.mat.RNA <-  function(obj,marker,Clu){
    Gene = marker$Gene
    Gene = intersect(Gene,rownames(obj@assays$RNA))
    dat = as.matrix(obj@assays$RNA@data[Gene,])
    dat = as.data.frame(t(dat))
    dat[,'Clusters'] = obj@meta.data[,Clu]
    data2 <- dat %>% group_by(Clusters) %>% summarise_all(list(mean)) %>% as.data.frame()
    dat_long2 <- melt(data2,id.vars =c('Clusters'),variable.name = 'Gene' ,value.name='Mean')
    dat_long2$IDD = paste0(dat_long2$Clusters,'_',dat_long2$Gene)
    
    
    dat[,'CellId'] = colnames(obj)
    dat_long <- melt(dat,id.vars =c('Clusters','CellId'),variable.name = 'Gene' ,value.name='Expr')
    head(dat_long)
    dat_long$IDD = paste0(dat_long$Clusters,'_',dat_long$Gene)
    
 
    plot.mat = merge(dat_long,dat_long2[,c('IDD','Mean')],by = 'IDD')
    plot.mat = merge(plot.mat,marker,by.x = 'Gene',by.y = 'Gene')
    return(plot.mat)
}


#'@description
#' cellType.Plot.RNA.violin
#' @param obj description
#' @param Clu description
#' @param Marker description
cellType.Plot.RNA.violin <- function(obj,Clu = 'seurat_clusters', Marker = NULL){
    Marker$Cluster = factor(Marker$Cluster,levels = unique(Marker$Cluster))
    mat = give.mat.RNA(obj = obj,marker = Marker,Clu =Clu)
    P1 = ggplot(mat,aes(x = Clusters,y = Expr,fill = Mean)) + 
        geom_violin(scale = "width")+
        theme_few()+
        scale_fill_gradientn( colours = c("#2166ac","#92c5de","#ffeda0" ,"#F8B195","#F67280","#C06C84","#6C5B7B")) +
        facet_grid( Cluster+Gene ~. ,scales = 'free_y')+
        theme(strip.background=element_rect(colour="black", fill="white"),axis.text.x = element_text(angle = 60, hjust = 1, vjust =1),
              
              strip.text = element_text(size = 7),strip.text.y = element_text(angle = 0) ) 
    return(P1)
}

#'@description
#' cellType.Plot.RNA.dot
#' @param obj description
#' @param Clu description
#' @param Marker description
cellType.Plot.RNA.dot <- function(obj,Clu = 'seurat_clusters',Marker = NULL,invert = F){
    Gene = unique(Marker$Gene)
    Gene = as.character(Gene)[length(Gene):1]
    if(invert){
        
        P2= DotPlot(object = obj,features = Gene,assay = 'RNA') +
            scale_color_gradientn(colours = c("#053061","#2166ac","#4393c3","#d1e5f0","#fddbc7",
                                              "#d6604d","#b2182b","#67001f"))+
            theme(axis.text.x=element_text(colour="black",angle = 90,hjust = 1,vjust = 0.5))
    }else{
        P2= DotPlot(object = obj,features = Gene,assay = 'RNA') +
            coord_flip()+
            scale_color_gradientn(colours = c("#053061","#2166ac","#4393c3","#d1e5f0","#fddbc7",
                                              "#d6604d","#b2182b","#67001f"))+
            theme(axis.text.x=element_text(colour="black",angle = 90,hjust = 1,vjust = 0.5))
    }
    
   
    return(P2)
}



#' #' @description
#' #' A short description...
#' #' @param obj description
#' #' @param features description
#' #' @param ncol description
#' VlnPlot.self2 = function(obj,
#'                          features = NULL,
#'                          ncol = 1){
#'     give.m = function(x) {
#'         la.df = data.frame(y = mean(x) + max(x)/5,label = paste0(round(mean(x), 2)));
#'         return(la.df)
#'     }
#'     P = lapply(features,FUN = function(x){
#'         VlnPlot(obj,features =x,cols = cb_pattern,pt.size = 0,raster=FALSE) +
#'             stat_summary(fun = mean, geom = "point", col = "black") +
#'             stat_summary(fun.data = give.m,size = 4,position = 'identity',
#'                          geom = "text",
#'                          col = "black")+NoLegend()
#'     })
#'     plot_grid_args <- c(P[c(1:length(features))],ncol = ncol)
#'     do.call(plot_grid, plot_grid_args)
#' }



#' @description
#' Marker.Gene.Heatmap
#' @param obj seurat obj
#' @param top 选择前几个基因进行分析
Marker.Gene.Heatmap <- function(obj,top = 5,markerGene = markerGene,assay = NULL,
                                downsample = NULL,plot_peaks = F){
    
   
    top5MG = markerGene %>%dplyr::filter(avg_log2FC > 0.5)%>% group_by(cluster) %>% top_n(n = -top, wt = p_val_adj)
    new_obj = ScaleData(obj,assay = assay, features = top5MG$gene)
    if(!is.null(downsample)){
        new_obj = subset(new_obj, downsample = downsample)
    }
    p1 = DoHeatmap(new_obj,assay = assay, features = top5MG$gene,size =5) + scale_fill_gradientn(colors = c("#045a8d","white","#a50f15")) 
    
    if(plot_peaks){
        library(tidyr)
        top5MG_Gene = ClosestFeature(obj,regions = top5MG$gene)
        gene_df <- top5MG_Gene[,c("gene_name","closest_region")]
        gene_df$positi <- c(dim(gene_df)[1]:1)
        gene_df <- gather(gene_df,key = 'geneType', value = 'ddd',-positi)
        p2 = ggplot()+
            geom_text(data = gene_df,aes(x = geneType, y=positi,label = ddd)) +
            theme_void()
        p1+p2
    }else{
        p1
    }
    
}



library(viridis)
library(ggpubr)


#' @description
#' A short description...
#' 
Motif.Enrichment <- function(obj_addMotif,da_peaks,assay = NULL,Cluster.by = 'seurat_clusters'){
    #da_peaks = subset(da_peaks,p_val_adj<0.01)
    open.peaks <- AccessiblePeaks(obj_addMotif)
    meta.feature <- GetAssayData(obj_addMotif, assay = assay, layer = "meta.features")
    
    celltype = unique(obj_addMotif@meta.data[,Cluster.by])
    celltype = as.character(celltype)
    celltype = intersect(celltype,unique(da_peaks$cluster))
    motif.mat = data.frame()
    
    for(i in celltype){
        diff.peaks = subset(da_peaks,cluster == i)
        peaks.name = diff.peaks[,'gene']
        peaks.matched <- MatchRegionStats(
            meta.feature = meta.feature[open.peaks, ],
            query.feature = meta.feature[peaks.name, ],
            n = 50000
        )
        
        # find motifs
        enriched.motif = FindMotifs(object = obj_addMotif,features = peaks.name,background = peaks.matched,verbose = F)
        enriched.motif[,'Cluster'] = i
        motif.mat = rbind(motif.mat,enriched.motif)
    }
    return(motif.mat)
}
    
pheatmap.motif <- function(enriched.motifs,n = 10,fc = 2,p.cut = 10,MotifRef = NULL){
    motif.mat = enriched.motifs[,c('motif','fold.enrichment','p.adjust','Cluster', 'motif.name')]
    motif.mat$`-log10P` <- -log10(motif.mat$p.adjust)
    motif.mat[which(motif.mat$`-log10P`  > 100),'-log10P'] <- 100
    sig.motifMat = subset(motif.mat,`-log10P` >=p.cut&fold.enrichment>fc)
    length(unique(sig.motifMat$motif.name))
    
    top.motifMat <- sig.motifMat %>% group_by(motif.name) %>% slice_max(n=1,order_by = `-log10P`,with_ties = F)
    
    motif = top.motifMat%>%group_by(Cluster)%>%slice_max(n=n ,order_by = `-log10P`,with_ties = F)
    motifName = unique(motif$motif.name)
    if(is.null(MotifRef)){
        motifName = motifName
    }else{
        motifName = intersect(motifName,MotifRef)
    }
    
    motif.mat = motif.mat[motif.mat$motif.name%in%motifName,]
    
    # long to width 
    motif.mat.w = dcast(data = motif.mat,motif+motif.name~Cluster,value.var = '-log10P',fill = 0)
    row.names(motif.mat.w) = motif.mat.w$motif.name
 
    
    motif.mat.w = motif.mat.w[,3:ncol(motif.mat.w)]
   
    # pheatmap
    mat = as.matrix(t(motif.mat.w))

    return(mat)
}

    
