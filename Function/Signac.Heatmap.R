library(ComplexHeatmap)
#' @description
#' A short description...
#' @param object description
#' @param assay description
#' @param scMarker description
#' @param plotGene description
#' @param topN description
#' @param groupBy description
#' @param downsample description
#' @param peak.distance description
#' @return description
plotCellHeatmapSelf <- function(object = NULL,
                                assay = NULL,
                                scMarker = NULL,
                                plotGene = F,
                                topN = 5,
                                groupBy = NULL,
                                downsample = NULL,
                                peak.distance = NULL) {
    
    #colnames(sig.peaks) = c("p_val","avg_log2FC" , "pct.1"   ,  "pct.2"   ,  "p_val_adj" ,"cluster" ,  "gene"  )
    #sig.peaks = subset(scMarker,p_val_adj <= p.adjust&avg_log2FC >= log2FC )
    
    if(!is.null(groupBy)){
        Idents(object) = object@meta.data[,groupBy]
    }
    
    
    if(is.null(assay)){
        assay = DefaultAssay(object)
    }else{
        DefaultAssay(object) = assay
    }
    
    if(!is.null(downsample)){
        object = subset(object, downsample = downsample)
    }
    
    
    if(plotGene){
        
        uniquePeaks = unique(scMarker$gene)
        closest_genes <- ClosestFeature(object, regions = uniquePeaks)
        
        # filter
        if(!is.null(peak.distance)){
            closest_genes = subset(closest_genes,distance <= peak.distance)
        }else {
            closest_genes = closest_genes
        }
        
        # set gene name
        freq = table(closest_genes$gene_name)
        freq_name = rep(names(freq),freq)
        freq_time = unlist(lapply(freq,function(x){ seq(1,x,1)}))
        UniqeGene.mat = data.frame(gene_name = freq_name,gene_id = paste0(freq_name,'-',freq_time) )
        
        closest_genes = closest_genes[order(closest_genes$gene_name),]
        closest_genes$gene_id = UniqeGene.mat$gene_id 
        
        scMarker = merge(scMarker,closest_genes,by.x = 'gene',by.y = 'query_region')
        scMarker = scMarker[order(scMarker$cluster),]
        
        #write.table(sig.peaks,file = paste0(outPath,'/differentPeaksAnnotationClosetGene.txt'),quote = F,sep = '\t')
        ## change object
        gene_conne = unique(scMarker[,c('gene','gene_id')])
        new_obj = subset(object,features = gene_conne$gene)
        rownames(new_obj@assays[[assay]]$counts) = gene_conne$gene_id
        
        new_obj2 <- CreateSeuratObject(counts = new_obj@assays[[assay]]$counts,meta.data = new_obj@meta.data)
        Idents(new_obj2) = Idents(new_obj)
        topMP = scMarker %>% group_by(cluster) %>% top_n(n = topN, wt = avg_log2FC) 
        
        
        a = unique(topMP$gene_id)
        # 找到b中的元素在a中的位置
        match_positions <- match(a, rownames(new_obj2))
        # 根据匹配的位置对a进行排序
        sorted_a <- a[order(match_positions)]
        new_obj2 = NormalizeData(new_obj2)
        new_obj2 = ScaleData(new_obj2,features = sorted_a)
        
        if (exists("cb_pattern")) {
            if(length(cb_pattern)> length(unique(Idents(new_obj))) ){
                DoHeatmap(new_obj2, features = topMP$gene_id,size =5,angle = 0,group.colors = cb_pattern) + scale_fill_gradientn(colors = c("#045a8d","white","#a50f15")) 
            }
        }else{
            DoHeatmap(new_obj2, features = topMP$gene_id,size =5,angle = 0) + scale_fill_gradientn(colors = c("#045a8d","white","#a50f15")) 
        }
        
        
    }else{
        topMP = scMarker %>% group_by(cluster) %>% top_n(n = topN, wt = avg_log2FC)
        
        a = unique(topMP$gene)
        # 找到b中的元素在a中的位置
        match_positions <- match(a, rownames(object))
        # 根据匹配的位置对a进行排序
        sorted_a <- a[order(match_positions)]
        
        object = ScaleData(object,features = sorted_a )
        
        if (exists("cb_pattern")) {
            if(length(cb_pattern)> length(unique(Idents(object))) ){
                DoHeatmap(object,assay = assay, features = topMP$gene,size =5,angle = 0,group.colors = cb_pattern) + scale_fill_gradientn(colors = c("#045a8d","white","#a50f15")) 
            }
        }else{
            DoHeatmap(object, assay = assay,slot = 'data', features = topMP$gene,size =5,angle = 0) + scale_fill_gradientn(colors = c("#045a8d","white","#a50f15")) 
        }
        
    }
}

#' @description
#' A short description...
#' @param object seurat object
#' @param assays description
#' @param groupBy description
#' @param topN 默认为NULL，如果不为空，则显示p值最小的marker

 
plotClusterHeatmapSelf <- function(object = NULL,
                                       assays = NULL,
                                       groupBy = NULL,
                                       topN = NULL,
                                       sig.peaks = NULL,
                                       toScale = T,
                                       toLog2Acc = F,
                                       log2Norm = T, 
                                       scaleSize = 10^4,
                                       scaleRows = T,
                                       toLimit = T,
                                       limits = c(-2, 2),
                                       zscore = T,
                                       show.ClosestGene = F,
                                       nLabel = NULL,
                                       PrintMarker = NULL,
                                       pal = NULL,
                                       invert= F,
                                       column_title = NULL,
                                       show_row_names = T,
                                       show_column_names = F,
                                       show_row_dend = T,
                                       show_column_dend = F,
                                       name = 'Zscore'
){
    
    if(is.null(assays)){
        assays = DefaultAssay(object)
    }else{
        DefaultAssay(object) = assays
    }
    
    if(!is.null(groupBy)){
        Idents(object) = object@meta.data[,groupBy]
    }
    
    if(is.null(sig.peaks)){
        stop('请输入计算得到的差异基因矩阵/向量')
    }
    
    
    
    
    if(!is.character(sig.peaks)){
        #colnames(sig.peaks) = c("p_val","avg_log2FC" , "pct.1"   ,  "pct.2"   ,  "p_val_adj" ,"cluster" ,  "gene"  )
        
        if(!is.null(topN)){
            sig.peaks = sig.peaks %>% group_by(cluster)%>%top_n(n = topN,wt = avg_log2FC )
        }
        UniquePeaks = unique(sig.peaks$gene)
    }else{
        UniquePeaks = unique(sig.peaks)
    }
    
   
    if(is.null(PrintMarker)){
        if(!is.null(nLabel)){
            PrintLabel = sig.peaks %>% group_by(cluster)%>%top_n(n = nLabel,wt = avg_log2FC)
            PrintLabel = unique(PrintLabel$gene)
        }
        
    }else{
        PrintLabel = PrintMarker
    }
    
    
    #print(length(PrintLabel))
    
    
    clusters =  unique(Idents(object))
    clusters =  clusters[order(clusters)]
    
    mat =  object@assays[[assays]]@counts
    
    if(assays=='chromvar'){
        mat =  object@assays[[assays]]@data
        
    }
    mat = mat[UniquePeaks,]
 
    #object = subset(object,assays = assays,features = UniquePeaks)
    
    
    meanM = .getClusterMean(object = object,mat = mat, clusters = clusters)
     
    meanM = as.matrix(meanM)
    
    if(toScale){
        meanM = .toscale(mat = meanM,log2Norm = log2Norm,scaleSize = scaleSize,scaleRows = scaleRows,toLimit = toLimit,limits = limits)
    }else if(toLog2Acc) {
        meanM = .log2Acc(mat =meanM)
    }
    
    meanM = as.matrix(t(meanM))
 
    if(zscore){
        meanM = .rowZscores(meanM,limit = T)
    }
    
 
    if(exists('PrintLabel')){
        if(show.ClosestGene){
            closeGene = ClosestFeature(object = object,regions = PrintLabel)
            PrintLabel = closeGene$gene_name
            names(PrintLabel) = closeGene$query_region
            ht1Anno <- HeatmapAnnotation(link = anno_mark(at = which(colnames(meanM) %in%names(PrintLabel)),labels = PrintLabel, labels_gp = gpar(fontsize = 10)))
        }else{
            PrintLabel = intersect(colnames(meanM),PrintLabel)
            ht1Anno <- HeatmapAnnotation(link = anno_mark(at = which(colnames(meanM) %in% PrintLabel),labels = PrintLabel, labels_gp = gpar(fontsize = 10)))
        }
    }else{
        ht1Anno = NULL
    }
   
    #return(meanM)
    plot.oo = .signacHeatmap(mat = meanM,name = name,pal = pal,ht1Anno = ht1Anno,invert= invert,column_title = column_title,show_row_names=show_row_names,show_column_names = show_column_names)
    return(plot.oo)
    
}  

 
.getClusterMean <- function(object = NULL,
                            mat = NULL,
                            clusters = NULL){
    j = 1
    for (i in clusters) {
        #i = clusters[2]
        subBarcode = WhichCells(object = object,idents = i)
        sub_mat = mat[,subBarcode] 
        mean.c = rowMeans(sub_mat)
        if(is.numeric(i)){
            Clu = paste0('C',i)
        }else{
            Clu = i
        }
       
        if (j == 1) {
            eval(parse(text = paste0('mean.mat = data.frame(`',Clu,'`= mean.c)')))
        }else{
            mean.mat[,Clu] = mean.c
        }
        j = j+1
    }
    return(mean.mat)
}


.log2Acc <- function(mat){
    mat = log2(mat+1)
    return(mat)
}

.toscale <- function(mat,
                     log2Norm = T, 
                     scaleSize = 10^4,
                     scaleRows = T,
                     toLimit = T,
                     limits = c(-2, 2)){
    if(log2Norm){
        mat <- log2(t(t(mat)/colSums(mat)) * scaleSize + 1)
    }
    if(scaleRows){
        mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    }
    if(toLimit){
        mat[mat > max(limits)] <- max(limits)
        mat[mat < min(limits)] <- min(limits)
    }
    return(mat)
}


.rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
    z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
    if(limit){
        z[z > max] <- max
        z[z < min] <- min
    }
    return(z)
}


.signacHeatmap <- function(mat,
                           name,
                           pal = NULL,
                           ht1Anno = NULL,
                           cluster_columns = F,
                           cluster_rows = F,
                           show_column_names  = T,
                           show_row_names = F,
                           show_column_dend = F,
                           show_row_dend = T,
                           column_title = NULL,
                           invert= F ){
    
    if(is.null(pal)){
        pal = paletteContinuous(set = "solarExtra", n = 100)
    }
    
    if(invert){
        mat = as.matrix(t(mat))
        plot.H <- Heatmap(mat,name = name,col = pal,show_column_names = show_column_names,show_row_names = show_row_names,cluster_rows = F,right_annotation = ht1Anno, column_title =column_title,
                          cluster_columns = cluster_columns,clustering_method_columns = "ward.D2",show_column_dend =show_column_dend,show_row_dend = show_row_dend)
    }else{
        plot.H <- Heatmap(mat,name = name,col = pal,show_column_names = show_column_names,cluster_rows = cluster_rows,column_title =column_title,
                          #Annotation
                          top_annotation = ht1Anno, 
                          show_row_names = show_row_names,cluster_columns = cluster_columns,clustering_method_columns = "ward.D2",show_column_dend =show_column_dend,show_row_dend = show_row_dend)

    }
 
    return(plot.H)
    
}

#'@description
#'A short description...
#' @param name description
#' @param name description
GetCombineATAC_GeneActivity = function(obj,
                                       assays = 'ATAC',
                                       topN = 20,
                                       marker.peaks,
                                       groupBy= 'CellType',
                                       pal1 = NULL,pal2 = NULL){
    DefaultAssay(obj) = assays;
    
    if (is.null(pal1)) {
        pal1 = paletteContinuous(set = "solarExtra", n = 100)
    }
    if(is.null(pal2)){
        pal2 = paletteContinuous(set = "blueYellow", n = 100)
    }
    
    if(!is.null(groupBy)){
        Idents(obj) = obj@meta.data[[groupBy]] 
    }
    
    clusters =  unique(Idents(obj))
    clusters =  clusters[order(clusters)]
    
    if(!is.null(topN)){
        sig.peaks = marker.peaks %>% group_by(cluster)%>%top_n(n = topN,wt = avg_log2FC)
    }else{
        sig.peaks = subset(marker.peaks,avg_log2FC >=1.5&p_val_adj <=0.01)
    }
    ClosestGene = ClosestFeature(object = obj,regions = unique(sig.peaks$gene))
    ClosestGene = subset(ClosestGene,distance<=20000)
    ClosestGene = ClosestGene[ClosestGene$gene_name%in%rownames(obj@assays$RNA),]
    UniquePeaks = unique(ClosestGene$query_region)
    
    mat =  obj@assays[[assays]]@counts
    mat = mat[UniquePeaks,]
    
    meanM = .getClusterMean(object = obj,mat = mat, clusters = clusters)
    meanM = as.data.frame(meanM)
    meanM[,'region'] = rownames(meanM)
    meanM <- merge(meanM,ClosestGene[,c('gene_name','query_region')],by.x = 'region',by.y = 'query_region')
    meanM = meanM[,-1]
    
    meanM = aggregate(meanM[,1:(ncol(meanM)-1)],list(gene_name = meanM$gene_name),sum)
    rownames(meanM ) =meanM$gene_name
    meanM = meanM[,-1]
    meanM = as.matrix(meanM)
    meanM = .toscale(mat = meanM,log2Norm = T,scaleSize = 10000,scaleRows = T,toLimit = T,limits = c(-2,2))
    
    Peak.heatmap = Heatmap(meanM,  show_column_names = T,cluster_rows = T,name = 'PeakMatrix',col = pal1,
                           show_row_names = T,cluster_columns = F,clustering_method_columns = "ward.D2",show_column_dend =F)
    
    gene_mat = obj@assays$RNA@counts[rownames(meanM),]
    gene_meanM = .getClusterMean(object = obj,mat = gene_mat, clusters = clusters)
    gene_meanM = .toscale(mat = gene_meanM,log2Norm = T,scaleSize = 10000,scaleRows = T,toLimit = T,limits = c(-2,2))
    geneActi.heatmap = Heatmap(gene_meanM,  show_column_names = T,cluster_rows = T,name = 'GeneActivity',col = pal2,
                               show_row_names = T,cluster_columns = F,clustering_method_columns = "ward.D2",show_column_dend =F)
    
    
    draw(Peak.heatmap+geneActi.heatmap)
    
}


#' @description
#' A short description...
#' @param name description
#' @param name description
#' 
CellMarkerEnrich<- function(obj,
                            topMP = topMP,assay = 'peak',
                            ref.geneList =mouse_marker ){
    cluster = unique(topMP$cluster)
    closetGene <- lapply(cluster, function(x){
        #x = cluster[1]
        uniquePeak = topMP[which(topMP$cluster==x),]$gene
        if(assay == 'peak'){
            closest_genes <- ClosestFeature(obj, regions = uniquePeak)
            closest_genes = unique(closest_genes$gene_name) 
        }else{
            closest_genes = unique(uniquePeak)
        }
        
        Enrich <- ref.geneList[ref.geneList$marker%in%closest_genes,]
        if(dim(Enrich)[1]>0){
            Enrich$cluster  = x
        }
        return(Enrich)
    })
    marker.dat <- do.call(rbind,closetGene) 
    return(marker.dat)
}
