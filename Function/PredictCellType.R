library(ggsci)
Predict_cellType_by_Enrichment <- function(FindAllMarkers_res, ref_list, topn=20, scale=TRUE){
    toplist <- FindAllMarkers_res %>% group_by(cluster) %>% top_n(topn,wt = avg_log2FC)
    cell_types = names(ref_list)
    cluster = sort(unique(toplist$cluster))
    inter_mat = matrix(0, ncol = length(cluster), nrow = length(cell_types))
    
    for (i in 1:length(cell_types)) {
        for (j in 1:length(cluster)) {
            pm1 = ref_list[[i]]
            pm2 = toplist[toplist$cluster==cluster[j], "gene_name", drop=TRUE]
            inter_mat[i,j] = length(intersect(pm1, pm2))
        }
    }
    rownames(inter_mat) = c(cell_types)
    inter_mat <- inter_mat[!rowSums(inter_mat)==0,]
    
    if (is.null(rownames(inter_mat))) {
        print("没有富集到任何细胞类型")
        return(NULL)
    }else{
        norma = sapply(rownames(inter_mat), function(x) length(ref_list[[x]])) / rowSums(inter_mat)
        enri = apply(inter_mat,2,function(x)x/norma)
        enri = apply(enri, 2, function(x) x/sum(x))
        
        if(scale){
            enri = t(scale(t(enri)))
        }
        rownames(enri) = rownames(inter_mat)
        colnames(enri) = cluster
        enri[which(is.nan(enri))] = 0
        #print(enri)
        if(sum(enri)==0){
            print("没有富集到任何细胞类型")
            return(NULL)
        }else{
            pheatmap::pheatmap(enri,border_color = 'white')
            return(enri)
        }
        
    }
    
}

#' @description
#' A short description...
#' @param predict_data description

Extract.predict.cell.type <- function(predict_data){
    predict_data[which(predict_data==0)] <- NA
    m <- list()
    n <- list()
    cluster <- colnames(predict_data)
    ref.celltype <- rownames(predict_data)
    for(ce in cluster) {
        #ce = cluster[7]
        m[[ce]] <- ref.celltype[which.max(x =predict_data[,ce] )]
        if (length(m[[ce]])== 0) {
            m[[ce]] = paste0('unknow_',ce)
            #cat(ce,m[[ce]],'\n')
        }else{
            m[[ce]] <- substr(m[[ce]],1,nchar(m[[ce]]))
            #cat(ce,m[[ce]],'\n')
        }
    }
    mm <- unlist(m)
    return(mm)
}