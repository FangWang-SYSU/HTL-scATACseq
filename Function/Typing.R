library(ggplot2)
library(Seurat)
library(Signac)
library(ggsci)
library(ggthemes)
library(ArchR)
library(ggpubr)
library(dplyr)

cb_pattern = c("#4dbbd5ff" ,"#E64b35ff", "#00a087ff" ,"#a65628", "#FF95A8FF","#BCBD22FF", "#fdbf6f", "#3c5488ff",
               "#f39b7fff", "#b09c85ff", "#7876b1ff", "#377eb8",  "#4daf4a","#97A1A7FF" ,
               "#984ea3" , "#ff7f00",  "#f781bf", "#b2df8a", "#5050FFFF", "#82581FFF" , "#E5614CFF",
               "#F0E685FF", "#D595A7FF", "#CDDEB7FF","#612A79FF" ,"#AE1F63FF","#1B1919FF",
               "#99CC00FF","#CC9900FF" ,"#9467BDFF", "#EFD500FF" , "#ADE2D0FF",pal_igv()(50))


 

.percemtHeatmap <- function(object,index1,index2,invert = F,name = NULL ){
    cM <- confusionMatrix(object@meta.data[,index1], object@meta.data[,index2])
    cM <- (cM / Matrix::rowSums(cM))*100
    if(invert){
        mat = as.matrix(t(cM))
    }else{
        mat = as.matrix(cM)
    }
    Pl = Heatmap(mat,
                 col =viridis(100),
                 # width = unit(18, "cm"),
                 # height = unit(8, "cm"),
                 
                 cell_fun = function(j, i, x, y, width, height, fill){
                     grid.text(sprintf("%.1f", mat[i, j]), 
                               x, y, gp = gpar(fontsize =6))},
                 rect_gp = gpar(col = "#525252",lwd = 0.8),show_row_dend = F,show_column_dend = F, 
                 name = name )
    return(Pl)
}


plot.percent <- function(object = NULL,type1,type2,name1,name2,title){
    p1 = .percemtHeatmap(object,index1 = type1,index2 = type2,name = name1)
    p2 = .percemtHeatmap(object,index1 = type2,index2 = type1,invert = T,name = name2 )
    draw(p1+p2,column_title = title)
    
}




######
# @author        : xinwang 
# @version       : 1.0.0
# @desc          : observe_data / expected_data (Roe)
# @param obj     : Seurat object of T cell. 
# @param group_by1     : column name in seurat, such as 'seurat_clusters'.
# @param group_by2     : column name in seurat, such as 'tissue'.
# @param show.sign     : heatmap display string, (sign, value, FALSE)
# @return        : A list contains (Roe matrix, pvalue, display_string). 
######

cell_distribution <- function(obj, group_by1, group_by2, show.sign=FALSE, angle_col=45,order=NULL){
    
    
    metadata = obj@meta.data
    
    # Roe for each cluster
    # observe_data1 = metadata %>% 
    #     group_by_(group_by1, group_by2) %>% 
    #     summarise(cell_num=n())
    observe_data1 = as.data.frame(table(metadata[,group_by1],metadata[,group_by2]))
    
    colnames(observe_data1) = c('group_by1', 'group_by2', 'cell_num')
    
    observe_data1 = observe_data1 %>% tidyr::spread(key =group_by1, value = cell_num, fill=0) %>% as.data.frame()
    
    rownames(observe_data1) = observe_data1[,1]
    observe_data1 = observe_data1[, 2:ncol(observe_data1)]
    observe_data2 = rowSums(observe_data1) - observe_data1
    
    expected_data = observe_data1
    pvalue = c()
    for(i in 1:ncol(observe_data1)){
        chisq_table = matrix(c(observe_data1[,i], observe_data2[,i]), nrow = 2, byrow = TRUE)
        chisq_res = chisq.test(chisq_table)
        expected_data[, i] = chisq_res$expected[1,]
        pvalue = c(pvalue, chisq_res$p.value)
    }
    Reo = observe_data1 / expected_data
    ####
    # plot
    
    if(show.sign=='signif'){
        annote.heatmap = t(Reo)
        annote.heatmap[annote.heatmap>1] = '+++'
        annote.heatmap[annote.heatmap<=1&annote.heatmap>0.8] = '++'
        annote.heatmap[annote.heatmap<=0.8&annote.heatmap>=0.2] = '+'
        annote.heatmap[annote.heatmap<0.2&annote.heatmap>0] = '+/-'
        annote.heatmap[annote.heatmap==0] = '-'
    }else if(show.sign=='value'){
        annote.heatmap = round(t(Reo),2)
    }else{
        annote.heatmap = t(Reo)
        annote.heatmap[annote.heatmap<Inf] = '' 
    }
    bk <- c(seq(0,0.99,by=0.01),seq(1.01,2,by=0.01))
    if(is.null(order)){
        order = colnames(Reo)
    }
    annote.heatmap = annote.heatmap[order,]
    pheatmap::pheatmap(t(Reo[,order]), 
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       treeheight_row = 0,
                       color = c(colorRampPalette(c("#053061","#2166ac","#4393c3","white"))(length(bk)/2),colorRampPalette(c("white","#fddbc7","#d6604d","#b2182b","#67001f"))(length(bk)/2)),
                       breaks = bk,
                       display_numbers = annote.heatmap,
                       cellwidth = 20, cellheight = 20,
                       fontsize = 10,
                       border_color = '#ffffff',
                       angle_col=angle_col,
                       main = 'Observed / Expected (Roe) '
    )
    res = list('Reo'= t(Reo), 'pvalue'=pvalue, 'annote'=annote.heatmap)
    return(res)
}



