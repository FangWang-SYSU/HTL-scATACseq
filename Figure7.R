library(Signac)
library(Seurat)
library(glue)
library(EnsDb.Mmusculus.v79)
library(ggsci)
library(ComplexHeatmap)
base_dir='./placeholder_project/code/'
CodePath = glue('{base_dir}/ATAC_pipeline_zyr/Code/')
source(file.path(CodePath,'Function/PipelineFuncHtml.R'))
source(file.path(CodePath,'Function/LoadsingleSample.func.R'))
source(file.path(CodePath,'Function/MotifAnalysis.func.R'))
source(file.path(CodePath,'Function/Signac.Heatmap.R'))
source(file.path(CodePath,'Function/ClusterAndTypeing.func.R'))
source(file.path(CodePath,'Function/cCREsStatistics.func.R'))
source(file.path(CodePath,'Function/Signac.Heatmap.R'))
source(file.path(CodePath,'Function/MergeSampleObj.R'))
source(file.path(CodePath,'Function/PredictCellType.R'))

setwd('./placeholder_project/code/')
source('./placeholder_project/code/jupyter_R/myFun.R')
source('./placeholder_project/code/jupyter_R/markerlist.R')
my_process_srt <- function(object, assay){
  DefaultAssay(object) <- assay
  object <- RunTFIDF(object)
  object <- FindTopFeatures(object, min.cutoff = 'q5',verbose = F)
  object <- RunSVD(object,verbose = F)
  object <- RunUMAP(object = object, reduction = 'lsi', dims = 2:30)
  object <- FindNeighbors(object = object, reduction = 'lsi', dims = 2:30 )
  object <- FindClusters(object = object, verbose = FALSE, algorithm = 3)
  return(object)
}
myRowScale<-function(m = NULL, min = -2, max = 2, limit = FALSE){
    z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
    if(limit){
        z[z > max] <- max
        z[z < min] <- min
    }
    return(z)
}
myAveragePeak <- function(object,groupby, assay=NULL, slot = 'data'){
    if(!is.null(assay)){
        DefaultAssay(object) = assay
    }else{
        assay = DefaultAssay(object)
    }
    Idents(object) = object@meta.data[[groupby]]
    clusters =  unique(Idents(object))
    mat =  GetAssayData(object,assay = assay, slot)
    
    j = 1
    for (i in clusters) {
        #i = clusters[2]
        subBarcode = WhichCells(object = object,idents = i)
        sub_mat = mat[,subBarcode,drop=F] 
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
base_path='./placeholder_analysis/round_cluster02/cCRE/'

expanded_colors <- c(
  "#D7191C", "#E6804D", "#FF6F61", "#F18A99", "#E6AB02",
  "#D95F02", "#FFB347", "#E4A8A3", "#EBEBEB", "#FADFC3",
  "#A6CD48", "#A0D6F7", "#94D3B6", "#FDAE61", "#FEE08B",
  "#A6D96A", "#1F78B4", "#7570B3", "#A6761D", "#D73027",
  "#ABDDA4", "#2B83BA", "#F46D43", "#D8B365", "#BDBDBD",
  "#A6BD8B", "#C9948B", "#D6616B", "#7B9E23", "#B69DAE",
  "#73C2A3", "#969696", "#6A9AC4", "#B97CAA", "#7AAE3C",
  "#D5B8E5", "#4DBBD5FF", "#E64B35FF", "#00A087FF", "#A65628",
  "#FF95A8FF", "#BCBD22FF", "#FDBF6F", "#3C5488FF", "#F39B7FFF",
  "#B09C85FF", "#7876B1FF", "#377EB8", "#4DAF4A", "#97A1A7FF",
  "#984EA3", "#FF7F00", "#F781BF", "#B2DF8A", "#5050FFFF",
  "#82581FFF", "#E5614CFF", "#F0E685FF", "#D595A7FF", "#CDDEB7FF",
  "#612A79FF", "#AE1F63FF", "#1B1919FF", "#99CC00FF", "#CC9900FF",
  "#9467BDFF", "#EFD500FF", "#ADE2D0FF", "#CE3D32FF", "#749B58FF",
  "#466983FF", "#BA6338FF", "#5DB1DDFF", "#802268FF", "#6BD76BFF",
  "#924822FF", "#837B8DFF", "#C75127FF", "#D58F5CFF", "#7A65A5FF",
  "#E4AF69FF", "#3B1B53FF", "#E7C76FFF", "#5A655EFF", "#33CC00FF",
  "#00CC33FF", "#00CC99FF", "#0099CCFF", "#0A47FFFF", "#4775FFFF",
  "#FFC20AFF", "#FFD147FF", "#990033FF", "#991A00FF", "#996600FF",
  "#809900FF", "#339900FF", "#00991AFF", "#009966FF", "#008099FF",
  "#003399FF", "#1A0099FF", "#660099FF", "#990080FF", "#D60047FF",
  "#FF1463FF", "#00D68FFF",
  # 50 new colors
  "#C22E29", "#F68712", "#FCD34D", "#16A34A", "#10B981",
  "#0EA5E9", "#6366F1", "#8B5CF6", "#EC4899", "#F43F5E",
  "#FF9F43", "#FFDA77", "#A3E635", "#4ADE80", "#2DD4BF",
  "#38BDF8", "#818CF8", "#A78BFA", "#E879F9", "#F87171",
  "#EF4444", "#FB923C", "#FACC15", "#84CC16", "#22D3EE",
  "#60A5FA", "#818CF8", "#C084FC", "#F472B6", "#FB7185",
  "#DC2626", "#EA580C", "#FBBF24", "#4D7C0F", "#064E3B",
  "#0F766E", "#075985", "#4338CA", "#6D28D9", "#DB2777",
  "#9D174D", "#7C2D12", "#78350F", "#365314", "#14532D",
  "#164E63", "#1E3A8A", "#312E81", "#581C87", "#701A75"
)
expanded_colors2 = c("#4dbbd5ff" ,"#E64b35ff", "#00a087ff" ,"#a65628", "#FF95A8FF","#BCBD22FF", "#fdbf6f", "#3c5488ff",
                   "#f39b7fff", "#b09c85ff", "#7876b1ff", "#377eb8",  "#4daf4a","#97A1A7FF" ,
                   "#984ea3" , "#ff7f00",  "#f781bf", "#b2df8a", "#5050FFFF", "#82581FFF" , "#E5614CFF",
                   "#F0E685FF", "#D595A7FF", "#CDDEB7FF","#612A79FF" ,"#AE1F63FF","#1B1919FF",
                   "#99CC00FF","#CC9900FF" ,"#9467BDFF", "#EFD500FF" , "#ADE2D0FF",pal_igv()(50))

pt_to_mm <- function(pt) {
  pt /2.13
}
label_size <-function(x){
    x/2.8453
}
mytheme =     theme_void()+
      theme(panel.border = element_rect(linewidth=pt_to_mm(0.25),fill=NA),
            plot.title = element_text(hjust = 0.5,face = "plain", size=8, margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt")),

            axis.title = element_text(size=7),
            axis.title.y= element_text(angle = 90,hjust=0.5, margin = margin(t = 0, r = 2, b = 0, l = 0, unit = "pt")),
            axis.title.x= element_text(margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt")),

            axis.text = element_text(size=6),
            axis.text.x = element_text(margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "pt")),
            axis.text.y = element_text(margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "pt")),

            axis.ticks = element_line(size = pt_to_mm(0.2)),
            axis.ticks.length =  unit(1, 'pt'),
            
            legend.spacing = unit(2, 'pt'),
            legend.key.width = unit(3, "mm"),  # Legend color bar width
            legend.key.height = unit(3, "mm"),    # Legend color bar height
            legend.key.spacing = unit(2, 'pt'),
            
            legend.text = element_text(size=6),
            legend.title = element_text(size=6,face = "plain"),
            legend.margin = margin(0, 0, 0, 0)
            )

umap = read.csv(glue('./placeholder_analysis/round_cluster02/merge/cell_meta.csv'))

new_celltype = read.csv('./placeholder_analysis/round_cluster02/merge/celltype_modify.csv')

head(new_celltype)

output_dir='./placeholder_output/raw_figure_output/Figure4/'

library(chromVAR)
library(motifmatchr)
library(TFBSTools)

library(BSgenome.Mmusculus.UCSC.mm10)
jaspar <- JASPAR2024::JASPAR2024()
pfm <- getMatrixSet(
  x = jaspar@db,opts = list(collection = "CORE",tax_group='vertebrates' , all_versions = FALSE))
hg38_mm10_TF_mapping = read.table('../hg38_mm10_TF_mapping.txt', sep='\t', header = TRUE)
keep_TF = unique(hg38_mm10_TF_mapping$tf_name)
keep_TF_upper = toupper(keep_TF)

dup_tf = keep_TF_upper[duplicated(keep_TF_upper)]

keep_TF = setdiff(keep_TF,dup_tf)
length(keep_TF)
keep_TF_motifname = c()
for(i in pfm@listData){
    if(i@name %in% keep_TF){
        keep_TF_motifname = c(keep_TF_motifname,i@ID)
    }
}
pfm_filter = pfm[names(pfm)%in%keep_TF_motifname]
pfm_filter


base_path='./placeholder_analysis/round_cluster02/cCRE/'
srt = readRDS(glue('{base_path}/celltype_rds/Epi_sub.rds'))
srt


srt <- AddMotifs(object = srt,
                 genome = BSgenome.Mmusculus.UCSC.mm10,
                 pfm = pfm_filter,verbose = F)

motif_meta = c()
for(m in pfm_filter@listData){
    motif_meta = rbind(motif_meta, c(m@ID, paste0(m@tags$family,collapse = '::'), m@name))
}
motif_meta = as.data.frame(motif_meta)
colnames(motif_meta) = c('motif', 'family', 'TF')
rownames(motif_meta) = motif_meta$motif

srt@meta.data[,'new_group'] = paste0(srt$Celltype_round4,':',srt$Time)

Idents(srt) = srt$new_group

srt_markers_new = c()
for(i in unique(srt$Celltype_round4)){
    message(i)
    tmp_srt = subset(srt, Celltype_round4==i)
    Idents(tmp_srt) = tmp_srt$new_group
    tmp_srt_markers = FindAllMarkers(tmp_srt, 
                                     test.use = 'wilcox', 
                                     min.pct = 0.01,
                                     logfc.threshold = 0.1,
                                     max.cells.per.ident=500)
    srt_markers_new = rbind(srt_markers_new, tmp_srt_markers)
}

head(srt_markers_new)

saveRDS(srt_markers_new, './placeholder_analysis/round_cluster02/cCRE/Epi_sub_marker_new.rds')

srt_markers = readRDS('./placeholder_analysis/round_cluster02/cCRE/Epi_sub_marker_new.rds')

load_peak<-function(name){
    base_path='./placeholder_analysis/round_cluster02/cCRE/'
    # Load differential peak matrix
    mat <- Matrix::readMM(glue("{base_path}/{name}_aggr_sparse.mtx"))
    row_names <- read.csv(glue("{base_path}/{name}_aggr_sparse_index.csv"), header = T)[,1]
    col_names <- read.csv(glue("{base_path}/{name}_aggr_sparse_columns.csv"), header = T)[,1]
    row_names = gsub(':','-',row_names)
    rownames(mat) <- row_names
    colnames(mat) <- col_names
    return(mat)
}

cell_mat = load_peak('all_celltime_markers')

cell_mat = as.matrix(cell_mat)

srt_markers_filter = srt_markers %>% 
    filter(avg_log2FC>0) %>% # ,p_val_adj<0.05
    group_by(cluster) %>%
    dplyr::top_n(100, avg_log2FC)

tmp_mat = cell_mat[,unique(srt_markers_filter$cluster)]

cell_mat_scale = t(apply(tmp_mat, 1, scale))
cell_mat_scale[is.na(cell_mat_scale)]=0
colnames(cell_mat_scale) = colnames(tmp_mat)

srt_markers_filter$Celltype = sapply(as.vector(srt_markers_filter$cluster), function(x)strsplit(x,':')[[1]][1])
srt_markers_filter$Time = sapply(as.vector(srt_markers_filter$cluster), function(x)strsplit(x,':')[[1]][2])
srt_markers_filter

celltime_order = unique(srt_markers_filter$Celltype)

srt_markers_filter$Time = factor(srt_markers_filter$Time, levels = time_levels)
srt_markers_filter$Celltype = factor(srt_markers_filter$Celltype, levels = celltime_order)
srt_markers_filter = srt_markers_filter %>% arrange(Celltype, Time)

celltime_order2 = srt_markers_filter[,c('Celltype','cluster', 'Time')] %>% distinct()

round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')

new_celltype = read.csv('./placeholder_analysis/round_cluster02/merge/celltype_modify.csv')
rownames(new_celltype) = new_celltype$Celltype_round4
new_celltype$modify_name = new_celltype$Celltype_round4_new
new_celltype$new_name = paste0(new_celltype$modify_name,'.', new_celltype$cluster_num)

celltime_order2$new_name = new_celltype[as.vector(celltime_order2$Celltype),'new_name']

unique(as.vector(celltime_order2$Celltype))

tmp_curr_celltype_color <- c(
  "#E60000", # bright red
  "#006CD1", # deep blue
  "#1A9900", # forest green
  "#FFA500", # bright orange
  "#752075", # deep purple
  "#40E0D0", # turquoise
  "#FFD700", # bright yellow
  "#004080", # navy blue
  "#FF69B4", # bright pink
  "#A0522D", # brown
  "#808000", # olive green
  "#8C9CBE", # gray blue
  "#8A2BE2", # vivid purple
  "#525252"  # dark gray
)
names(tmp_curr_celltype_color) = c('AT1','AT1_Nek10','AT1_Agt','AT1_Tldc2','AT1_Mfrp','AT2','AT2_Gfra1',
    'AT2_Sgta','Airway_Cilitated_Lypd2','Airway_Cilitated_Muc16','Airway_Acp7','Airway_Club_Sftpb',
    'Airway_Club_Scgb1a1','Airway_Spock3')

tmp_curr_celltype_color2 = tmp_curr_celltype_color
names(tmp_curr_celltype_color2) = new_celltype[names(tmp_curr_celltype_color2),'new_name']

celltime_order2$new_name = factor(celltime_order2$new_name, levels = unique(celltime_order2$new_name))

peak_order = srt_markers_filter$gene
left_ann = rowAnnotation(df=celltime_order2[,c('new_name', 'Time')],
                        col=list('new_name'=tmp_curr_celltype_color2, 'Time'=time_color2),
                        simple_anno_size=unit(0.1,'cm'))
cell_mat_scale_ht = t(cell_mat_scale[peak_order, as.vector(celltime_order2$cluster)])

ht=Heatmap(cell_mat_scale_ht,
          col=circlize::colorRamp2(c(-1,-0.5,0,0.5,1),c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")),
          #col=circlize::colorRamp2(c(0,0.1,.2),c("white","#fdae61","#d7191c")),
          cluster_rows = F, cluster_columns = F,
          show_row_names = F, show_column_names = F,
          left_annotation=left_ann,
          use_raster=T,
          border = TRUE,
          border_gp = gpar(col = "black",lwd = 0.2) ,
          column_names_gp = gpar(fontsize=7),
)

pdf(glue('{output_dir}/Epi_CellTime_cCREs_ht.pdf'), width=8, height=6)
draw(ht)
dev.off()


srt_markers_filter = srt_markers %>% 
    filter(avg_log2FC>0) %>% # ,p_val_adj<0.05
    group_by(cluster) %>%
    dplyr::top_n(2000, avg_log2FC)

dim(srt_markers_filter)


srt_close_genes = ClosestFeature(srt, StringToGRanges(srt_markers_filter$gene))

head(srt_close_genes)

length(unique(srt_close_genes$gene_name))

mm10_hg38_homoGene = read.csv('./placeholder_project/genome_annotation_file/mm10_hg38_homoGene_ensembl.csv')
#background_gene = unique(na.omit(mm10_hg38_homoGene$Human.gene.name))

close_hg38 = mm10_hg38_homoGene[mm10_hg38_homoGene$Gene.name%in%srt_close_genes$gene_name,]

close_hg38 = close_hg38[,c('Human.gene.name', 'Gene.name')] %>% unique()
close_hg38 = close_hg38[!duplicated(close_hg38$Gene.name),]
rownames(close_hg38) = close_hg38$Gene.name

length(unique(close_hg38$Human.gene.name))

background_gene = unique(na.omit(close_hg38$Human.gene.name))

library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)

background_gene_id = bitr(background_gene, 
                          fromType = 'SYMBOL', 
                          toType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')

srt_close_genes$hg38gene = close_hg38[srt_close_genes$gene_name,'Human.gene.name']

cluster_top_genes = list()
for(i in unique(srt_markers_filter$cluster)){
    tmp_genes = srt_markers_filter[srt_markers_filter$cluster==i,'gene', drop=T]
    tmp_genes = srt_close_genes[srt_close_genes$query_region%in%tmp_genes,]
    tmp_genes = unique(na.omit(tmp_genes$hg38gene))
    tmp_genes = tmp_genes[tmp_genes!='']
    if(length(tmp_genes)>0){
        tmpgene_id = bitr(tmp_genes,fromType = 'SYMBOL', toType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
        cluster_top_genes[[i]] = unique(na.omit(tmpgene_id$ENTREZID))
    }

}

head(cluster_top_genes)

srt_markers_filter2 = srt_markers %>% 
    filter(p_val_adj<0.05, avg_log2FC>0) %>% #
    group_by(cluster) %>%
    dplyr::top_n(2000, avg_log2FC)

cluster_top_genes2 = list()
for(i in unique(srt_markers_filter2$cluster)){
    tmp_genes = srt_markers_filter2[srt_markers_filter2$cluster==i,'gene', drop=T]
    tmp_genes = srt_close_genes[srt_close_genes$query_region%in%tmp_genes,]
    tmp_genes = unique(na.omit(tmp_genes$hg38gene))
    tmp_genes = tmp_genes[tmp_genes!='']
    cluster_top_genes2[[i]] = unique(na.omit(tmp_genes))
}

cluster_enrich_Disease = lapply(cluster_top_genes, function(x){
    disease_enrich = enrichDGN(
        x,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        universe=background_gene_id$ENTREZID,
        #minGSSize = 1,
        #maxGSSize = 5000,
        qvalueCutoff = 1,
        readable = FALSE)
})

all_disease = do.call(rbind, lapply(names(cluster_enrich_Disease), function(x){
    xx = cluster_enrich_Disease[[x]]@result
    xx$cluster = x
    return(xx)
}))

all_disease[, c('Celltype', 'Time')] = t(sapply(all_disease$cluster, function(x)strsplit(x,':')[[1]]))

all_disease_sig = all_disease %>%
    filter(p.adjust<0.05) %>%
    filter(Count>5)

head(all_disease_sig)

saveRDS(all_disease, './placeholder_analysis/round_cluster02/cCRE/Epi_sub_marker_new_Disease.rds')

all_disease = readRDS('./placeholder_analysis/round_cluster02/cCRE/Epi_sub_marker_new_Disease.rds')

lung_res = all_disease[grepl('lung',all_disease$Description), ]
# lung_res[, c('Celltype', 'Time')] = t(sapply(lung_res$cluster, function(x)strsplit(x,':')[[1]]))

new_celltype = read.csv('./placeholder_analysis/round_cluster02/merge/celltype_modify.csv')
rownames(new_celltype) = new_celltype$Celltype_round4
new_celltype$modify_name = new_celltype$Celltype_round4_new
new_celltype$new_name = paste0(new_celltype$modify_name,'.', new_celltype$cluster_num)

new_celltype[new_celltype$new_name=='AT1.116', 'new_name'] = 'AT1_Spock2.116'
new_celltype[new_celltype$new_name=='AT2.124', 'new_name'] = 'AT2_Sftpc.124'

lung_res$new_name = new_celltype[lung_res$Celltype,'new_name']

tmp_curr_celltype_color =  c(
  "#E41A1C", # bright red
  "#377EB8", # bright blue
  "#4DAF4A", # bright green
  "#984EA3", # deep purple
  "#FF7F00", # bright orange
  "#FFD92F", # bright yellow
  "#A6CEE3", # lake blue
  "#1F78B4", # deep blue
  "#B2DF8A", # light green
  "#33A02C", # deep green
  "#FB9A99", # light red
  "#CAB2D6", # light purple
  "#6A3D9A", # dark purple
  "#FF1493"  # vivid pink
)
tmp_curr_celltype_color <- c(
  "#E60000", # bright red
  "#006CD1", # deep blue
  "#1A9900", # forest green
  "#FFA500", # bright orange
  "#752075", # deep purple
  "#40E0D0", # turquoise
  "#FFD700", # bright yellow
  "#004080", # navy blue
  "#FF69B4", # bright pink
  "#A0522D", # brown
  "#808000", # olive green
  "#8C9CBE", # gray blue
  "#8A2BE2", # vivid purple
  "#525252"  # dark gray
)
names(tmp_curr_celltype_color) = c('AT1','AT1_Nek10','AT1_Agt','AT1_Tldc2','AT1_Mfrp','AT2','AT2_Gfra1',
    'AT2_Sgta','Airway_Cilitated_Lypd2','Airway_Cilitated_Muc16','Airway_Acp7','Airway_Club_Sftpb',
    'Airway_Club_Scgb1a1','Airway_Spock3')

tmp_curr_celltype_color2 = tmp_curr_celltype_color
names(tmp_curr_celltype_color2) = new_celltype[names(tmp_curr_celltype_color2),'new_name']

library(seriation)

register_DendSer()

tmp = lung_res %>% 
    filter(pvalue<0.01,Count>3) %>%
    mutate(cell_disease=paste0(new_name,':', Description)) %>%
    mutate(FoldEnrichment = log2(FoldEnrichment)) %>%
    reshape2::acast(cell_disease~Time, value.var = 'FoldEnrichment', fill=0)
tmp = tmp[, time_levels]
o = seriate(tmp,method='Heatmap', 
            seriation_method = 'DendSer_BAR',dist_fun=dist)

tmp = tmp[get_order(o,1),]
disease_order = rownames(tmp)

stage = cutree(o[[1]], 15)

ht = Heatmap(tmp,
        cluster_rows=F, 
        cluster_columns=F,
        show_row_names=T,
        show_column_dend=F,
        show_row_dend = F,
        row_split = factor(stage[rownames(tmp)],c(1,13,12,6,15,2,10,14,7,11,3,4,9,8,5))
       )
draw(ht)

disease_order = rownames(tmp)[unlist(row_order(ht))]

tmp_lung_res_gene = lung_res %>%
    filter(pvalue<0.01,Count>5) %>%
    #filter(Description%in%keep_trait$label) %>%
    arrange(Description, new_name)
gene_list = sapply(tmp_lung_res_gene$geneID, function(x){
    background_gene_id[background_gene_id$ENTREZID%in%strsplit(x,'/')[[1]], 'SYMBOL']
})
names(gene_list) = paste0(tmp_lung_res_gene$new_name, ':',tmp_lung_res_gene$Description)

names(gene_list)

gene_list['Airway_Cilitated_Muc16:Aplasia/Hypoplasia of the lungs']

length(unique(unlist(gene_list)))

#keep_trait = all_trait_sq %>% filter(Time>0 , Celltype>0)
a=lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    #filter(Description%in%keep_trait$label) %>%
    arrange(Description, Celltype) %>%
    mutate(cell_disease=paste0(new_name,':', Description)) %>%
    mutate(cell_disease=factor(cell_disease, rev(disease_order))) %>%
    mutate(Time=factor(Time, time_levels)) %>%
    ggplot(aes(x=Time, y=cell_disease,fill=log2(FoldEnrichment)))+
    #geom_tile(color='black', size=0.1)+
    geom_point(aes(size=-log10(pvalue)),shape=22)+
    #geom_point(aes(size=-log10(Coefficient_P_value)), shape=22, stroke=0.5)+
    #geom_tile()+
    scale_size_continuous(range = c(1,4))+
    #facet_wrap(~Trait, ncol=3)+
    #scale_fill_gradientn(colours = c("#08519c","#deebf7",'white', 'red', 'darkred'))+
    scale_fill_gradientn(colours = c("#deebf7","#08519c"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1),
          axis.text = element_text(size=6))

options(repr.plot.width=10, repr.plot.height=10)
a

ggsave(glue('{output_dir}/Epi_function2(cell and time dependent).pdf'), a,  
       width=180, height = 150, units='mm', dpi=600, bg='transparent')

head(lung_res)

library(ggbeeswarm)
round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')

# number of cell types affected by diseases at each time point
options(repr.plot.width=12, repr.plot.height=6)
a_top=lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    #filter(Description%in%keep_trait$label) %>%
    arrange(Description, Celltype) %>%
    group_by(Time) %>%
    summarise(n=length(unique(new_name))) %>%
    mutate(Time=factor(Time, time_levels)) %>%
    ggplot(aes(x=Time,y=n, fill=Time))+
    geom_bar(stat='identity', width=0.6)+
    scale_fill_manual(values = time_color2)+
    labs(y='Celltype_num')+
    theme_bw()+
    theme(legend.position = 'none', panel.grid = element_blank())

# number of cell types affected by diseases at each time point
options(repr.plot.width=12, repr.plot.height=6)
b_top=lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    filter(!duplicated(paste0(Description, Time)))%>%
    #filter(Description%in%keep_trait$label) %>%
    arrange(Description, Celltype) %>%
    mutate(Time=factor(Time, time_levels)) %>%
    ggplot(aes(x=Time, fill=Time))+
    geom_bar(stat='count', width=0.6)+
    scale_fill_manual(values = time_color2)+
    labs(y='Disease_num')+
    theme_bw()+
    theme(legend.position = 'none', panel.grid = element_blank())
a_top/b_top
ggsave(glue('{output_dir}/Epi_function2(cell and time dependent-count time points).pdf'), a_top/b_top,  
       width=6, height = 4, units='in', dpi=600, bg='transparent')

# number of cells affected per disease
options(repr.plot.width=12, repr.plot.height=6)
a_right=lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    filter(!duplicated(paste0(Description, Celltype)))%>%
    #filter(Description%in%keep_trait$label) %>%
    arrange(Description, Celltype) %>%
    mutate(Time=factor(Time, time_levels)) %>%
    ggplot(aes(x=Description, fill=Description))+
    geom_bar(stat='count', width=0.6)+
    scale_fill_igv()+
    scale_y_continuous(breaks=0:12)+
    labs(y='Celltype_num')+
    theme_bw()+coord_flip()+
    theme(legend.position = 'none', panel.grid = element_blank())

# number of time points affected per disease
options(repr.plot.width=12, repr.plot.height=6)
b_right=lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    filter(!duplicated(paste0(Description, Time)))%>%
    #filter(Description%in%keep_trait$label) %>%
    arrange(Description, Celltype) %>%
    mutate(Time=factor(Time, time_levels)) %>%
    ggplot(aes(x=Description, fill=Description))+
    geom_bar(stat='count', width=0.6)+
    scale_fill_igv()+
      scale_y_continuous(breaks=0:12)+
    labs(y='Time_num')+
    theme_bw()+coord_flip()+
    theme(legend.position = 'none', panel.grid = element_blank())
a_right+b_right
ggsave(glue('{output_dir}/Epi_function2(cell and time dependent-count diseases).pdf'), a_right+b_right,  
       width=10, height = 4, units='in', dpi=600, bg='transparent')

xx=lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    filter(!duplicated(paste0(Description, Celltype)))%>%
    group_by(Description) %>%
    summarise(num=length(unique(Celltype))) %>%
    #filter(Description%in%keep_trait$label) %>%
    arrange(num)

length(unique(lung_res$Celltype))

a_right=lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    filter(!duplicated(paste0(Description, Celltype)))%>%
    #filter(Description%in%keep_trait$label) %>%
    mutate(Description=factor(Description, xx$Description)) %>%
    ggplot(aes(x=Description, fill=Celltype))+
    geom_bar(stat='count', width=0.6)+
    scale_fill_manual(values = tmp_curr_celltype_color)+
    scale_y_continuous(breaks=0:12)+
    labs(y='Celltype_num')+
    theme_bw()+coord_flip()+
    theme(legend.position = 'none', panel.grid = element_blank())
options(repr.plot.width=6, repr.plot.height=6)
a_right

# number of cells affected per disease
options(repr.plot.width=12, repr.plot.height=6)
a=lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    filter(!duplicated(paste0(Description, Celltype)))%>%
    group_by(Celltype,new_name) %>%
    summarise(num=length(unique(Description))) %>%
    as.data.frame() %>%
    arrange(num) %>%
    mutate(new_name=factor(new_name, new_name)) %>%
    ggplot(aes(x=new_name, y=num, fill=new_name))+
    geom_bar(stat='identity', width=0.6)+
    scale_fill_manual(values = tmp_curr_celltype_color2)+
    labs(y='Disease_num')+
    theme_bw()+coord_flip()+
    theme(legend.position = 'none', panel.grid = element_blank())

# number of time points affected per disease
options(repr.plot.width=12, repr.plot.height=6)
b=lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    filter(!duplicated(paste0(Description, Celltype, Time)))%>%
    group_by(Celltype,new_name) %>%
    summarise(num=length(unique(Time))) %>%
    as.data.frame() %>%
    arrange(num) %>%
    mutate(new_name=factor(new_name, new_name)) %>%
    ggplot(aes(x=new_name, y=num, fill=new_name))+
    geom_bar(stat='identity', width=0.6)+
    scale_fill_manual(values = tmp_curr_celltype_color2)+
    labs(y='Time_num')+
    theme_bw()+coord_flip()+
    theme(legend.position = 'none', panel.grid = element_blank())
a+b
ggsave(glue('{output_dir}/Epi_function2(cell and time dependent-count cell types).pdf'), a+b,  
       width=10, height = 4, units='in', dpi=600, bg='transparent')

library(scatterpie)

pie_data = lung_res %>%
    filter(pvalue<0.01,Count>3) %>%
    #filter(Description%in%keep_trait$label) %>%
    arrange(Description, new_name) %>%
    dcast(Description+Time~new_name, value.var = 'Count', fill=0) 
pie_data_cols = colnames(pie_data)[-(1:2)]
pie_data$group = 1:nrow(pie_data)

options(repr.plot.width=12, repr.plot.height=12)
a2=pie_data %>%
    mutate(Description=factor(Description, xx$Description)) %>%
    mutate(Time=factor(Time, time_levels)) %>%
    mutate(Time=as.numeric(Time)) %>%
    mutate(Description=as.numeric(Description)) %>%
    ggplot(aes(x=Time, y=Description))+
    geom_scatterpie(
        aes(x=Time, y=Description,group=group, r=0.4),   # r = radius
        cols=pie_data_cols, color=NA
    ) +
    coord_equal()+
    scale_x_continuous(breaks=1:15, labels = time_levels)+
    scale_y_continuous(breaks=1:length(unique(xx$Description)), labels = xx$Description)+
    scale_fill_manual(values = tmp_curr_celltype_color2)+
    # #scale_fill_igv()+
    theme_bw()
a2

library(igraph)
round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')

tmp_lung_res0 =lung_res%>% filter(pvalue<0.01,Count>3) %>% mutate(Celltype = new_name) 

write.table(tmp_lung_res0, glue('{output_dir}/Lung_enrich.txt'), quote=F, sep='\t')

tmp_lung_res0 = read.table(glue('{output_dir}/Lung_enrich.txt'), header=T, sep='\t')

tmp_lung_res0_P7 = tmp_lung_res0 %>% filter(Time=='P7')

options(repr.plot.width=8, repr.plot.height=8)
#tmp_time = unlist(mult_dis[mult_dis$Celltype==tmp_celltype, 'Time'])

tmp_lung_res = tmp_lung_res0 %>% filter(Description%in%tmp_lung_res0_P7$Description, Celltype%in%tmp_lung_res0_P7$Celltype) 
#tmp_lung_res$Time = paste0(tmp_lung_res$Time, tmp_lung_res$Celltype)
#other_time = setdiff(time_levels,unique(tmp_lung_res$Time))
#tmp_net = data.frame('Celltype'=tmp_celltype, 'Time'=other_time)

net_work_data = as.matrix(tmp_lung_res[,c('Description', 'Time')])
net_work_data = rbind(net_work_data, as.matrix(tmp_lung_res[,c('Time', 'Celltype')]))
net_work_data = as.data.frame(net_work_data)
#net_work_data = rbind(net_work_data, tmp_net)
net_work_data = net_work_data %>% unique()

g <- graph_from_data_frame(net_work_data , directed = T)


nodes_color = sapply(names(V(g)), function(x){
    if(x %in% net_work_data$Time){
        return(time_color2[x])
    }else if(x %in% net_work_data$Celltype){
        return(tmp_curr_celltype_color2[x])
    }else{
        return('gray')
    }
})
nodes_size = sapply(names(V(g)), function(x){
    if(x %in% tmp_lung_res$Time){
        return(12)
    }else if(x %in% tmp_lung_res$Celltype){
        return(12)
    }else{
        return(12)
    }
})

set.seed(1234)
#g_order = c(tmp_celltype, time_levels, unique(tmp_lung_res$Description))
#    order_idx = match(g_order, V(g)$name)

ly = layout_as_tree(g)#[order_idx,]
ly_roate = cbind(-ly[,2], ly[,1])
plot(g, edge.arrow.size = 0.3, vertex.color = nodes_color, vertex.label.color = "black",
    vertex.size=nodes_size,
        vertex.label.dist=2, 
        vertex.label.degree=0,
        layout=ly_roate)#

write.table(net_work_data, glue('{output_dir}/Lung_enrich_network.txt'), quote=F, sep='\t', row.names=F, col.names=T)

pdf(glue('{output_dir}/P7relationship network.pdf'), width=12, height=6)

options(repr.plot.width=12, repr.plot.height=10)
plot(g, edge.arrow.size = 0.3, vertex.color = nodes_color, vertex.label.color = "black",
    vertex.size=nodes_size,
        vertex.label.dist=2, 
        vertex.label.degree=0,
        layout=ly_roate)#

dev.off()

options(repr.plot.width=6, repr.plot.height=4)
a=tmp_lung_res0 %>%
    group_by(Celltype, Time) %>%
    dplyr::summarise(num=n()) %>%
    mutate(num = ifelse(num>2, 2, num)) %>%
    ggplot(aes(x=num))+
    geom_bar(stat='count', width=0.6)+
    scale_x_continuous(breaks = c(1,2), labels = c(1,'2+'))+
    labs(x='Disease num', y='Celltype+Time Count')+
    theme_bw()
a
ggsave(glue('{output_dir}/Epi_function(sameCsameTdifferentD).pdf'), a,  
       width=40, height = 50, units='mm', dpi=600, bg='transparent')

# affect multiple diseases
mult_dis = tmp_lung_res0 %>%
    group_by(Celltype, Time) %>%
    dplyr::summarise(num=n()) #%>% 
    #filter(num>1)

pdf(glue('{output_dir}/sameCsameTdifferentDrelationship network.pdf'), width=12, height=6)

options(repr.plot.width=12, repr.plot.height=10)
par(mfrow=c(2,4))
for(tmp_celltype in unique(mult_dis$Celltype)){
    # tmp_celltype = as.character(mult_dis[i, 'Celltype'])
    tmp_time = unlist(mult_dis[mult_dis$Celltype==tmp_celltype, 'Time'])
    
    tmp_lung_res = tmp_lung_res0 %>% filter(Celltype==tmp_celltype, Time%in%tmp_time) 
    #other_time = setdiff(time_levels,unique(tmp_lung_res$Time))
    #tmp_net = data.frame('Celltype'=tmp_celltype, 'Time'=other_time)

    net_work_data = as.matrix(tmp_lung_res[,c('Celltype', 'Time')])
    net_work_data = rbind(net_work_data, as.matrix(tmp_lung_res[,c('Time', 'Description')]))
    net_work_data = as.data.frame(net_work_data)
    #net_work_data = rbind(net_work_data, tmp_net)
    net_work_data = net_work_data %>% unique()

    g <- graph_from_data_frame(net_work_data , directed = T)


    nodes_color = sapply(names(V(g)), function(x){
        if(x %in% net_work_data$Time){
            return(time_color2[x])
        }else if(x %in% net_work_data$Celltype){
            return(tmp_curr_celltype_color2[x])
        }else{
            return('gray')
        }
    })
    nodes_size = sapply(names(V(g)), function(x){
        if(x %in% tmp_lung_res$Time){
            return(12)
        }else if(x %in% tmp_lung_res$Celltype){
            return(18)
        }else{
            return(6)
        }
    })

    set.seed(1234)
    #g_order = c(tmp_celltype, time_levels, unique(tmp_lung_res$Description))
    #    order_idx = match(g_order, V(g)$name)

    ly = layout_as_tree(g)#[order_idx,]
    ly_roate = cbind(-ly[,2], ly[,1])
    plot(g, edge.arrow.size = 0.3, vertex.color = nodes_color, vertex.label.color = "black",
        vertex.size=nodes_size,
         vertex.label.dist=2, 
         vertex.label.degree=0,
         layout=ly_roate)
}

dev.off()

# example

celltye_time_disease_genes_1 = list()

for(k in 1:nrow(mult_dis)){
    tc = mult_dis$Celltype[k]
    tt = mult_dis$Time[k]
    tmp_lung_res = tmp_lung_res0 %>% filter(Celltype==tc, Time==tt) 

    tmp_same_gene = NULL
    tmp_unique_gene = list()
    for(j in tmp_lung_res$Description){
        i = tmp_lung_res[tmp_lung_res$Description==j, 'geneID']
        tmp_gene = strsplit(i, '/')[[1]]
        if(is.null(tmp_same_gene)){
            tmp_same_gene = tmp_gene
        }else{
            tmp_same_gene = intersect(tmp_same_gene, tmp_gene)
        }
        tmp_unique_gene[[j]] = tmp_gene

    }
    tmp_unique_gene_name = list()
    for(i in names(tmp_unique_gene)){
        tmp_gene = setdiff(tmp_unique_gene[[i]], tmp_same_gene)
        tmpgene_id = bitr(tmp_gene,toType = 'SYMBOL', fromType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
        tmp_unique_gene_name[[i]] = tmpgene_id$SYMBOL
    }
    #tmp_unique_gene_name

    tmpgene_id = bitr(tmp_same_gene,toType = 'SYMBOL', fromType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
    #tmpgene_id$SYMBOL
    celltye_time_disease_genes_1[[paste0(tc,':',tt)]] = list('unique'=tmp_unique_gene_name,
                                                             'same'=tmpgene_id$SYMBOL)
}

celltye_time_disease_genes_1[['AT1_Nek10.112:P0']]



options(repr.plot.width=6, repr.plot.height=4)
a=tmp_lung_res0 %>%
    group_by(Celltype, Description) %>%
    dplyr::summarise(num=n()) %>%
    mutate(num = ifelse(num>2, 2, num)) %>%
    ggplot(aes(x=num))+
    geom_bar(stat='count', width=0.6)+
    scale_x_continuous(breaks = c(1,2), labels = c(1,'2+'))+
    labs(x='Disease num', y='Celltype+Time Count')+
    theme_bw()
a
ggsave(glue('{output_dir}/Epi_function(sameCdifferentTsameD).pdf'), a,  
       width=40, height = 50, units='mm', dpi=600, bg='transparent')

# affect single disease
single_dis = tmp_lung_res0 %>%
    group_by(Celltype, Description) %>%
    dplyr::summarise(num=n()) %>% 
    filter(num>1)

single_dis

pdf(glue('{output_dir}/sameCdifferentTsameDrelationship network.pdf'), width=12, height=6)
options(repr.plot.width=12, repr.plot.height=10)
par(mfrow=c(2,4))
for(tmp_celltype in unique(single_dis$Celltype)){

    # tmp_celltype = as.character(mult_dis[i, 'Celltype'])
    tmp_time = unlist(single_dis[single_dis$Celltype==tmp_celltype, 'Description'])
    #tmp_celltype='AT1_Tldc2'
    #tmp_time = 'stage, non-small cell lung cancer'
    tmp_lung_res = tmp_lung_res0 %>% filter(Celltype==tmp_celltype, Description%in%tmp_time) 
    #other_time = setdiff(time_levels,unique(tmp_lung_res$Time))
    #tmp_net = data.frame('Celltype'=tmp_celltype, 'Time'=other_time)
    #print(tmp_lung_res)
    net_work_data = as.matrix(tmp_lung_res[,c('Celltype', 'Time')])
    net_work_data = rbind(net_work_data, as.matrix(tmp_lung_res[,c('Time', 'Description')]))
    net_work_data = as.data.frame(net_work_data)
    #net_work_data = rbind(net_work_data, tmp_net)
    net_work_data = net_work_data %>% unique()

    g <- graph_from_data_frame(net_work_data , directed = T)


    nodes_color = sapply(names(V(g)), function(x){
        if(x %in% net_work_data$Time){
            return(time_color2[x])
        }else if(x %in% net_work_data$Celltype){
            return(round4_color[x])
        }else{
            return('gray')
        }
    })
    nodes_size = sapply(names(V(g)), function(x){
        if(x %in% tmp_lung_res$Time){
            return(12)
        }else if(x %in% tmp_lung_res$Celltype){
            return(18)
        }else{
            return(6)
        }
    })

    set.seed(1234)
    #g_order = c(tmp_celltype, time_levels, unique(tmp_lung_res$Description))
    #    order_idx = match(g_order, V(g)$name)

    ly = layout_as_tree(g)#[order_idx,]
    ly_roate = cbind(-ly[,2], ly[,1])
    plot(g, edge.arrow.size = 0.3, vertex.color = nodes_color, vertex.label.color = "black",
        vertex.size=nodes_size,
         vertex.label.dist=2, 
         vertex.label.degree=0,
         layout=ly_roate)
}
dev.off()



celltye_time_disease_genes_2 = list()

for(k in 1:nrow(single_dis)){
    tc = single_dis$Celltype[k]
    tt = single_dis$Description[k]
    tmp_lung_res = tmp_lung_res0 %>% filter(Celltype==tc, Description==tt) 

    tmp_same_gene = NULL
    tmp_unique_gene = list()
    for(j in tmp_lung_res$Time){
        i = tmp_lung_res[tmp_lung_res$Time==j, 'geneID']
        tmp_gene = strsplit(i, '/')[[1]]
        tmp_gene = bitr(tmp_gene,toType = 'SYMBOL', fromType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
        tmp_gene = tmp_gene$SYMBOL
        #tmp_gene = intersect(tmp_gene, cluster_top_genes2[[paste0(tc,':',j)]])
        if(is.null(tmp_same_gene)){
            tmp_same_gene = tmp_gene
        }else{
            tmp_same_gene = intersect(tmp_same_gene, tmp_gene)
        }
        tmp_unique_gene[[j]] = tmp_gene

    }
#     tmp_unique_gene_name = list()
#     for(i in names(tmp_unique_gene)){
#         tmp_gene = setdiff(tmp_unique_gene[[i]], tmp_same_gene)
        
#         tmp_unique_gene_name[[i]] = tmpgene_id$SYMBOL
#     }
    #tmp_unique_gene_name

    # tmpgene_id = bitr(tmp_same_gene,toType = 'SYMBOL', fromType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
    #tmpgene_id$SYMBOL
    celltye_time_disease_genes_2[[paste0(tc,':',tt)]] = list('unique'=tmp_unique_gene,
                                                             'same'=tmp_same_gene)
}

celltye_time_disease_genes_2

options(repr.plot.width=6, repr.plot.height=4)
a=tmp_lung_res0 %>%
    group_by(Description) %>%
    dplyr::summarise(num=length(unique(Celltype))) %>% 
    mutate(num = ifelse(num>2, 2, num)) %>%
    ggplot(aes(x=num))+
    geom_bar(stat='count', width=0.6)+
    scale_x_continuous(breaks = c(1,2), labels = c(1,'2+'))+
    labs(x='Disease num', y='Celltype+Time Count')+
    theme_bw()
a
ggsave(glue('{output_dir}/Epi_function(differentCsameD).pdf'), a,  
       width=40, height = 50, units='mm', dpi=600, bg='transparent')

tmp_lung_res0 %>%
    group_by(Description) %>%
    dplyr::summarise(num=length(unique(Celltype)))
# affect single disease
single_dis2 = tmp_lung_res0 %>%
    group_by(Description) %>%
    dplyr::summarise(num=length(unique(Celltype))) #%>% 
    #filter(num>1)
single_dis2

tmp_lung_res0[tmp_lung_res0$Description=='Adenocarcinoma of lung, stage I',]

pdf(glue('{output_dir}/differentCsameDrelationship network.pdf'), width=12, height=12)
options(repr.plot.width=12, repr.plot.height=10)
par(mfrow=c(3,3))
for(tmp_Description in unique(single_dis2$Description)){
    # tmp_celltype = as.character(mult_dis[i, 'Celltype'])
    # tmp_time = unlist(single_dis2[single_dis2$Celltype==tmp_celltype, 'Description'])
    
    tmp_lung_res = tmp_lung_res0 %>% filter(Description%in%tmp_Description) 
    tmp_lung_res$Celltype = paste0(tmp_lung_res$Celltype, ':', 1:nrow(tmp_lung_res))
    tmp_lung_res$Time = paste0(tmp_lung_res$Time, ':', 1:nrow(tmp_lung_res))
    
    #other_time = setdiff(time_levels,unique(tmp_lung_res$Time))
    #tmp_net = data.frame('Celltype'=tmp_celltype, 'Time'=other_time)

    net_work_data = as.matrix(tmp_lung_res[,c('Celltype', 'Time')])
    net_work_data = rbind(net_work_data, as.matrix(tmp_lung_res[,c('Time', 'Description')]))
    net_work_data = as.data.frame(net_work_data)
    #net_work_data = rbind(net_work_data, tmp_net)
    net_work_data = net_work_data %>% unique()

    g <- graph_from_data_frame(net_work_data , directed = T)


    nodes_color = sapply(names(V(g)), function(x){
        
        if(x %in% net_work_data$Time){
            x = strsplit(x, ':')[[1]][1]
            return(time_color2[x])
        }else if(x %in% net_work_data$Celltype){
            x = strsplit(x, ':')[[1]][1]
            return(tmp_curr_celltype_color2[x])
        }else{
            return('gray')
        }
    })
    nodes_size = sapply(names(V(g)), function(x){
        x = strsplit(x, ':')[[1]][1]
        if(x %in% tmp_lung_res$Time){
            return(12)
        }else if(x %in% tmp_lung_res$Celltype){
            return(18)
        }else{
            return(6)
        }
    })

    set.seed(1234)
    #g_order = c(tmp_celltype, time_levels, unique(tmp_lung_res$Description))
    #    order_idx = match(g_order, V(g)$name)

    ly = layout_as_tree(g)#[order_idx,]
    ly_roate = cbind(-ly[,2], ly[,1])
    plot(g, edge.arrow.size = 0.3, vertex.color = unname(nodes_color), vertex.label.color = "black",
        vertex.size=nodes_size,
         #vertex.label.dist=2, 
         #vertex.label.degree=0,
         layout=ly_roate)
}
dev.off()

celltye_time_disease_genes_3 = list()

for(k in 1:nrow(single_dis2)){
    tmp_Description = single_dis2$Description[k]
    tmp_lung_res = tmp_lung_res0 %>% filter(Description%in%tmp_Description) 
    
    tmp_same_gene = NULL
    tmp_unique_gene = list()

    for(j in 1:nrow(tmp_lung_res)){
        tc = tmp_lung_res$Celltype[j]
        tt = tmp_lung_res$Time[j]
        
        i = tmp_lung_res[j, 'geneID']
        tmp_gene = strsplit(i, '/')[[1]]
        if(is.null(tmp_same_gene)){
            tmp_same_gene = tmp_gene
        }else{
            tmp_same_gene = intersect(tmp_same_gene, tmp_gene)
        }
        tmp_unique_gene[[paste0(tc,':',tt)]] = tmp_gene

    }
    tmp_unique_gene_name = list()
    for(i in names(tmp_unique_gene)){
        tmp_gene = setdiff(tmp_unique_gene[[i]], tmp_same_gene)
        tmpgene_id = bitr(tmp_gene,toType = 'SYMBOL', fromType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
        tmp_unique_gene_name[[i]] = tmpgene_id$SYMBOL
    }
    #tmp_unique_gene_name
    tmpgene_id = bitr(tmp_same_gene,toType = 'SYMBOL', fromType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
    #tmpgene_id$SYMBOL
    celltye_time_disease_genes_3[[tmp_Description]] = list('unique'=tmp_unique_gene_name,
                                                             'same'=tmpgene_id$SYMBOL)
}

celltye_time_disease_genes_3

options(repr.plot.width=6, repr.plot.height=4)
a=tmp_lung_res0 %>%
    group_by(Time, Description) %>%
    dplyr::summarise(num=n()) %>%
    mutate(num = ifelse(num>2, 2, num)) %>%
    ggplot(aes(x=num))+
    geom_bar(stat='count', width=0.6)+
    scale_x_continuous(breaks = c(1,2), labels = c(1,'2+'))+
    labs(x='Disease num', y='Celltype+Time Count')+
    theme_bw()
a
ggsave(glue('{output_dir}/Epi_function(sameTdifferentCsameD).pdf'), a,  
       width=40, height = 50, units='mm', dpi=600, bg='transparent')

# affect single disease
single_dis3 = tmp_lung_res0 %>%
    group_by(Time, Description) %>%
    dplyr::summarise(num=n()) #%>% 
    #filter(num>1)

length(unique(single_dis3$Time))

pdf(glue('{output_dir}/sameTdifferentCsameDrelationship network.pdf'), width=12, height=12)
options(repr.plot.width=12, repr.plot.height=16)
par(mfrow=c(4,4))
for(tmp_celltype in unique(single_dis3$Time)){
    # tmp_celltype = as.character(mult_dis[i, 'Celltype'])
    tmp_time = unlist(single_dis3[single_dis3$Time==tmp_celltype, 'Description'])
    
    tmp_lung_res = tmp_lung_res0 %>% filter(Time==tmp_celltype, Description%in%tmp_time) 
    #other_time = setdiff(time_levels,unique(tmp_lung_res$Time))
    #tmp_net = data.frame('Celltype'=tmp_celltype, 'Time'=other_time)

    net_work_data = as.matrix(tmp_lung_res[,c('Time','Celltype')])
    net_work_data = rbind(net_work_data, as.matrix(tmp_lung_res[,c('Celltype', 'Description')]))
    net_work_data = as.data.frame(net_work_data)
    #net_work_data = rbind(net_work_data, tmp_net)
    net_work_data = net_work_data %>% unique()

    g <- graph_from_data_frame(net_work_data , directed = T)


    nodes_color =  sapply(names(V(g)), function(x){
     
        if(x %in% names(time_color2)){
            return(time_color2[x])
        }else if(x %in% names(tmp_curr_celltype_color2)){
            return(tmp_curr_celltype_color2[x])
        }else{
            return('gray')
        }
    })
    nodes_size = sapply(names(V(g)), function(x){
        if(x %in% tmp_lung_res$Time){
            return(12)
        }else if(x %in% tmp_lung_res$Celltype){
            return(18)
        }else{
            return(6)
        }
    })

    set.seed(1234)
    #g_order = c(tmp_celltype, time_levels, unique(tmp_lung_res$Description))
    #    order_idx = match(g_order, V(g)$name)

    ly = layout_as_tree(g)#[order_idx,]
    ly_roate = cbind(-ly[,2], ly[,1])
    plot(g, edge.arrow.size = 0.3, vertex.color = nodes_color, vertex.label.color = "black",
        vertex.size=nodes_size,
         vertex.label.dist=2, 
         vertex.label.degree=0,
         layout=ly_roate)
}
dev.off()

single_dis3

celltye_time_disease_genes_4 = list()

for(k in 1:nrow(single_dis3)){
    tmp_Description = single_dis3$Description[k]
    tt = single_dis3$Time[k]
    tmp_lung_res = tmp_lung_res0 %>% filter(Description%in%tmp_Description) %>% filter(Time%in%tt) 
    
    tmp_same_gene = NULL
    tmp_unique_gene = list()

    for(j in 1:nrow(tmp_lung_res)){
        tc = tmp_lung_res$Celltype[j]
        i = tmp_lung_res[j, 'geneID']
        tmp_gene = strsplit(i, '/')[[1]]
        if(is.null(tmp_same_gene)){
            tmp_same_gene = tmp_gene
        }else{
            tmp_same_gene = intersect(tmp_same_gene, tmp_gene)
        }
        tmp_unique_gene[[paste0(tmp_Description,':',tc)]] = tmp_gene

    }
    tmp_unique_gene_name = list()
    for(i in names(tmp_unique_gene)){
        tmp_gene = setdiff(tmp_unique_gene[[i]], tmp_same_gene)
        tmpgene_id = bitr(tmp_gene,toType = 'SYMBOL', fromType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
        tmp_unique_gene_name[[i]] = tmpgene_id$SYMBOL
    }
    #tmp_unique_gene_name
    tmpgene_id = bitr(tmp_same_gene,toType = 'SYMBOL', fromType = c('ENTREZID'), OrgDb = 'org.Hs.eg.db')
    #tmpgene_id$SYMBOL
    celltye_time_disease_genes_4[[tt]] = list('unique'=tmp_unique_gene_name,
                                              'same'=tmpgene_id$SYMBOL)
}

celltye_time_disease_genes_4

myfragments=list()
for(x in Fragments(srt)){
    if(length(x@cells)>0){
        xname = x@path
        if(xname %in% names(myfragments)){
            myfragments[[xname]]@cells = c(myfragments[[xname]]@cells, x@cells)
        }else{
            myfragments[[xname]] = x
        }
    }
}
myfragments = unname(myfragments)

Fragments(srt) = NULL

Fragments(srt) = myfragments

srt_markers_filter2 = srt_markers %>% 
    filter( avg_log2FC>0) %>% #p_val_adj<0.05,
    group_by(cluster) %>%
    dplyr::top_n(2000, avg_log2FC)

cluster_top_genes2 = list()
for(i in unique(srt_markers_filter2$cluster)){
    tmp_genes = srt_markers_filter2[srt_markers_filter2$cluster==i,'gene', drop=T]
    tmp_genes = srt_close_genes[srt_close_genes$query_region%in%tmp_genes,]
    tmp_genes = unique(na.omit(tmp_genes$hg38gene))
    tmp_genes = tmp_genes[tmp_genes!='']
    cluster_top_genes2[[i]] = unique(na.omit(tmp_genes))
}

celltye_time_disease_genes_1[['AT1_Mfrp.103:P1']]

celltye_time_disease_genes_1[['AT1_Nek10.112:P0']]

#options(jupyter.plot_mimetypes = c("image/png", "image/jpeg", "image/svg+xml"))
xx = celltye_time_disease_genes_1[['AT1_Nek10.112:P0']]$unique
xx$`Adenosquamous cell lung cancer` = c(xx$`Adenosquamous cell lung cancer`, 'EGFR')
xx$`Non-small cell lung cancer stage I` = c(xx$`Non-small cell lung cancer stage I`, 'EGFR')
xx$`stage, non-small cell lung cancer` = c(xx$`stage, non-small cell lung cancer`, 'EGFR')
xx$`Progression of non-small cell lung cancer` = c(xx$`Progression of non-small cell lung cancer`, 'EGFR')
xx$`Congenital hypoplasia of lung` = c(xx$`Congenital hypoplasia of lung`, 'EGFR')

lt = list_to_matrix(xx)
m1 = make_comb_mat(lt)
pdf(glue('{output_dir}/AT1_nek10_geneset.pdf'), width = 12,height = 4)
draw(UpSet(m1))
dev.off()

aa = intersect(xx$`Congenital hypoplasia of lung`, xx$`Non-small cell lung cancer stage I`)
intersect(aa, xx$`stage, non-small cell lung cancer`)

aa = intersect(xx$`Congenital hypoplasia of lung`, xx$`stage, non-small cell lung cancer`)
intersect(aa, xx$`Adenosquamous cell lung cancer`)

aa = intersect(xx$`Non-small cell lung cancer stage I`, xx$`stage, non-small cell lung cancer`)
intersect(aa, xx$`Progression of non-small cell lung cancer`)

aa = intersect(xx$`Non-small cell lung cancer stage I`, xx$`stage, non-small cell lung cancer`)
intersect(aa, xx$`Adenosquamous cell lung cancer`)

intersect(xx$`Congenital hypoplasia of lung`, xx$`Non-small cell lung cancer stage I`)

intersect(xx$`Congenital hypoplasia of lung`, xx$`stage, non-small cell lung cancer`)

intersect(xx$`Non-small cell lung cancer stage I`, xx$`stage, non-small cell lung cancer`)

intersect(xx$`Non-small cell lung cancer stage I`, xx$`Progression of non-small cell lung cancer`)

xx


tmp_gene=c('Egfr')
cur_cre = c('chr11-16769901-16770402','chr18-60714271-60714772')
#cur_cre = tmp_gene_cre2 %>% filter(cluster%in%tmp_cluster) %>% pull(gene)

tmp_srt = subset(srt, Celltype_round4=='AT1_Nek10')

tmp_srt$Time = factor(tmp_srt$Time, time_levels)

tmp_srt@assays$CRE@meta.features[,'peak_group'] = 'false'
tmp_srt@assays$CRE@meta.features[cur_cre,'peak_group'] = 'true'

options(repr.plot.width=10, repr.plot.height=10)

a=CoveragePlot(tmp_srt, 
               region = tmp_gene,#c('chr11-16769901-16770402'),#'chr7-144963831-144964332',
               group.by = 'Time',
               ymax='q50', 
               peaks.group.by='peak_group',
               window=200, 
               extend.upstream = 0,
               extend.downstream = 0,
              )&  scale_fill_manual(values = time_color2)
a
ggsave(glue('{output_dir}/AT1_nek10_Egfr_CoveragePlot.pdf'), a,  
       width=100, height = 120, units='mm', dpi=600, bg='transparent')

options(repr.plot.width=10, repr.plot.height=10)

a=CoveragePlot(tmp_srt, 
               region = 'Ndst1',#c('chr11-16769901-16770402'),#'chr7-144963831-144964332',
               group.by = 'Time',
               #ymax='q50', 
               peaks.group.by='peak_group',
               window=200, 
               extend.upstream = 0,
               extend.downstream = 3000,
               cols=time_color2
              )&  scale_fill_manual(values = time_color2)
a
ggsave(glue('{output_dir}/AT1_nek10_Ndst1_CoveragePlot.pdf'), a,  
       width=100, height = 120, units='mm', dpi=600, bg='transparent')

options(repr.plot.width=10, repr.plot.height=10)
a=CoveragePlot(tmp_srt, 
               region = cur_cre,#c('chr11-16769901-16770402'),#'chr7-144963831-144964332',
               group.by = 'Time',
               ymax='q95', 
               window=100, 
               extend.upstream = 250,
               extend.downstream = 250
              )&  scale_fill_manual(values = time_color2)
a
ggsave(glue('{output_dir}/AT1_nek10_Egfr_Ndst1_cre_CoveragePlot.pdf'), a,  
       width=100, height = 120, units='mm', dpi=600, bg='transparent')

p0_srt = srt

p0_srt$new_name = new_celltype[as.vector(p0_srt$Celltype_round4),'new_name']

umap = read.csv(glue('./placeholder_analysis/round_cluster02/round1/Endothelium//cell_meta_ann.csv'))
rownames(umap) = umap$X

same_bc = intersect(umap$X,colnames(p0_srt))

head(unlist(umap[same_bc, c('UMAP_1')]))
p0_srt@reductions$umap@cell.embeddings[same_bc, 1] = unlist(umap[same_bc, c('UMAP_1')])
p0_srt@reductions$umap@cell.embeddings[same_bc, 2] = unlist(umap[same_bc, c('UMAP_2')])

tmp_umap = p0_srt@meta.data
tmp_umap[,c('UMAP_1', 'UMAP_2')] = p0_srt@reductions$umap@cell.embeddings

a = DimPlot(p0_srt, group.by = 'new_name', cols = tmp_curr_celltype_color2, raster = T)
a
ggsave(glue('{output_dir}/Epi_umap_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/Epi_umap.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

a = myDimPlot(tmp_umap, groupby='new_name', label=FALSE, group_color=tmp_curr_celltype_color2, point_size=0.4)
a
ggsave(glue('{output_dir}/Epi_umap_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/Epi_umap.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

# ggsave(glue('{output_dir}/Epi_umap.pdf'), a,  
#        width=140, height = 100, units='mm', dpi=600, bg='transparent')

a = DimPlot(p0_srt, group.by = 'Time', cols = time_color2, raster = T)
a
ggsave(glue('{output_dir}/Epi_umap(Time)_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/Epi_umap(Time).pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

tmp_lung_res0 = read.table(glue('{output_dir}/Lung_enrich.txt'), header=T, sep='\t')

head(tmp_lung_res0)

xlim_manual = c(min(p0_srt@reductions$umap@cell.embeddings[,1]),
                 max(p0_srt@reductions$umap@cell.embeddings[,1]))
xlim_manual

ylim_manual =c(min(p0_srt@reductions$umap@cell.embeddings[,2]),
                 max(p0_srt@reductions$umap@cell.embeddings[,2]))
ylim_manual

p0_srt_sub = subset(p0_srt, Time=='P0')
tmp_umap_df = p0_srt_sub@meta.data
tmp_umap_df = tmp_umap_df[, c('new_name'), drop=F]
tmp_umap_df[, c('UMAP_1', 'UMAP_2')] = p0_srt_sub@reductions$umap@cell.embeddings

tmp_lung_res0_P0 = tmp_lung_res0 %>% filter(Time=='P0')

for(i in unique(tmp_lung_res0_P0$Description)){
    tmp_genes = tmp_lung_res0_P0 %>% filter(Description==i) %>% pull(Celltype)
    tmp_umap_df[, i] = 0
    tmp_umap_df[tmp_umap_df$new_name%in%tmp_genes, i] = 1
}

library(mascarade)

maskTable <- generateMask(dims=p0_srt@reductions$umap@cell.embeddings, 
                          clusters=p0_srt$new_name,
                          minDensity = 15,smoothSigma = 0.1)



options(repr.plot.width=24, repr.plot.height=6)
a =melt(tmp_umap_df, id.vars = c('UMAP_1', 'UMAP_2', 'new_name'))  %>%
  filter(value>0) %>%
  # filter(new_name!='Airway_Spock3.4') %>%
  #filter(variable=='Adenosquamous cell lung cancer') %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) + 
  #stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F, adjust = 10) + 
  geom_point(aes(color=new_name),size = .5, shape=16) + 
  facet_wrap(~variable, nrow = 1) + 
  scale_color_manual(values = tmp_curr_celltype_color2) +
  geom_path(data=maskTable, aes(x=umap_1,y=umap_2,group=group),linewidth=0.2,linetype = 2, colour ="gray") +
  #scale_fill_viridis(option="magma") + 
  #galaxyTheme_black() + 
  #lims(x=xlim_manual, y=ylim_manual) +
  theme_void() +
  theme(strip.text = element_text(size=6),  # remove facet title
        legend.position = "none" ,
        panel.border = element_rect(fill=NA, color='black'))   # remove legend
a
ggsave(glue('{output_dir}/P0_nek10_density2.pdf'), a,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

p0_srt_sub = subset(p0_srt, Time=='P7')
tmp_umap_df = p0_srt_sub@meta.data
tmp_umap_df = tmp_umap_df[, c('new_name'), drop=F]
tmp_umap_df[, c('UMAP_1', 'UMAP_2')] = p0_srt_sub@reductions$umap@cell.embeddings

tmp_lung_res0_P0 = tmp_lung_res0 %>% filter(Time=='P7')
for(i in unique(tmp_lung_res0_P0$Description)){
    tmp_genes = tmp_lung_res0_P0 %>% filter(Description==i) %>% pull(Celltype)
    #print(tmp_genes)
    if(tmp_genes[1]=='AT1_Spock2.116') {
        tmp_genes='AT1.116'
    }
    tmp_umap_df[, i] = 0
    tmp_umap_df[tmp_umap_df$new_name%in%tmp_genes, i] = 1
}



options(repr.plot.width=24, repr.plot.height=6)
a =melt(tmp_umap_df, id.vars = c('UMAP_1', 'UMAP_2', 'new_name'))  %>%
  filter(value>0) %>%
  # filter(new_name!='Airway_Spock3.4') %>%
  #filter(variable=='Adenosquamous cell lung cancer') %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) + 
  #stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F, adjust = 10) + 
  geom_point(aes(color=new_name),size = 0.5, shape=16) + 
  facet_wrap(~variable, nrow = 1) + 
  scale_color_manual(values = tmp_curr_celltype_color2) +
  geom_path(data=maskTable, aes(x=umap_1,y=umap_2,group=group),linewidth=0.2,linetype = 2, colour ="gray") +
  #scale_fill_viridis(option="magma") + 
  #galaxyTheme_black() + 
  #lims(x=xlim_manual, y=ylim_manual) +
  theme_void() +
  theme(strip.text = element_text(size=6),  # remove facet title
        legend.position = "none" ,
        panel.border = element_rect(fill=NA, color='black'))   # remove legend
a
ggsave(glue('{output_dir}/P7_density2.pdf'), a,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

mm10_hg38_homoGene = read.csv('./placeholder_project/genome_annotation_file/mm10_hg38_homoGene_ensembl.csv')
#background_gene = unique(na.omit(mm10_hg38_homoGene$Human.gene.name))
head(mm10_hg38_homoGene)

p0_srt_sub = subset(p0_srt, Time=='P0')

tmp_all_gene = unname(unlist(celltye_time_disease_genes_1[['AT1_Nek10.112:P0']]))
tmp_all_gene = sapply(tmp_all_gene, function(x){
    paste0(substr(x,1,1), tolower(substr(x,2,nchar(x))))
    }
)
tmp_all_gene = unname(tmp_all_gene)

gene_activity = GeneActivity(p0_srt,features = tmp_all_gene)

p0_srt[['RNA']] = CreateAssayObject(counts=gene_activity)
DefaultAssay(p0_srt) = 'RNA'
p0_srt = NormalizeData(p0_srt, assay = 'RNA', normalization.method = 'LogNormalize',
                    scale.factor = median(p0_srt$nCount_RNA))

umap = read.csv(glue('./placeholder_analysis/round_cluster02/round1/Endothelium//cell_meta_ann.csv'))
rownames(umap) = umap$X

same_bc = intersect(umap$X,colnames(p0_srt))

head(unlist(umap[same_bc, c('UMAP_1')]))
p0_srt@reductions$umap@cell.embeddings[same_bc, 1] = unlist(umap[same_bc, c('UMAP_1')])
p0_srt@reductions$umap@cell.embeddings[same_bc, 2] = unlist(umap[same_bc, c('UMAP_2')])

density_colormap1 =c(
  "#313695",  # blue
  "#74ADD1",  # light blue
  "#FDAE61",  # orange
  "#D73027",  # red
  "#A50026"   # deep red
)
density_colormap2 =c(
  "#440154",  # deep purple / low value
  "#3B528B",  # blue
  "#21908C",  # cyan
  "#5DC863",  # green
  "#FDE725"   # high value / yellow
)

options(repr.plot.width=15, repr.plot.height=5)
p0_srt_sub = p0_srt#subset(p0_srt, Time=='P0')
a=DimPlot(p0_srt_sub, group.by = 'new_name', raster = F)+
    #theme(legend.position = 'none')+
    labs(title='')+
    scale_color_manual(values = tmp_curr_celltype_color2)
b=scCustomize::Plot_Density_Custom(seurat_object = p0_srt_sub, 
    features = c("Egfr"), viridis_palette= "viridis", pt.size=0.5)&scale_color_gradientn(colours = density_colormap2)
c=scCustomize::Plot_Density_Custom(seurat_object = p0_srt_sub, 
    features = c('Ndst1'), viridis_palette= "viridis", pt.size=0.5)&scale_color_gradientn(colours = density_colormap1)
#b=FeaturePlot(p0_srt, features = 'Egfr', cols = c('lightgray', 'red'))
a+b+c

ggsave(glue('{output_dir}/AT1_egfr_ndst_umap.pdf'), a+b+c,  
       width=320, height = 100, units='mm', dpi=600, bg='transparent')



