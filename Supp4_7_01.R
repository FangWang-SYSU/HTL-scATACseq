setwd('./placeholder_project/code/')
source('./placeholder_project/code/jupyter_R/myFun.R')
source('./placeholder_project/code/jupyter_R/markerlist.R')
sample_info = read.csv('./placeholder_project//sample_information/sample_information.csv',row.names = 1)
#base_path='./placeholder_analysis/round_cluster02/round0/'

myHeatmapPeaks <- function (marker.peaks, obj = NULL, group = "seurat_clusters", 
    top = 10, assay = "ATAC", min_log2FC = 0.25, max_pval = 0.05, 
    color = circlize::colorRamp2(c(-2, 0, 2), c("#045a8d", "white", 
        "#a50f15")), row_fontsize = 3, column_fontsize = 3, cell_meta = NULL, 
    mat = NULL, filter_pseudogene = FALSE, lwd = 0.05) 
{
    if (!is.null(obj)) {
        cell_meta = obj@meta.data
        mat = obj@assays[[assay]]$data
    }
    if (filter_pseudogene) {
        marker.peaks = filter_pseudoGene(marker.peaks)
    }
    top_peaks = marker.peaks %>% dplyr::filter(avg_log2FC >= 
        min_log2FC & p_val_adj <= max_pval) %>% group_by(cluster) %>% 
        top_n(top, avg_log2FC) %>% arrange(cluster, -avg_log2FC)
    cell_meta = cell_meta[order(cell_meta[, group]), ]
    top_mat = mat[top_peaks$gene, rownames(cell_meta)]
    top_mat = t(apply(top_mat, 1, scale))
    ht = Heatmap(as.matrix(top_mat), col = color, cluster_rows = F, 
        cluster_columns = F, show_row_names = F, show_column_names = F, 
        column_split = cell_meta[, group], row_split = top_peaks$cluster, 
        use_raster = TRUE, column_gap = unit(0, "mm"), row_gap = unit(0, 
            "mm"), border = TRUE, border_gp = gpar(col = "black", 
            lwd = lwd), row_title_gp = gpar(fontsize = row_fontsize), 
        column_title_gp = gpar(fontsize = column_fontsize))
    return(ht)
}

pseudoObj_peaks <- readRDS(glue('./placeholder_analysis/round_cluster02/round0/pseudoObj_peaks.rds'))
pseudoObj_cells <- readRDS(glue('./placeholder_analysis/round_cluster02/round0/pseudoObj_cells.rds'))

merge_path='./placeholder_analysis/round_cluster02/merge'


umap = read.csv(glue('{merge_path}/cell_meta.csv'))
rownames(umap) = umap$X
umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)

new_celltype = read.csv('./placeholder_analysis/round_cluster02/merge/celltype_modify.csv')
rownames(new_celltype) = new_celltype$Celltype_round4
new_celltype$modify_name = new_celltype$Celltype_round4_new
new_celltype$new_name = paste0(new_celltype$modify_name,'.', new_celltype$cluster_num)

pseudoObj_cells = subset(pseudoObj_cells, cells = umap$X)
pseudoObj_cells

pseudoObj_cells@meta.data[umap$X, colnames(umap)] = umap
Idents(pseudoObj_cells) = pseudoObj_cells$Celltype_round1
pseudoObj_cells@meta.data[,'nCount_ATAC'] = pseudoObj_cells$n_fragment
pseudoObj_cells@meta.data[,'nFeature_ATAC'] = pseudoObj_cells$n_peaks

pseudoObj_cells@meta.data$Celltype_round1_new = new_celltype[pseudoObj_cells$Celltype_round4, 'new_name']

pseudoObj_cells[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(pseudoObj_cells@meta.data[,c('UMAP_1','UMAP_2')]),
  key = "UMAP_",
  assay = DefaultAssay(pseudoObj_cells)
)
pseudoObj_cells[["tsne"]] <- CreateDimReducObject(
  embeddings = as.matrix(pseudoObj_cells@meta.data[,c('TSNE_1','TSNE_2')]),
  key = "UMAP_",
  assay = DefaultAssay(pseudoObj_cells)
)

a=DimPlot(pseudoObj_cells, reduction = 'umap', group.by = 'Celltype_round1')
b=DimPlot(pseudoObj_cells, reduction = 'tsne', group.by = 'Celltype_round1')
a+b

merge_path='./placeholder_analysis/round_cluster02_gene/'
mat <- Matrix::readMM(glue("{merge_path}/gene_X_sparse.mtx"))
row_names <- read.csv(glue("{merge_path}/gene_X_sparse_index.csv"), header = T)[,1]
col_names <- read.csv(glue("{merge_path}/gene_X_sparse_columns.csv"), header = T)[,1]
row_names = gsub(':','-',row_names)
rownames(mat) <- row_names
colnames(mat) <- col_names

gene_activate_srt = CreateSeuratObject(mat, min.cells = 0, min.features = 0, project='GeneActivity')

umap2=umap[umap$X%in%colnames(gene_activate_srt),]
gene_activate_srt@meta.data[umap2$X, colnames(umap2)] = umap2

gene_activate_srt[['umap']]=CreateDimReducObject(embeddings = as.matrix(gene_activate_srt@meta.data[,c('UMAP_1', 'UMAP_2')]),
                                                 key='UMAP_')
gene_activate_srt[['tsne']]=CreateDimReducObject(embeddings = as.matrix(gene_activate_srt@meta.data[,c('TSNE_1', 'TSNE_2')]),
                                                 key='TSNE_')

a=DimPlot(gene_activate_srt, reduction = 'umap', group.by = 'Celltype_round1')
b=DimPlot(gene_activate_srt, reduction = 'tsne', group.by = 'Celltype_round1')
a+b

saveRDS(gene_activate_srt, './placeholder_analysis/round_cluster02/merge/gene_activate_srt.rds')
saveRDS(pseudoObj_cells, './placeholder_analysis/round_cluster02/merge/pseudoObj_srt.rds')

pseudoObj_peaks <- readRDS(glue('./placeholder_analysis/round_cluster02/round0/pseudoObj_peaks.rds'))
pseudoObj_cells <- readRDS(glue('./placeholder_analysis/round_cluster02/merge/pseudoObj_srt.rds'))

gene_activate_srt = readRDS('./placeholder_analysis/round_cluster02/merge/gene_activate_srt.rds')

drop_cells = rownames(gene_activate_srt@meta.data[is.na(gene_activate_srt$Celltype_round1),])

gene_activate_srt = subset(gene_activate_srt, cells=drop_cells, invert=TRUE)

gene_activate_srt

# Load differential peak matrix
tmp_dir='./placeholder_analysis/round_cluster02/merge/'
peak_mat <- Matrix::readMM(glue("{tmp_dir}/peak_cluster_markers_sparse.mtx"))
row_names <- read.csv(glue("{tmp_dir}/peak_cluster_markers_sparse_index.csv"), header = T)[,1]
col_names <- read.csv(glue("{tmp_dir}/peak_cluster_markers_sparse_columns.csv"), header = T)[,1]
row_names = gsub(':','-',row_names)
rownames(peak_mat) <- row_names
colnames(peak_mat) <- col_names

output_dir='./placeholder_output/raw_figure_output/Supp01/'

round_num='Round1'
type_name='All'
#output_dir='../pycode/round_cluster_02_after/celltype_marker/'
marker = read.csv(glue('../pycode/round_cluster_02_after/celltype_marker//{round_num}_All_marker_know.csv'))

colnames(marker)[c(2,5,6,7)] = c('cluster', 'avg_log2FC', 'p_val','p_val_adj')
marker$gene_name=marker$gene

head(marker)

marker2 = marker%>%filter(FindDiffMethod=='gene')
marker2 = marker2[marker2$gene%in%rownames(gene_activate_srt),]

marker2$avg_log2FC = marker2$scores
marker2=marker2[!duplicated(marker2$gene_name),]

top_peaks = marker2 %>% 
    dplyr::filter(avg_log2FC >=0.25&p_val_adj <=0.05) %>%
    group_by(cluster) %>% 
    top_n(10,avg_log2FC) %>%
    arrange(cluster, -avg_log2FC)

dim(gene_activate_srt)
dim(gene_activate_srt@meta.data)

options(repr.plot.width=6, repr.plot.height=10)
a=DotPlot(gene_activate_srt, 
        features = unique(top_peaks$gene),
        group.by='Celltype_round1')+
    scale_color_gradientn(colors = c('#565ca2','#cce5a6','#f7d84f','#972245'))+
    scale_size_continuous(range=c(0,5), limits = c(0,50))+
    #coord_flip()+
    theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1,vjust=1),
        axis.title = element_blank(),
        strip.background = element_blank(),
        legend.key.height = unit(3,'mm'),
        legend.key.size =  unit(2,'mm'))
a

ggsave(glue('{output_dir}/Round1_dotplot.pdf'), a,
       width=12, height=3, units='in', dpi=600, bg='transparent')

genes_id = top_peaks[top_peaks$cluster=='Immune','gene',drop=T]
print(genes_id)
ncol = length(genes_id)
suppressMessages(suppressWarnings({
    plot_res =myCoveragePlotMultiple(pseudoObj_cells,pseudoObj_peaks,region = genes_id,
                                 ncol=ncol,
                                 ymax='q95', window=500, heights = c(20,1,1),
                                 extend.upstream = 2000, 
                                 extend.downstream = 2000)
    a=myCoverageModify(plot_res, ncol=ncol,region = genes_id,
                 group_color=NULL, genes_color=c('black', 'black'), peaks_color='red',
                 line_color_with_group=F)
}))
options(repr.plot.width=2*ncol, repr.plot.height=4)
a

ggsave(glue('{output_dir}/Round1_immune_coverageplot.pdf'), a,
       width=12, height=3, units='in', dpi=600, bg='transparent')

genes_id = top_peaks[top_peaks$cluster=='Endothelium','gene',drop=T]
print(genes_id)
ncol = length(genes_id)
suppressMessages(suppressWarnings({
    plot_res =myCoveragePlotMultiple(pseudoObj_cells,pseudoObj_peaks,region = genes_id,
                                 ncol=ncol,
                                 ymax='q95', window=500, heights = c(20,1,1),
                                 extend.upstream = 2000, 
                                 extend.downstream = 2000)
    a=myCoverageModify(plot_res, ncol=ncol,region = genes_id,
                 group_color=NULL, genes_color=c('black', 'black'), peaks_color='red',
                 line_color_with_group=F)
}))
options(repr.plot.width=2*ncol, repr.plot.height=4)
a
ggsave(glue('{output_dir}/Round1_endo_coverageplot.pdf'), a,
       width=12, height=3, units='in', dpi=600, bg='transparent')

genes_id = top_peaks[top_peaks$cluster=='Epithelium','gene',drop=T]
print(genes_id)
ncol = length(genes_id)
suppressMessages(suppressWarnings({
    plot_res =myCoveragePlotMultiple(pseudoObj_cells,pseudoObj_peaks,region = genes_id,
                                 ncol=ncol,
                                 ymax='q95', window=500, heights = c(20,1,1),
                                 extend.upstream = 2000, 
                                 extend.downstream = 2000)
    a=myCoverageModify(plot_res, ncol=ncol,region = genes_id,
                 group_color=NULL, genes_color=c('black', 'black'), peaks_color='red',
                 line_color_with_group=F)
}))
options(repr.plot.width=2*ncol, repr.plot.height=4)
a
ggsave(glue('{output_dir}/Round1_epi_coverageplot.pdf'), a,
       width=12, height=3, units='in', dpi=600, bg='transparent')

genes_id = top_peaks[top_peaks$cluster=='Stroma','gene',drop=T]
print(genes_id)
ncol = length(genes_id)
suppressMessages(suppressWarnings({
    plot_res =myCoveragePlotMultiple(pseudoObj_cells,pseudoObj_peaks,region = genes_id,
                                 ncol=ncol,
                                 ymax='q95', window=500, heights = c(20,1,1),
                                 extend.upstream = 2000, 
                                 extend.downstream = 2000)
    a=myCoverageModify(plot_res, ncol=ncol,region = genes_id,
                 group_color=NULL, genes_color=c('black', 'black'), peaks_color='red',
                 line_color_with_group=F)
}))
options(repr.plot.width=2*ncol, repr.plot.height=4)
a
ggsave(glue('{output_dir}/Round1_Stromal_coverageplot.pdf'), a,
       width=12, height=3, units='in', dpi=600, bg='transparent')

head(marker %>% 
    filter(FindDiffMethod=='peak'))

tmp_markers = marker %>% 
    filter(FindDiffMethod=='peak') %>%
    group_by(cluster) %>% 
    top_n(100,scores) %>%
    mutate(names=gsub(':','-',names)) %>%
    mutate(gene=names)

same_peaks = intersect(tmp_markers$names, rownames(peak_mat))

tmp_mat2=peak_mat[same_peaks, colnames(gene_activate_srt)]

tmp_markers = tmp_markers %>%
    filter(names%in%same_peaks)

x=myHeatmapPeaks(mat=tmp_mat2, 
                 cell_meta=gene_activate_srt@meta.data, 
                 marker.peaks=tmp_markers,
                 group='Celltype_round1',
                 top=30,assay='ATAC',
                 min_log2FC=0.25,max_pval=0.05,
                 color=circlize::colorRamp2(c(-2,0,2),c("#045a8d","white","#a50f15")),
                 row_fontsize=12, column_fontsize=12)
options(repr.plot.width=16, repr.plot.height=6)
draw(x)

output_dir

pdf(glue('{output_dir}/Round1_peak_heatmap.pdf'), width = 10,height = 6)
draw(x)
dev.off()

new_celltype = read.csv('./placeholder_analysis/round_cluster02/merge/celltype_modify.csv')
rownames(new_celltype) = new_celltype$Celltype_round4
new_celltype$modify_name = new_celltype$Celltype_round4_new
new_celltype$new_name = paste0(new_celltype$modify_name,'.', new_celltype$cluster_num)

tmp_srt_gene = subset(gene_activate_srt, Celltype_round1=='Immune')
tmp_srt_peak = subset(pseudoObj_cells, Celltype_round1=='Immune')

options(repr.plot.width=16, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'Celltype_round1')
b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round1')
a+b

tmp_umap = read.csv(glue('./placeholder_analysis/round_cluster02/round1/Immune/cell_meta.csv'))
rownames(tmp_umap) = tmp_umap$X
tmp_same_bc = intersect(rownames(tmp_umap), colnames(tmp_srt_peak))

tmp_srt_peak = subset(tmp_srt_peak, cells=tmp_same_bc)
tmp_umap = tmp_umap[tmp_same_bc,]

tmp_srt_peak@meta.data$UMAP_1 = NULL
tmp_srt_peak@meta.data$UMAP_2 = NULL

tmp_srt_peak[['umap']]=NULL
tmp_srt_peak[['tsne']]=NULL

tmp_srt_peak[['umap']] <- CreateDimReducObject(
  embeddings = as.matrix(tmp_umap[colnames(tmp_srt_peak),c('UMAP_1','UMAP_2')]),
  key = "UMAP_",
  assay = DefaultAssay(tmp_srt_peak)
)
tmp_srt_peak[['tsne']] <- CreateDimReducObject(
  embeddings = as.matrix(tmp_umap[colnames(tmp_srt_peak),c('TSNE_1','TSNE_2')]),
  key = "TSNE_",
  assay = DefaultAssay(tmp_srt_peak)
)

round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')
new_celltype$color = round4_color[new_celltype$Celltype_round4]
write.table(new_celltype, './placeholder_analysis/round_cluster02/merge/celltype_modify_color.csv', sep=',', row.names = F, quote = F)

round4_num = readRDS('./placeholder_analysis/round_cluster02/merge/round4_cluster_num.rds')
round4_color2 = round4_color
names(round4_color2) = round4_num[names(round4_color2)]

names(round4_color) = new_celltype[names(round4_color), 'new_name']

tmp_srt_peak$new_name = as.vector(new_celltype[as.vector(tmp_srt_peak$Celltype_round4), 'new_name'])
tmp_srt_peak$name_num = as.vector(new_celltype[as.vector(tmp_srt_peak$Celltype_round4), 'cluster_num'])
tmp_srt_gene$new_name = as.vector(new_celltype[as.vector(tmp_srt_gene$Celltype_round4), 'new_name'])
tmp_srt_gene$name_num = as.vector(new_celltype[as.vector(tmp_srt_gene$Celltype_round4), 'cluster_num'])

round4_color2['2']

options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'new_name', raster = T)+
    scale_colour_manual(values = round4_color)
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a


ggsave(glue('{output_dir}/Immune_round4_umap.pdf'), 
       a,
       width=14, height=6, units='in', dpi=600, bg='transparent')

options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'name_num',
     raster = T, label = T, label.box=T,repel=T,label.color = "white", #alpha=0.5,
     cols=round4_color2)+NoLegend()
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a


ggsave(glue('{output_dir}/Immune_round4_umap(label).pdf'), 
       a,
       width=6, height=6, units='in', dpi=600, bg='transparent')

tmp_srt_gene <- NormalizeData(tmp_srt_gene, normalization.method = "LogNormalize", scale.factor = 10000)
tmp_srt_gene <- FindVariableFeatures(tmp_srt_gene, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tmp_srt_gene)
tmp_srt_gene <- ScaleData(tmp_srt_gene, features = all.genes)
tmp_srt_gene

Idents(tmp_srt_gene) = tmp_srt_gene$Celltype_round4
gene_markers = FindAllMarkers(tmp_srt_gene, max.cells.per.ident = 300)

output_dir

write.table(gene_markers, 
            file=glue('{output_dir}/Immune_round4_gene_markers.txt'),
            sep='\t', quote=F, row.names=F)

gene_markers = read.table(glue('{output_dir}/Immune_round4_gene_markers.txt'), header = T)
gene_markers$cluster = new_celltype[gene_markers$cluster, 'cluster_num']

#tmp_srt_gene$Celltype_round4 = gsub('Macorphage_','Macrophage_', tmp_srt_gene$Celltype_round4)
#gene_markers$cluster = gsub('Macorphage_','Macrophage_', gene_markers$cluster)


top_peaks = gene_markers %>% 
    dplyr::filter(avg_log2FC >=0.25&p_val_adj <=0.05) %>%
    group_by(cluster) %>% 
    top_n(4,avg_log2FC) %>%
    arrange(cluster, -avg_log2FC)
cluster_order = sort(as.vector(unique(tmp_srt_gene$name_num)))
top_peaks$cluster = factor(top_peaks$cluster, levels=cluster_order)
top_peaks = top_peaks %>% arrange(cluster, -avg_log2FC)
tmp_srt_gene$name_num = factor(tmp_srt_gene$name_num, levels=cluster_order)
options(repr.plot.width=12, repr.plot.height=18)
a=DotPlot(tmp_srt_gene, 
        features = unique(top_peaks$gene),
        group.by='name_num')+
    scale_color_gradientn(colors = c('#565ca2','#cce5a6','#f7d84f','#972245'))+
    scale_size_continuous(range=c(0,5), limits = c(0,50))+
    #coord_flip()+
    theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1,vjust=1),
        axis.title = element_blank(),
        strip.background = element_blank(),
        legend.key.height = unit(3,'mm'),
        legend.key.size =  unit(2,'mm'))+
    coord_flip()
a

ggsave(glue('{output_dir}/Immune_round4_gene_dot.pdf'), 
       a,
       width=10, height=26, units='in', dpi=600, bg='transparent')



tmp_srt_gene = subset(gene_activate_srt, Celltype_round1=='Epithelium')
tmp_srt_peak = subset(pseudoObj_cells, Celltype_round1=='Epithelium')

options(repr.plot.width=16, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'Celltype_round1')
b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round1')
a+b

tmp_umap = read.csv(glue('./placeholder_analysis/round_cluster02/round1/Endothelium/cell_meta.csv'))
rownames(tmp_umap) = tmp_umap$X

tmp_same_bc = intersect(rownames(tmp_umap), colnames(tmp_srt_peak))

length(tmp_same_bc)

tmp_srt_peak = subset(tmp_srt_peak, cells=tmp_same_bc)
tmp_umap = tmp_umap[tmp_same_bc,]

tmp_srt_peak@meta.data$UMAP_1 = NULL
tmp_srt_peak@meta.data$UMAP_2 = NULL

tmp_srt_peak[['umap']] <- CreateDimReducObject(
  embeddings = as.matrix(tmp_umap[colnames(tmp_srt_peak),c('UMAP_1','UMAP_2')]),
  key = "UMAP_",
  assay = DefaultAssay(tmp_srt_peak)
)
tmp_srt_peak[['tsne']] <- CreateDimReducObject(
  embeddings = as.matrix(tmp_umap[colnames(tmp_srt_peak),c('TSNE_1','TSNE_2')]),
  key = "TSNE_",
  assay = DefaultAssay(tmp_srt_peak)
)

round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')
round4_num = readRDS('./placeholder_analysis/round_cluster02/merge/round4_cluster_num.rds')
round4_color2 = round4_color
names(round4_color2) = round4_num[names(round4_color2)]

names(round4_color) = new_celltype[names(round4_color), 'new_name']

tmp_srt_peak$new_name = as.vector(new_celltype[as.vector(tmp_srt_peak$Celltype_round4), 'new_name'])
tmp_srt_peak$name_num = as.vector(new_celltype[as.vector(tmp_srt_peak$Celltype_round4), 'cluster_num'])
tmp_srt_gene$new_name = as.vector(new_celltype[as.vector(tmp_srt_gene$Celltype_round4), 'new_name'])
tmp_srt_gene$name_num = as.vector(new_celltype[as.vector(tmp_srt_gene$Celltype_round4), 'cluster_num'])

options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'Celltype_round4', raster = T)+
    scale_colour_manual(values = round4_color)
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a

ggsave(glue('{output_dir}/Epi_round4_umap.pdf'), 
       a,
       width=8, height=6, units='in', dpi=600, bg='transparent')

options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'new_name', raster = T)+
    scale_colour_manual(values = round4_color)
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a
ggsave(glue('{output_dir}/Epi_round4_umap.pdf'), 
       a,
       width=8, height=6, units='in', dpi=600, bg='transparent')


options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'name_num',
     raster = T, label = T, label.box=T,repel=T,label.color = "white", #alpha=0.5,
     cols=round4_color2)+NoLegend()
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a
ggsave(glue('{output_dir}/Epi_round4_umap(label).pdf'), 
       a,
       width=6, height=6, units='in', dpi=600, bg='transparent')

tmp_srt_gene <- NormalizeData(tmp_srt_gene, normalization.method = "LogNormalize", scale.factor = 10000)
tmp_srt_gene <- FindVariableFeatures(tmp_srt_gene, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tmp_srt_gene)
tmp_srt_gene <- ScaleData(tmp_srt_gene, features = all.genes)
tmp_srt_gene

Idents(tmp_srt_gene) = tmp_srt_gene$Celltype_round4
gene_markers = FindAllMarkers(tmp_srt_gene, max.cells.per.ident = 300)

write.table(gene_markers, 
            file=glue('{output_dir}/Epi_round4_gene_markers.txt'),
            sep='\t', quote=F, row.names=F)

gene_markers = read.table(glue('{output_dir}/Epi_round4_gene_markers.txt'), header = T)
gene_markers$cluster = new_celltype[gene_markers$cluster, 'cluster_num']

top_peaks = gene_markers %>% 
    dplyr::filter(avg_log2FC >=0.25&p_val_adj <=0.05) %>%
    group_by(cluster) %>% 
    top_n(4,avg_log2FC) %>%
    arrange(cluster, -avg_log2FC)
cluster_order = sort(as.vector(unique(tmp_srt_gene$name_num)))
top_peaks$cluster = factor(top_peaks$cluster, levels=cluster_order)
top_peaks = top_peaks %>% arrange(cluster, -avg_log2FC)
tmp_srt_gene$name_num = factor(tmp_srt_gene$name_num, levels=cluster_order)
options(repr.plot.width=12, repr.plot.height=18)
a=DotPlot(tmp_srt_gene, 
        features = unique(top_peaks$gene),
        group.by='name_num')+
    scale_color_gradientn(colors = c('#565ca2','#cce5a6','#f7d84f','#972245'))+
    scale_size_continuous(range=c(0,5), limits = c(0,50))+
    #coord_flip()+
    theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1,vjust=1),
        axis.title = element_blank(),
        strip.background = element_blank(),
        legend.key.height = unit(3,'mm'),
        legend.key.size =  unit(2,'mm'))+
    coord_flip()
a

ggsave(glue('{output_dir}/Epithelial_round4_gene_dot.pdf'), 
       a,
       width=7, height=12, units='in', dpi=600, bg='transparent')

tmp_srt_gene = subset(gene_activate_srt, Celltype_round1=='Stroma')
tmp_srt_peak = subset(pseudoObj_cells, Celltype_round1=='Stroma')

options(repr.plot.width=16, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'Celltype_round1')
b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round1')
a+b

tmp_umap = read.csv(glue('./placeholder_analysis/round_cluster02/round1/Stroma/cell_meta.csv'))
rownames(tmp_umap) = tmp_umap$X

tmp_same_bc = intersect(rownames(tmp_umap), colnames(tmp_srt_peak))

length(tmp_same_bc)

tmp_srt_peak = subset(tmp_srt_peak, cells=tmp_same_bc)
tmp_umap = tmp_umap[tmp_same_bc,]

tmp_srt_peak@meta.data$UMAP_1 = NULL
tmp_srt_peak@meta.data$UMAP_2 = NULL

tmp_srt_peak[['umap']] <- CreateDimReducObject(
  embeddings = as.matrix(tmp_umap[colnames(tmp_srt_peak),c('UMAP_1','UMAP_2')]),
  key = "UMAP_",
  assay = DefaultAssay(tmp_srt_peak)
)
tmp_srt_peak[['tsne']] <- CreateDimReducObject(
  embeddings = as.matrix(tmp_umap[colnames(tmp_srt_peak),c('TSNE_1','TSNE_2')]),
  key = "TSNE_",
  assay = DefaultAssay(tmp_srt_peak)
)

round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')
round4_num = readRDS('./placeholder_analysis/round_cluster02/merge/round4_cluster_num.rds')
round4_color2 = round4_color
names(round4_color2) = round4_num[names(round4_color2)]

names(round4_color) = new_celltype[names(round4_color), 'new_name']

tmp_srt_peak$new_name = as.vector(new_celltype[as.vector(tmp_srt_peak$Celltype_round4), 'new_name'])
tmp_srt_peak$name_num = as.vector(new_celltype[as.vector(tmp_srt_peak$Celltype_round4), 'cluster_num'])
tmp_srt_gene$new_name = as.vector(new_celltype[as.vector(tmp_srt_gene$Celltype_round4), 'new_name'])
tmp_srt_gene$name_num = as.vector(new_celltype[as.vector(tmp_srt_gene$Celltype_round4), 'cluster_num'])

options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'Celltype_round4', raster = T)+
    scale_colour_manual(values = round4_color)
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a

ggsave(glue('{output_dir}/Stroma_round4_umap.pdf'), 
       a,
       width=12, height=6, units='in', dpi=600, bg='transparent')

options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'new_name', raster = T)+
    scale_colour_manual(values = round4_color)
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a
ggsave(glue('{output_dir}/Stroma_round4_umap.pdf'), 
       a,
       width=8, height=6, units='in', dpi=600, bg='transparent')
options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'name_num',
     raster = T, label = T, label.box=T,repel=T,label.color = "white", #alpha=0.5,
     cols=round4_color2)+NoLegend()
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a
ggsave(glue('{output_dir}/Stroma_round4_umap(label).pdf'), 
       a,
       width=6, height=6, units='in', dpi=600, bg='transparent')

tmp_srt_gene <- NormalizeData(tmp_srt_gene, normalization.method = "LogNormalize", scale.factor = 10000)
tmp_srt_gene <- FindVariableFeatures(tmp_srt_gene, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tmp_srt_gene)
tmp_srt_gene <- ScaleData(tmp_srt_gene, features = all.genes)
tmp_srt_gene

Idents(tmp_srt_gene) = tmp_srt_gene$Celltype_round4
gene_markers = FindAllMarkers(tmp_srt_gene, max.cells.per.ident = 300)

write.table(gene_markers, 
            file=glue('{output_dir}/Stroma_round4_gene_markers.txt'),
            sep='\t', quote=F, row.names=F)

gene_markers = read.table(glue('{output_dir}/Stroma_round4_gene_markers.txt'), header = T)
gene_markers$cluster = new_celltype[gene_markers$cluster, 'cluster_num']

top_peaks = gene_markers %>% 
    dplyr::filter(avg_log2FC >=0.25&p_val_adj <=0.05) %>%
    group_by(cluster) %>% 
    top_n(4,avg_log2FC) %>%
    arrange(cluster, -avg_log2FC)
cluster_order = sort(as.vector(unique(tmp_srt_gene$name_num)))
top_peaks$cluster = factor(top_peaks$cluster, levels=cluster_order)
top_peaks = top_peaks %>% arrange(cluster, -avg_log2FC)
tmp_srt_gene$name_num = factor(tmp_srt_gene$name_num, levels=cluster_order)
options(repr.plot.width=12, repr.plot.height=18)
a=DotPlot(tmp_srt_gene, 
        features = unique(top_peaks$gene),
        group.by='name_num')+
    scale_color_gradientn(colors = c('#565ca2','#cce5a6','#f7d84f','#972245'))+
    scale_size_continuous(range=c(0,5), limits = c(0,50))+
    #coord_flip()+
    theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1,vjust=1),
        axis.title = element_blank(),
        strip.background = element_blank(),
        legend.key.height = unit(3,'mm'),
        legend.key.size =  unit(2,'mm'))+
    coord_flip()
a

ggsave(glue('{output_dir}/Stroma_round4_gene_dot.pdf'), 
       a,
       width=8, height=14, units='in', dpi=600, bg='transparent')



tmp_srt_gene = subset(gene_activate_srt, Celltype_round1=='Endothelium')
tmp_srt_peak = subset(pseudoObj_cells, Celltype_round1=='Endothelium')

options(repr.plot.width=16, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'Celltype_round1')
b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round1')
a+b

tmp_umap = read.csv(glue('./placeholder_analysis/round_cluster02/round1/Epithelium/cell_meta.csv'))
rownames(tmp_umap) = tmp_umap$X

tmp_same_bc = intersect(rownames(tmp_umap), colnames(tmp_srt_peak))

length(tmp_same_bc)

tmp_srt_peak = subset(tmp_srt_peak, cells=tmp_same_bc)
tmp_umap = tmp_umap[tmp_same_bc,]

tmp_srt_peak@meta.data$UMAP_1 = NULL
tmp_srt_peak@meta.data$UMAP_2 = NULL

tmp_srt_peak[['umap']]=NULL
tmp_srt_peak[['tsne']]=NULL

tmp_srt_peak[['umap']] <- CreateDimReducObject(
  embeddings = as.matrix(tmp_umap[colnames(tmp_srt_peak),c('UMAP_1','UMAP_2')]),
  key = "UMAP_",
  assay = DefaultAssay(tmp_srt_peak)
)
tmp_srt_peak[['tsne']] <- CreateDimReducObject(
  embeddings = as.matrix(tmp_umap[colnames(tmp_srt_peak),c('TSNE_1','TSNE_2')]),
  key = "TSNE_",
  assay = DefaultAssay(tmp_srt_peak)
)

round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')
round4_num = readRDS('./placeholder_analysis/round_cluster02/merge/round4_cluster_num.rds')
round4_color2 = round4_color
names(round4_color2) = round4_num[names(round4_color2)]

names(round4_color) = new_celltype[names(round4_color), 'new_name']

tmp_srt_peak$new_name = as.vector(new_celltype[as.vector(tmp_srt_peak$Celltype_round4), 'new_name'])
tmp_srt_peak$name_num = as.vector(new_celltype[as.vector(tmp_srt_peak$Celltype_round4), 'cluster_num'])
tmp_srt_gene$new_name = as.vector(new_celltype[as.vector(tmp_srt_gene$Celltype_round4), 'new_name'])
tmp_srt_gene$name_num = as.vector(new_celltype[as.vector(tmp_srt_gene$Celltype_round4), 'cluster_num'])

options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'new_name', raster = T)+
    scale_colour_manual(values = round4_color)
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a
ggsave(glue('{output_dir}/Endo_round4_umap.pdf'), 
       a,
       width=8, height=6, units='in', dpi=600, bg='transparent')
options(repr.plot.width=12, repr.plot.height=6)
a=DimPlot(tmp_srt_peak, reduction = 'umap', group.by = 'name_num',
     raster = T, label = T, label.box=T,repel=T,label.color = "white", #alpha=0.5,
     cols=round4_color2)+NoLegend()
#b=DimPlot(tmp_srt_peak, reduction = 'tsne', group.by = 'Celltype_round2', raster = T)
a
ggsave(glue('{output_dir}/Endo_round4_umap(label).pdf'), 
       a,
       width=6, height=6, units='in', dpi=600, bg='transparent')

tmp_srt_gene <- NormalizeData(tmp_srt_gene, normalization.method = "LogNormalize", scale.factor = 10000)
tmp_srt_gene <- FindVariableFeatures(tmp_srt_gene, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tmp_srt_gene)
tmp_srt_gene <- ScaleData(tmp_srt_gene, features = all.genes)
tmp_srt_gene

Idents(tmp_srt_gene) = tmp_srt_gene$Celltype_round4
gene_markers = FindAllMarkers(tmp_srt_gene, max.cells.per.ident = 300)

write.table(gene_markers, 
            file=glue('{output_dir}/Endo_round4_gene_markers.txt'),
            sep='\t', quote=F, row.names=F)

gene_markers = read.table(glue('{output_dir}/Endo_round4_gene_markers.txt'), header = T)
gene_markers$cluster = new_celltype[gene_markers$cluster, 'cluster_num']

top_peaks = gene_markers %>% 
    dplyr::filter(avg_log2FC >=0.25&p_val_adj <=0.05) %>%
    group_by(cluster) %>% 
    top_n(4,avg_log2FC) %>%
    arrange(cluster, -avg_log2FC)
cluster_order = sort(as.vector(unique(tmp_srt_gene$name_num)))
top_peaks$cluster = factor(top_peaks$cluster, levels=cluster_order)
top_peaks = top_peaks %>% arrange(cluster, -avg_log2FC)
tmp_srt_gene$name_num = factor(tmp_srt_gene$name_num, levels=cluster_order)
options(repr.plot.width=12, repr.plot.height=18)
a=DotPlot(tmp_srt_gene, 
        features = unique(top_peaks$gene),
        group.by='name_num')+
    scale_color_gradientn(colors = c('#565ca2','#cce5a6','#f7d84f','#972245'))+
    scale_size_continuous(range=c(0,5), limits = c(0,50))+
    #coord_flip()+
    theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1,vjust=1),
        axis.title = element_blank(),
        strip.background = element_blank(),
        legend.key.height = unit(3,'mm'),
        legend.key.size =  unit(2,'mm'))+
    coord_flip()
a

ggsave(glue('{output_dir}/Endo_round4_gene_dot.pdf'), 
       a,
       width=8, height=18, units='in', dpi=600, bg='transparent')

