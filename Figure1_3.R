setwd('./placeholder_project/code/')
source('./placeholder_project/code/jupyter_R/myFun.R')
source('./placeholder_project/code/jupyter_R/markerlist.R')
# Load required libraries
library(ggtree)
library(ape)
library(dplyr)
library(igraph)
library(forcats)
library(ggnewscale)
library(dplyr)
library(Seurat)
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
distribution_func<-function(meta, group_col, group_row, color=NULL, position='fill'){
  if(is.null(color)){
    color = c("#4dbbd5ff" ,"#E64b35ff", "#00a087ff" ,"#a65628", "#FF95A8FF","#BCBD22FF", "#fdbf6f", "#3c5488ff",
                       "#f39b7fff", "#b09c85ff", "#7876b1ff", "#377eb8",  "#4daf4a","#97A1A7FF" ,
                       "#984ea3" , "#ff7f00",  "#f781bf", "#b2df8a", "#5050FFFF", "#82581FFF" , "#E5614CFF",
                       "#F0E685FF", "#D595A7FF", "#CDDEB7FF","#612A79FF" ,"#AE1F63FF","#1B1919FF",
                       "#99CC00FF","#CC9900FF" ,"#9467BDFF", "#EFD500FF" , "#ADE2D0FF",pal_igv()(50))
  }
  a=meta %>%
    #dplyr::filter(sample_id!='0726_P1') %>%
    ggplot(aes_string(x=group_col))+
    geom_bar(aes_string(fill=group_row), stat='count', position = position)+
    scale_fill_manual(values = color)+
    theme_bw()+
    labs(y='Cells')+
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=60,hjust = 1,vjust=1))
  return(a)
}

distribution_ht_func<-function(meta, group_col, group_row, scale='none', show='count', size=6){
  orig_df = as.data.frame.matrix(table(meta[,group_row], meta[,group_col]))
  if(scale=='row'){
    df = orig_df / rowSums(orig_df)
  }else if(scale=='column'){
    df = t(t(orig_df) / colSums(orig_df))
  }else{
    df=orig_df
  }
  if(show=='count'){
    display_num=orig_df
  }else{
    display_num = round(df,2)
  }
  pheatmap::pheatmap(df, 
                     display_numbers = display_num,
                     border_color = 'white',
                     number_color = 'black',
                     fontsize = size,
                     color = grDevices::colorRampPalette(c('#2C7BB6','#ABD9E9','#FFFFFF','#FDAE61','#D7191C'))(100),
                     cluster_cols = F,cluster_rows = F)
}
distribution_Roe_func <- function(obj, group_by1, group_by2, show.sign=FALSE, 
                                  angle_col=45,order=NULL,obj_is_meta=FALSE,size=6){
  if(!obj_is_meta){
    metadata = obj@meta.data
  }else{
    metadata = obj
  }
  #metadata = obj@meta.data
  # Roe for each cluster
  observe_data1 = metadata %>% 
    dplyr::group_by_(group_by1, group_by2) %>% 
    dplyr::summarise(cell_num=n())
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
                     #color = grDevices::colorRampPalette(c(0,0.5,1,1.5,2),
                     #                                    c('#2C7BB6','#ABD9E9','#FFFFFF','#FDAE61','#D7191C'))(100),
                     color = c(colorRampPalette(c('#2C7BB6','#ABD9E9','#FFFFFF'))(length(bk)/2),
                               colorRampPalette(c('#FFFFFF','#FDAE61','#D7191C'))(length(bk)/2)),
                     breaks = bk,
                     display_numbers = annote.heatmap,
                     #cellwidth = 36, cellheight = 16, 
                     border_color = 'white',
                     number_color = 'black',
                     fontsize = size,
                     #border_color = '#ffffff',
                     angle_col=angle_col,
                     main = 'Cluster Dist (Roe)'
  )
  
  res = list('Reo'= t(Reo), 'pvalue'=pvalue, 'annote'=annote.heatmap)
  return(res)
}


method_names = c(
  "tenX","tenXarc","sci-ATAC","sciMAP-ATAC", "snATAC","dscATAC","s3ATAC","txciATAC",   
  "sciATAC3","SPATAC","scifi", "HD-ATAC",
  "AS240726", "AS240606", "AS240806", "AS240823", "AS240815", "AS241017", "AS240903", "AS241018", "AS240905", "AS240722",
  "AS241009", "AS240617", "AS240425", "AS240828", "AS240819", "AS240821", "AS241022", "AS240906", "AS240731", "AS240709"
)
method_color = c(
  "#1C6E8C", # Deep blue-green, muted cool tone
  "#4A7B33", # Deep olive green, low saturation
  "#A06723", # Deep brown-orange, muted warm tone
  "#5C3B8C", # Deep purple, reduced brightness
  "#357D74", # Dark teal, soft cool tone
  "#8C7B00", # Deep golden, darker warm tone
  "#5A2A2A", # Deep reddish brown, avoids conflict with the first color
  "#3D7FF0", # Dark cyan-green, lower saturation
  "#274B9F", # Deep blue, low brightness
  "#7A3D5C", # Deep pink, soft yet distinctive
  "#505050",  # Dark gray, used as a neutral color
  "#FF0000", # Bright red, used as the most prominent color
  "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
  "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000"
)
names(method_color) = method_names

celltype_round1_color = c('Immune'='#D7191C',
                 'Epithelium'='#1F78B4',
                 'Endothelium'='#FEE08B',
                 'Stroma'='#A6761D')
method_color_keep = c(
  'HD-ATAC'="#FF0000",  # Deep blue - classic and steady (suitable for control group)
  'SPATAC'="#59A14F",  # Natural green - soft but clear
  'CH-ATAC-seq'="#EDC948",  # Golden yellow - bright but not glaring
  'txciATAC'="#B07AA1",  # Violet - neutral and distinctive
  'dsciATACseq'="#4E79A7",  # Bright red - highlight color (key experimental group)
  'tenXarc'="#76B7B2",  # Blue-green - fresh and understated
    'tenX'="#76B7B2",  # Blue-green - fresh and understated
  'sciMAP-ATAC'="#FF9DA7",  # Light pink - soft contrast
  'snATAC'="#9C755F"   # Tan - neutral background color
)

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
            legend.key.width = unit(3, "mm"),  # width
            legend.key.height = unit(3, "mm"),    # height
            legend.key.spacing = unit(2, 'pt'),
            
            legend.text = element_text(size=6),
            legend.title = element_text(size=6,face = "plain"),
            legend.margin = margin(0, 0, 0, 0)
            )

output_dir='./placeholder_output/figure1/'

setwd('./placeholder_data/download_data/')

sample_info = read.csv('./placeholder_data/sample_information/sample_information.csv',row.names = 1)
# Replicates of in-house data
umap = read.csv(glue('./placeholder_analysis/round_cluster02/round0/cell_meta.csv'))

best_leiden = 'leiden_0.1_round1'
umap$leiden = umap[,best_leiden]

umap$leiden = factor(umap$leiden, levels = c(0:(length(unique(umap$leiden))-1)))

umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
rownames(umap) = umap$X

pseudoObj_peaks <- readRDS(glue('./placeholder_analysis/round_cluster02/round0/pseudoObj_peaks.rds'))
pseudoObj_cells <- readRDS(glue('./placeholder_analysis/round_cluster02/round0/pseudoObj_cells.rds'))
pseudoObj_cells@meta.data[umap$X, colnames(umap)] = umap

leiden_levels = 0:(length(unique(pseudoObj_cells@meta.data[,best_leiden]))-1)

pseudoObj_cells$leiden = factor(pseudoObj_cells@meta.data[,best_leiden], levels=leiden_levels)

Idents(pseudoObj_cells) = pseudoObj_cells$leiden

pseudoObj_cells@meta.data[,'nCount_ATAC'] = pseudoObj_cells$n_fragment
pseudoObj_cells@meta.data[,'nFeature_ATAC'] = pseudoObj_cells$n_peaks

a=FragmentHistogram(object = pseudoObj_cells, group.by = 'Time') & 
    scale_fill_manual(values = time_color2)
a

ggsave(glue('{output_dir}/method_comparison_FragmentSize1.pdf'),a,
       width=180, height=140, units='mm', dpi=600, bg='transparent')

pseudoObj_cells_tsse = TSSEnrichment(object = pseudoObj_cells,fast=F, n=5000)

a=TSSPlot(pseudoObj_cells_tsse, group.by = 'Time') & 
    scale_color_manual(values = time_color2)
a

ggsave(glue('{output_dir}/method_comparison_TSSe2.pdf'),a,
       width=180, height=140, units='mm', dpi=600, bg='transparent')

pseudoObj_cells[['group']] = 'HD-ATAC Aggregrate'
pseudoObj_cells@meta.data[pseudoObj_cells$sampleID=='AS240709','group'] = 'HD-ATAC P0 rep1'
pseudoObj_cells@meta.data[pseudoObj_cells$sampleID=='AS241009','group'] = 'HD-ATAC P0 rep2'
Idents(pseudoObj_cells) = pseudoObj_cells$group

tmp_srt = subset(pseudoObj_cells, group%in%c('HD-ATAC P0 rep1', 'HD-ATAC P0 rep2'))
fgs = Fragments(tmp_srt)
new_fgs = list()
for(i in 1:length(fgs)){
    if(length(fgs[[i]]@cells)>0){
        new_fgs = c(new_fgs,fgs[[i]])
    }
}
new_fgs

peak_matrix = FeatureMatrix(new_fgs, granges(pseudoObj_peaks))

saveRDS(peak_matrix, '../pycode/round_cluster02/merge/replication_matrix.rds')

peak_matrix = readRDS('../pycode/round_cluster02/merge/replication_matrix.rds')

set.seed(1234)
sample_peaks = sample(rownames(peak_matrix), 60000)

rep1_count = rowSums(peak_matrix[sample_peaks,
                               rownames(tmp_srt@meta.data[tmp_srt$group=='HD-ATAC P0 rep1',])])

rep2_count = rowSums(peak_matrix[sample_peaks,
                               rownames(tmp_srt@meta.data[tmp_srt$group=='HD-ATAC P0 rep2',])])


rep_count = data.frame('rep1'=log2((rep1_count)/sum(rep1_count)*1e6+1),
                       'rep2'=log2((rep2_count)/sum(rep2_count)*1e6+1)
                      )

head(rep_count)

a = rep_count %>%
    # head(10000) %>%
    filter(rep1<7, rep2<7) %>%
    ggplot(aes(x=rep1, y=rep2))+
    ggpointdensity::geom_pointdensity(aes(color=..ndensity..), size=0.1,method.args=list(n=1000))+
    #geom_abline(slope = 1, intercept = 0)+
    theme_classic()+
    scale_color_distiller(palette = "Spectral", direction = -1, breaks=c(0.0,0.5,1), limits=c(0,1))+
    stat_cor(size=label_size(6))+
    labs(x='Replication1(log2CPM)',y='Replication2(log2CPM)')+
    mytheme+
    coord_equal()+
    theme(legend.position=c(0.8,0.2),
          legend.key.width = unit(1.2, "mm"),  # width
            legend.key.height = unit(1.2, "mm"),    # height
         )

a

ggsave(glue('{output_dir}/method_comparison_reproducibility.pdf'), a,
       width=39, height=39, units='mm', dpi=600, bg='transparent')

# Relationship between barcode and fragments overlapping peaks/unique reads to confirm sample quality 

sample='AS241022'
tmp_sc = read.csv(glue('../aggr_cellranger/{sample}/outs/singlecell.csv'))
tmp_sc = tmp_sc[order(tmp_sc$peak_region_fragments, decreasing = T),]

tmp_sc$index = sapply(tmp_sc$barcode, function(x)strsplit(x,'-')[[1]][2])


tmp_sc2 = tmp_sc[tmp_sc$index=='4',]
thr= min(tmp_sc2[tmp_sc2$is__cell_barcode==1,'peak_region_fragments'])
tmb_barcode2 = tmp_sc2[tmp_sc2$is__cell_barcode==0&tmp_sc2$peak_region_fragments>thr,
                       'barcode']
tmp_sc2 = tmp_sc2[!tmp_sc2$barcode%in%tmb_barcode2,]
tmp_sc2$barcodes_num = 1:nrow(tmp_sc2)
tmp_sc2$is_cell = ifelse(tmp_sc2$is__cell_barcode, 'Cells', 'Non-cells')
a=tmp_sc2 %>%
  ggplot(aes(x=barcodes_num,y=peak_region_fragments))+
  geom_line(aes(color=as.character(is_cell)), size=1)+
  scale_color_manual(values = c('Cells'='#F39800', 'Non-cells'='gray'), name='is_cell')+
  geom_hline(yintercept = thr, linetype='dashed', color='gray')+
  #facet_wrap(~index, ncol=6, scales = 'free')+
  scale_y_log10()+
  scale_x_log10()+
  labs(x='Barcodes(log10)', y='Fragments Overlapping Peaks (log10)')+
  theme_classic()+
  mytheme+
    theme(legend.position=c(0.4,0.2),
      #legend.key.width = unit(1.2, "mm"),  # width
      #  legend.key.height = unit(1.2, "mm"),    # height
     )

options(repr.plot.width=12, repr.plot.height=8)
a 

ggsave(glue('{output_dir}/method_comparison_barcode_quality.pdf'), a,
       width=43, height=39, units='mm', dpi=600, bg='transparent')

keep_method = c('HD-ATAC', 'SPATAC',  'txciATAC','tenXarc', 'dsciATACseq')#,'CH-ATAC-seq')#,'sciMAP-ATAC', 'snATAC')'scifi','sciATAC3',

all_method_qc = readRDS('./placeholder_processed/all.rds')

table(all_method_qc$method)


tsse_plot=all_method_qc %>%
  filter(!is.na(TSSe)) %>%
  filter(method%in%keep_method) %>%
  mutate(method=fct_reorder(method, TSSe)) %>%
  ggplot(aes(x=method, y=TSSe))+
  geom_boxplot(aes(color=method),outlier.shape = NA, size=0.5)+
  scale_color_manual(values = method_color_keep)+
  stat_compare_means(comparisons=list(c('HD-ATAC','txciATAC')), label = "p.format", label.y = 12, size=label_size(6))+
  lims(y=c(0,18))+
  theme_bw()+
   mytheme+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 60, hjust=1,vjust=1),
        legend.position = 'none')+
    NoLegend()
options(repr.plot.width=6, repr.plot.height=4)
tsse_plot
# all_method_qc2 = all_method_qc
# all_method_qc2[all_method_qc2$method=='HD-ATAC', 'method']=all_method_qc2[all_method_qc2$method=='HD-ATAC', 'sample']
# all_method_qc2 %>%
#   filter(!is.na(TSSe)) %>%
#   mutate(method=fct_reorder(method, TSSe)) %>%
#   ggplot(aes(x=method, y=TSSe))+
#   geom_boxplot(aes(color=method),outlier.shape = NA, size=0.5)+
#   scale_color_manual(values = method_color)+
#   lims(y=c(0,30))+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.text.x=element_text(angle = 60, hjust=1,vjust=1),
#         legend.position = 'none')

ggsave(glue('{output_dir}/method_comparison_TSSe_box_plot.pdf'), tsse_plot,
       width=80, height=40, units='mm', dpi=600, bg='transparent')

Fritss_plot=all_method_qc %>%
  filter(!is.na(FRiTSS)) %>%
  filter(method%in%keep_method) %>%
  mutate(method=fct_reorder(method, FRiTSS)) %>%
  ggplot(aes(x=method, y=FRiTSS))+
  geom_boxplot(aes(color=method),outlier.shape = NA, size=0.5)+
  scale_color_manual(values = method_color_keep)+
   # stat_compare_means(comparisons=list(c('HD-ATAC','txciATAC')), label = "p.format",  size=label_size(6))+
  #lims(y=c(0,30))+
  theme_bw()+
    mytheme+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 60, hjust=1,vjust=1),
        legend.position = 'none')
Fritss_plot

ggsave(glue('{output_dir}/method_comparison_FRiTSS_box_plot.pdf'), Fritss_plot,
       width=60, height=40, units='mm', dpi=600, bg='transparent')

head(all_method_qc)

Counts_plot=all_method_qc %>%
  filter(!is.na(Unique)) %>%
  filter(method%in%keep_method) %>%
  mutate(method=fct_reorder(method, Unique)) %>%
  ggplot(aes(x=method, y=log10(Unique)))+
  geom_boxplot(aes(color=method),outlier.shape = NA, size=0.5)+
  scale_color_manual(values = method_color_keep)+
  #lims(y=c(0,30))+
  theme_bw()+
    mytheme+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 60, hjust=1,vjust=1),
        legend.position = 'none')
Counts_plot

ggsave(glue('{output_dir}/method_comparison_UniqueFrag_box_plot.pdf'), Counts_plot,
       width=60, height=40, units='mm', dpi=600, bg='transparent')

all_method_qc = readRDS('./placeholder_processed/all.rds')

Firp_plot = all_method_qc %>%
  filter(method%in%keep_method) %>%
  mutate(method=fct_reorder(method, FRiP)) %>%
  ggplot(aes(x=method, y=FRiP))+
  geom_boxplot(aes(color=method),outlier.shape = NA, size=0.5)+
  scale_color_manual(values = method_color_keep)+
   mytheme+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 60, hjust=1,vjust=1),
        legend.position = 'none')+

    NoLegend()
options(repr.plot.width=6, repr.plot.height=4)
Firp_plot

ggsave(glue('{output_dir}/method_comparison_FRiP_box_plot.pdf'), Firp_plot,
       width=60, height=40, units='mm', dpi=600, bg='transparent')

# -------comparison_different_depthsUnique fragments  --------
all_method_qc = readRDS('./placeholder_processed/all.rds')
all_method_qc[all_method_qc$method=='CH-ATAC-seq',
              'Total'] = all_method_qc[all_method_qc$method=='CH-ATAC-seq', 
                                       'Unique']

tmp_depth = all_method_qc %>% group_by(method) %>% dplyr::summarise(depth=mean(Unique))
all_nPeaks2 = all_method_qc %>% filter(Total<=50000)
all_nPeaks2 = all_nPeaks2 %>% left_join(tmp_depth)
all_nPeaks2$cut = cut(all_nPeaks2$Total, c(0, 1000, 5000, 10000, 20000,30000,40000, 50000),
                      labels = c(1000, 5000, 10000, 20000,30000,40000, 50000)/1000)

#all_nPeaks2$cut = cut(all_nPeaks2$Total, c(800, 1200, 4500,5500, 9500, 10500, 20000))
cell_num = all_nPeaks2 %>%
  #filter(method%in%c('txciATAC', 'scifi', 'SPATAC','HD-ATAC','sciATAC3')) %>%
  group_by(cut, method)%>%
  dplyr::summarise(num = n())
# FRiP
a=all_nPeaks2 %>%
  #filter(method%in%c('txciATAC', 'scifi', 'SPATAC','HD-ATAC','sciATAC3')) %>%
  group_by(cut, method)%>%
  dplyr::summarise(Unique = mean(FRiP)) %>%
  ggplot(aes(x=cut, y=Unique, color=method, group=method))+
  geom_point()+
  geom_line()+
  labs(x='Sequencing depth\n(x1000 reads per cell)', y='FRiP')+
  scale_color_manual(values = method_color)+
  theme_bw()+
  theme(panel.grid = element_blank())
# TSSe
b=all_nPeaks2 %>%
  #filter(method%in%c('txciATAC', 'scifi', 'SPATAC','HD-ATAC','sciATAC3')) %>%
  group_by(cut, method)%>%
  dplyr::summarise(Unique = mean(TSSe)) %>%
  ggplot(aes(x=cut, y=Unique, color=method, group=method))+
  geom_point()+
  geom_line()+
  labs(x='Sequencing depth\n(x1000 reads per cell)', y='TSSe')+
  scale_color_manual(values = method_color)+
  theme_bw()+
  theme(panel.grid = element_blank())

# FRiTSS
c=all_nPeaks2 %>%
  #filter(method%in%c('txciATAC', 'scifi', 'SPATAC','HD-ATAC','sciATAC3')) %>%
  group_by(cut, method)%>%
  dplyr::summarise(Unique = mean(FRiTSS)) %>%
  ggplot(aes(x=cut, y=Unique, color=method, group=method))+
  geom_point()+
  geom_line()+
  labs(x='Sequencing depth\n(x1000 reads per cell)', y='FRiTSS')+
  scale_color_manual(values = method_color)+
  theme_bw()+
  theme(panel.grid = element_blank())

options(repr.plot.width=18, repr.plot.height=4)
a+b+c

ggsave(glue('{output_dir}/method_comparison_different_depths_line_plot.pdf'), a+b+c,
       width=330, height=80, units='mm', dpi=600, bg='transparent')

head(all_nPeaks2)

nPeaks_plot = all_method_qc %>%
  filter(method%in%keep_method) %>%
  mutate(method=fct_reorder(method, nPeaks)) %>%
  ggplot(aes(x=method, y=log10(nPeaks)))+
  geom_boxplot(aes(color=method),outlier.shape = NA, size=0.5)+
  scale_color_manual(values = method_color_keep)+
  lims(y=c(1.5,4.5))+
   mytheme+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 60, hjust=1,vjust=1),
        legend.position = 'none')+

    NoLegend()
options(repr.plot.width=6, repr.plot.height=4)
nPeaks_plot

ggsave(glue('{output_dir}/method_comparison_nPeaks_box_plot.pdf'), nPeaks_plot,
       width=60, height=40, units='mm', dpi=600, bg='transparent')

my_cellnumber = read.table('../pycode/output/qc_recorde.txt', header = T)
sample_info = read.csv('./placeholder_data/sample_information/sample_information.csv',row.names = 1)
my_cellnumber$Time = sample_info[my_cellnumber$sample, 'Time']


# my_cellnumber = my_cellnumber %>% filter(total>15000)
cell_number = my_cellnumber %>% group_by(Time) %>% dplyr::summarise(cells=sum(after_doublet)) %>% as.data.frame()
cell_number$method='HD-ATAC'
cell_number = rbind(cell_number, data.frame(Time='tenX',cells=c(7729,8067,7749,5297,9374,5537,6639), 
                                            method='tenX'))
tmp_meta = readRDS('./placeholder_processed/scifi_seq2.RDS')
#tmp_meta = as.data.frame(tmp_meta[tmp_meta$tech=='scifi',])
cell_number = rbind(cell_number, data.frame(Time='scifi',cells=as.vector(table(tmp_meta$library)), 
                                            method='scifi'))

#tmp_srt = fread('./GSE215901/GSM6645255_SL1_Cell_Ranger_metadata.csv.gz')
#tmp_srt = tmp_srt[(tmp_srt$is_GRCh38_cell_barcode+tmp_srt$is_mm10_cell_barcode)==1,]
#tmp_srt = fread('./GSE215901/GSE215901_RAW/GSM6645256_SL2_Cell_Ranger_metadata.csv.gz')
#tmp_srt = tmp_srt[(tmp_srt$is_GRCh38_cell_barcode+tmp_srt$is_mm10_cell_barcode)==1,]
meta1 = fread('./GSE231708/GSE231708_phased_txci.cc16ko_wt_mm10_lung_metadata.txt.gz')
meta2 = fread('./GSE231708/GSE231708_standard_txci.cc16ko_wt_mm10_lung_metadata.txt.gz')
meta3 = fread('./GSE231708/GSE231708_standard_txci.native_hg38_lung_metadata.txt.gz')
meta4 = fread('./GSE231708/GSE231708_standard_txci.native_mm10_liver_metadata.txt.gz')
meta5 = fread('./GSE231708/GSE231708_standard_txci.native_mm10_lung_metadata.txt.gz')

cell_number = rbind(cell_number, data.frame(Time='txciATAC',cells=c(as.vector(table(meta1$Replicate)),
                                                                 as.vector(table(meta2$Replicate)),
                                                                 as.vector(table(meta3$Replicate)),
                                                                 as.vector(table(meta4$Replicate)),
                                                                 as.vector(table(meta5$Replicate))
), 
                                            method='txciATAC'))
tmp_meta = readRDS('./placeholder_processed/sciATAC_seq3.RDS')
cell_number = rbind(cell_number, data.frame(Time='sciATAC3',cells=as.vector(table(tmp_meta$donor_id)), 
                                            method='sciATAC3'))


cell_number = rbind(cell_number, data.frame(Time='SPATAC',cells=c(9696,10203,9244,22487,23257,
                                                                  15751,17723,33078,16138,47518,
                                                                  73396,33079,72331,111471,60271,
                                                                  59929,43232,40175,101793,87645), 
                                            method='SPATAC'))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(cell_number, varname="cells", 
                    groupnames=c("method"))


# Convert dose to a factor variable
a = cell_number %>%
  group_by(method) %>%
  dplyr::summarise(cells=round(mean(cells))) %>%
  mutate(method = factor(method, levels=c('HD-ATAC','SPATAC','txciATAC','scifi','sciATAC3',  'tenX'))) %>%
  ggplot(aes(x=method, y=cells))+
  geom_bar(aes(fill=method), color='black',stat='identity')+
  geom_text(aes(y=cells-3000,label=cells), color='white')+
  geom_errorbar(aes(x=method,ymin=cells,ymax=cells+sd),data=df2, width=0.2, inherit.aes = F)+
  scale_fill_manual(values = method_color)+
  theme_bw()+
  theme(panel.grid = element_blank())

# b = cell_number %>%
#   mutate(Time=fct_reorder(Time,cells)) %>%
#   ggplot(aes(x=Time, y=cells))+
#   geom_bar(aes(fill=method), color='black',stat='identity', size=0.3)+
#   geom_hline(yintercept = 8000, linetype='dashed', color='black')+
#   geom_text(aes(x='tenX',y=8000+3000,label='8k'), inherit.aes = F, size=4)+
#   geom_hline(yintercept = 80000, linetype='dashed', color='black')+
#   geom_text(aes(x='tenX',y=80000+3000,label='80k'), inherit.aes = F, size=4)+
#   scale_fill_manual(values = method_color)+
#   
#   theme_bw()+
#   theme(panel.grid = element_blank())

options(repr.plot.width=6, repr.plot.height=4)
a

ggsave(glue('{output_dir}/method_comparison_throughput.pdf'), a,
       width=140, height=80, units='mm', dpi=600, bg='transparent')

df3 = df2 %>% filter(method%in%c('HD-ATAC','tenX'))

df3[2,1] = 'tenXarc'

cell_number2 = cell_number
cell_number2[cell_number2$method=='tenX', 'method'] = 'tenXarc'

# Convert dose to a factor variable
cell_number_plot = cell_number2 %>%
  group_by(method) %>%
  dplyr::summarise(cells=round(mean(cells))) %>%
  filter(method%in%c('HD-ATAC','tenXarc')) %>%
  mutate(method = factor(method, levels=c('HD-ATAC', 'tenXarc'))) %>%
  ggplot(aes(x=method, y=cells/10000))+
  geom_bar(aes(fill=method), color='black',stat='identity', size=0.1)+
  geom_text(aes(y=(cells-3000)/10000,label=cells), color='white', size=label_size(6))+
  geom_errorbar(aes(x=method,ymin=cells/10000,ymax=(cells+sd)/10000),data=df3, width=0.2, inherit.aes = F)+
  scale_fill_manual(values = method_color_keep)+
     mytheme+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 60, hjust=1,vjust=1),
        legend.position = 'none')+

    NoLegend()
cell_number_plot

ggsave(glue('{output_dir}/method_comparison_throughput2.pdf'), a,
       width=24, height=40, units='mm', dpi=600, bg='transparent')

options(repr.plot.width=8, repr.plot.height=4)
all_qc_plot = (cell_number_plot|tsse_plot|Firp_plot|nPeaks_plot|Fritss_plot|Counts_plot)+
    patchwork::plot_layout(widths = c(0.6,1,1,1,1,1))
all_qc_plot

ggsave(glue('{output_dir}/method_comparison_all.pdf'), all_qc_plot,
       width=160, height=35, units='mm', dpi=600, bg='transparent')

base_dir = './placeholder_collision/AS240118A1_Combine2/'

peaks.path  = glue('{base_dir}/filtered_peak_bc_matrix/peaks.bed')
peaks <- read.table(file = peaks.path,col.names = c("chr", "start", "end"))
gr <- makeGRangesFromDataFrame(peaks)
fragPath = glue('{base_dir}/fragments.tsv.gz')
singlecells = glue('{base_dir}/filtered_peak_bc_matrix/barcodes.tsv')
singlecells = read.csv(singlecells, header = F)
singlecells = singlecells[,1]
# Fragments
frags <- CreateFragmentObject(path = fragPath,cells = singlecells)
# Counts
obj.counts <- FeatureMatrix(fragments = frags,
                            features = gr,
                            cells = singlecells,process_n = 1e6)
# Assay
Frag.assay <- CreateChromatinAssay(obj.counts, fragments = frags)
# Object
merge.obj <- CreateSeuratObject(Frag.assay, assay = "ATAC")

merge.obj$group = ifelse(merge.obj$nCount_ATAC>median(merge.obj$nCount_ATAC), 'HighDepth', 'LowDepth')

table(merge.obj$group )

granges(merge.obj)

a=FragmentHistogram(object = merge.obj)
a

ggsave(glue('{output_dir}/method_comparison_FragmentSize_collision.pdf'),a,
       width=80, height=80, units='mm', dpi=600, bg='transparent')

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('GRCh38_chr', seqlevels(annotation))

Annotation(merge.obj) = annotation

merge.obj = TSSEnrichment(object =  merge.obj,fast=F, n=20000)
a=TSSPlot(merge.obj)
a

ggsave(glue('{output_dir}/method_comparison_TSSe_collision.pdf'),a,
       width=100, height=100, units='mm', dpi=600, bg='transparent')

base_dir = './placeholder_collision/AS240118A1_Combine2/'
combine_meta = read.table(glue('{base_dir}/singlecell.csv'), sep=',', header = T)
combine_meta_2 = combine_meta[(combine_meta[,'is_mm10_cell_barcode']+
                               combine_meta[,'is_GRCh38_cell_barcode'])>0,]
new_type = apply(combine_meta_2, 1, function(x){
    thr = 500
    hg_count = as.numeric(x['passed_filters_GRCh38'])
    ms_count = as.numeric(x['passed_filters_mm10'])
    if(min(ms_count,hg_count)>thr){
        return('Collisions')
    }else{
        if(hg_count>ms_count){
            return('GRCh38')
        }else{
            return('mm10')
        }
    }
})

combine_meta_2[, 'type'] = new_type


A549_meta = read.table('./placeholder_collision/AS240118A1_A549_hg38_mm10_singlecell.csv',
                       sep=',', header = T)
AML12_meta = read.table('./placeholder_collision/AS240118A1_AML12_hg38_mm10_singlecell.csv',
                       sep=',', header = T)

A549_meta = A549_meta[(A549_meta[,'is_mm10_cell_barcode']+
                               A549_meta[,'is_GRCh38_cell_barcode'])>0,]
AML12_meta = AML12_meta[(AML12_meta[,'is_mm10_cell_barcode']+
                               AML12_meta[,'is_GRCh38_cell_barcode'])>0,]

A549_meta$source='A549'
AML12_meta$source='AML12'
split_meta = rbind(A549_meta,AML12_meta)

split_meta[, 'type'] = 'GRCh38'
split_meta[split_meta[,'is_mm10_cell_barcode']>0, 'type'] = 'mm10'
split_meta[(split_meta[,'is_mm10_cell_barcode']+
                                    split_meta[,'is_GRCh38_cell_barcode'])==2, 'type'] = 'Collisions'

new_type = apply(split_meta, 1, function(x){
    thr = 500
    hg_count = as.numeric(x['passed_filters_GRCh38'])
    ms_count = as.numeric(x['passed_filters_mm10'])
    if(min(ms_count,hg_count)>thr){
        return('Collisions')
    }else{
        if(hg_count>ms_count){
            return('GRCh38')
        }else{
            return('mm10')
        }
    }
})

split_meta[, 'type'] = new_type

options(repr.plot.width=8, repr.plot.height=8)
options(repr.plot.width=16, repr.plot.height=8)
c=combine_meta_2 %>%
    ggplot(aes(x=passed_filters_GRCh38/1000, y=passed_filters_mm10/1000, color=type)) +
    geom_point(size=1)+
    lims(x=c(0,40), y=c(0,40))+
    labs(x='GRCh38 genome (reads x1,000)',y='mm10 genome (reads x1,000)')+
    scale_color_npg()+
    theme_bw()+
    theme(panel.grid=element_blank())+
    coord_equal()
a=combine_meta_2 %>%
    ggplot(aes(x=log10(passed_filters_GRCh38), y=log10(passed_filters_mm10), color=type)) +
    geom_point(size=1)+
    #lims(x=c(0,40), y=c(0,40))+
    labs(x='GRCh38 genome (log10)',y='mm10 genome (log10)')+
    scale_color_npg()+
    theme_bw()+
    theme(panel.grid=element_blank())+
    coord_equal()
a+c

options(repr.plot.width=16, repr.plot.height=8)
d=split_meta %>%
    ggplot(aes(x=(passed_filters_GRCh38/1000), y=(passed_filters_mm10/1000), color=type)) +
    geom_point(size=1)+
    lims(x=c(0,40), y=c(0,40))+
        labs(x='GRCh38 genome (reads x1,000)',y='mm10 genome (reads x1,000)')+
    scale_color_npg()+
    theme_bw()+
    theme(panel.grid=element_blank())+
    coord_equal()
b=split_meta %>%
    ggplot(aes(x=log10(passed_filters_GRCh38), y=log10(passed_filters_mm10), color=type)) +
    geom_point(size=1)+
    #lims(x=c(0,40), y=c(0,40))+
    labs(x='GRCh38 genome (log10)',y='mm10 genome (log10)')+
    scale_color_npg()+
    theme_bw()+
    theme(panel.grid=element_blank())+
    coord_equal()
b+d

options(repr.plot.width=16, repr.plot.height=16)
ggsave(glue('{output_dir}/method_comparison_collisions_scatter_plot.pdf'), a+b+c+d,
       width=200, height=200, units='mm', dpi=600, bg='transparent')

c=combine_meta_2 %>%
    mutate(type = factor(type, c('mm10', 'GRCh38', 'Collisions'))) %>%
    arrange(type) %>%
    ggplot(aes(x=passed_filters_GRCh38/1000, y=passed_filters_mm10/1000, color=type)) +
    geom_point(size=0.1)+
    lims(x=c(0,40), y=c(0,40))+
    labs(x='GRCh38 genome (reads x1,000)',y='mm10 genome (reads x1,000)')+
    scale_color_manual(values = c('mm10'='gray', 'GRCh38'='#9e9e9e', 'Collisions'='red'))+
    mytheme+
    coord_equal()
d=split_meta %>%
    mutate(type = factor(type, c('mm10', 'GRCh38', 'Collisions'))) %>%
    arrange(type) %>%
    ggplot(aes(x=(passed_filters_GRCh38/1000), y=(passed_filters_mm10/1000), color=type)) +
    geom_point(size=0.1)+
    lims(x=c(0,40), y=c(0,40))+
        labs(x='GRCh38 genome (reads x1,000)',y='mm10 genome (reads x1,000)')+
    scale_color_manual(values = c('mm10'='gray', 'GRCh38'='#9e9e9e', 'Collisions'='red'))+
    mytheme+
    coord_equal()

ggsave(glue('{output_dir}/method_comparison_collisions_scatter_plot2.pdf'), (c+d)+plot_layout(guides = 'collect'),
       width=100, height=50, units='mm', dpi=600, bg='transparent')

plot.df = as.data.frame.table(table(combine_meta_2$type))
colnames(plot.df) = c('Source', 'Freq')
plot.df$percent = plot.df$Freq/sum(plot.df$Freq)
plot.df <- plot.df %>% 
    arrange(desc(Source)) %>%
    mutate(ypos = cumsum(percent)- 0.5*percent)

pl1 = ggplot(plot.df, aes(x="", y=percent, fill=Source))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0) +
  #geom_text(aes(x=.5,y = ypos, label = paste0(round(percent,1),'%')), size=5)+
  geom_text_repel(aes(x=1,y = ypos, label = paste0(Source,'(',round(percent,4)*100,')%')), 
                  size=5, color='black')+
  theme_void()+
  scale_fill_npg()+
  scale_color_igv()+
  #labs(title = 'ChipSeeker') +
  theme( plot.title = element_text(hjust = 0.5),
       legend.position='none')
pl2 = ggplot(plot.df, aes(x="", y=percent, fill=Source))+
 geom_bar(width = 1, stat = "identity")+
 coord_polar("y", start=0) +
 #geom_text(aes(x=1.5,y = ypos, label = paste0(round(percent,1),'%')), size=5)+
 geom_text_repel(aes(x=1,y = ypos, label = paste0(Source,'(',round(percent,4)*100,')%')), 
                 size=5, color='black')+
 theme_void()+
 scale_fill_manual(values=c('mm10'='gray', 'GRCh38'='#9e9e9e', 'Collisions'='red'))+
 #scale_color_npg()+
 #labs(title = 'ChipSeeker') +
 theme( plot.title = element_text(hjust = 0.5),
      legend.position='none')



plot.df = as.data.frame.table(table(split_meta$type))
colnames(plot.df) = c('Source', 'Freq')
plot.df$percent = plot.df$Freq/sum(plot.df$Freq)
plot.df <- plot.df %>% 
    arrange(desc(Source)) %>%
    mutate(ypos = cumsum(percent)- 0.5*percent)

pl3 = ggplot(plot.df, aes(x="", y=percent, fill=Source))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0) +
  #geom_text(aes(x=.5,y = ypos, label = paste0(round(percent,1),'%')), size=5)+
  geom_text_repel(aes(x=1,y = ypos, label = paste0(Source,'(',round(percent,4)*100,')%')), 
                  size=5, color='white')+
  theme_void()+
  scale_fill_npg()+
  scale_color_igv()+
  #labs(title = 'ChipSeeker') +
  theme( plot.title = element_text(hjust = 0.5),
       legend.position='none')
pl4 = ggplot(plot.df, aes(x="", y=percent, fill=Source))+
 geom_bar(width = 1, stat = "identity")+
 coord_polar("y", start=0) +
 #geom_text(aes(x=1.5,y = ypos, label = paste0(round(percent,1),'%')), size=5)+
 geom_text_repel(aes(x=1,y = ypos, label = paste0(Source,'(',round(percent,4)*100,')%')), 
                 size=5, color='white')+
 theme_void()+
 scale_fill_manual(values=c('mm10'='gray', 'GRCh38'='#9e9e9e', 'Collisions'='red'))+
 #scale_color_npg()+
 #labs(title = 'ChipSeeker') +
 theme( plot.title = element_text(hjust = 0.5),
      legend.position='none')

options(repr.plot.width=6, repr.plot.height=6)
pl2+pl4

ggsave(glue('{output_dir}/method_comparison_collisions_pie_chart.pdf'), pl2+pl4,
       width=150, height=75, units='mm', dpi=600, bg='transparent')

sample_info = read.csv('./placeholder_data/sample_information/sample_information.csv',row.names = 1)

base_path='./placeholder_analysis/round_cluster02/'
merge_path='./placeholder_analysis/round_cluster02/merge'

umap = read.csv(glue('{merge_path}/cell_meta.csv'))
tsne = umap

umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)

tsne[, c('UMAP_1', 'UMAP_2')] = tsne[, c('TSNE_1', 'TSNE_2')]

dim(umap)

table(umap$Celltype_round1)

xx=round(table(umap$Celltype_round1)/dim(umap)[1],4)
xx

xx = table(umap$Celltype_round4)
xx

min(xx)
max(xx)
median(xx)
mean(xx)
length(xx)

options(repr.plot.width=9, repr.plot.height=8)
a=myDimPlot(umap, groupby='Celltype_round1', label=F, group_color=celltype_round1_color, point_size=0.1)
a

ggsave(glue('{output_dir}/major_cluster_umap_round1_legend.pdf'), get_legend(a),
       width=50, height=50, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/major_cluster_umap_round1_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

ggsave(glue('{output_dir}/major_cluster_umap_round1_cells.png'), a,
       width=180, height=160, units='mm', dpi=600, bg='transparent')

# round1 TSNE
a=myDimPlot(tsne, groupby='Celltype_round1', label=F, group_color=celltype_round1_color, point_size=0.1)

ggsave(glue('{output_dir}/major_cluster_tsne_round1_legend.pdf'), get_legend(a),
       width=50, height=50, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/major_cluster_tsne_round1_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

# round4 TSNE
a=myDimPlot(tsne, groupby='Celltype_round4', label=F, group_color=expanded_colors, point_size=0.1)
ggsave(glue('{output_dir}/major_cluster_tsne_round4_legend.pdf'), get_legend(a),
       width=450, height=250, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/major_cluster_tsne_round4_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

# Round4 UMAP
a=myDimPlot(umap, groupby='Celltype_round4', label=T, group_color=expanded_colors, point_size=0.1)
ggsave(glue('{output_dir}/major_cluster_umap_round4_legend.pdf'), get_legend(a),
       width=450, height=250, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/major_cluster_umap_round4_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

# UMAP Time
a = myDimPlot(umap, groupby='Time', label=FALSE, group_color=time_color2, point_size=0.1)
ggsave(glue('{output_dir}/major_cluster_umap_Time_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/major_cluster_umap_Time_cells.pdf'), a+NoLegend(),
       width=60, height=60, units='mm', dpi=600, bg='transparent')

xx = table(umap$Celltype_round4)
cluster_num = 1:129
names(cluster_num) = names(xx[order(xx)])

cluster_num[cluster_num==22]

umap$cluster_num = cluster_num[umap$Celltype_round4]

celltype_meta = umap %>% 
    group_by(Celltype_round1,Celltype_round2,Celltype_round3,Celltype_round4) %>%
    dplyr::summarise(n=n())
celltype_meta$cluster_num = cluster_num[celltype_meta$Celltype_round4]
write.table(celltype_meta, './placeholder_analysis/round_cluster02/merge/celltype_modify.csv',
            sep=',', row.names=F)

set.seed(1234)
sample_ix = sample(1:nrow(umap), 200000)

head(umap)

write.table(umap[sample_ix, 'X'], './placeholder_analysis/round_cluster02/merge/sample20wcells.txt',
           quote=F)

round4_color = expanded_colors[1:129]
names(round4_color) = names(xx)

saveRDS(cluster_num,'./placeholder_analysis/round_cluster02/merge/round4_cluster_num.rds')

saveRDS(round4_color,'./placeholder_analysis/round_cluster02/merge/round4_color.rds')



a=cluster_num[cluster_num==108]
round4_color[names(a)] = '#FFa90AFF'
a=cluster_num[cluster_num==77]
round4_color[names(a)] = '#136d33FF'
a=cluster_num[cluster_num==122]
round4_color[names(a)] = '#296223FF'

round4_color2 = round4_color
names(round4_color2) = 1:129

a=cluster_num[cluster_num%in%c(122,60)]
a
round4_color[names(a)] 

library(ggtext)
darken_color <- function(color, factor = 1.5) {
  col <- col2rgb(color)
  col <- col / factor  # 
  col <- rgb(t(col), maxColorValue = 255)
  return(col)
}

label_data = umap[sample_ix, ] %>% 
    group_by(cluster_num) %>% 
    dplyr::summarise(x=median(TSNE_1),y=median(TSNE_2), Celltype_round4=unique(Celltype_round4))

a=umap[sample_ix, ] %>%
  ggplot(aes(x=TSNE_1, y=TSNE_2))+
  geom_point(aes(,fill=Celltype_round4),size=0.3, shape=21, stroke=NA)+
  #geom_point(aes(x=x,y=y,label=cluster_num,color=Celltype_round4),
  #           data=label_data, shape=21, stroke=0.5, size=5,
  #           fill=alpha('#eeeeee', 0.4),)+
  #geom_text(aes(x=x,y=y,label=cluster_num), color='black',          
  #           data=label_data,size=label_size(4))+
  geom_label_repel(aes(x=x,y=y,label=cluster_num,color=Celltype_round4),data=label_data,
                   fill=alpha('#eeeeee', 0.5),
                   #color='black',
                   #label.color='black',segment.colour='black',
                   size=label_size(6), #label.r=unit(0.25, 'lines'),
                   #label.size=label_size(0.5),
                   label.padding=0.1,
                  segment.size=0.3,
                  max.overlaps=1000,
                   force=0.008,force_pull=0.008
                  )+
  scale_fill_manual(values = round4_color)+
  scale_color_manual(values = darken_color(round4_color, factor=1.15))+
  theme_void()+
  guides(fill = guide_legend(override.aes = list(size=4, stroke=NA,shape=21)))

ggsave(glue('{output_dir}/major_cluster_tsne_round4_cells_raster.pdf'), a+NoLegend(),
       width=100, height=100, units='mm', dpi=1200, bg='transparent')

a=umap[sample_ix, ] %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2,fill=Time))+
  geom_point(size=0.3, shape=21, stroke=NA)+
  scale_fill_manual(values = time_color2)+
  theme_void()+
  guides(fill = guide_legend(override.aes = list(size=4, stroke=NA,shape=21)))

ggsave(glue('{output_dir}/major_cluster_umap_time_cells_raster.pdf'), a+NoLegend(),
       width=60, height=60, units='mm', dpi=1200, bg='transparent')

a=umap[sample_ix, ] %>%
  ggplot(aes(x=TSNE_1, y=TSNE_2,fill=Time))+
  geom_point(size=0.2, shape=21, stroke=NA)+
  scale_fill_manual(values = time_color2)+
  theme_void()+
  guides(fill = guide_legend(override.aes = list(size=4, stroke=NA,shape=21)))

ggsave(glue('{output_dir}/major_cluster_tsne_time_cells_raster.pdf'), a+NoLegend(),
       width=60, height=60, units='mm', dpi=1200, bg='transparent')

# TSNE Time
a = myDimPlot(tsne, groupby='Time', label=FALSE, group_color=time_color2, point_size=0.1)
ggsave(glue('{output_dir}/major_cluster_tsne_Time_legend.pdf'), get_legend(a),
       width=150, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/major_cluster_tsne_Time_cells.png'), a+NoLegend(),
       width=60, height=60, units='mm', dpi=600, bg='transparent')

# -----sample_qc_bar_plot----------------------------------------------------------
qc_res = read.table('../pycode/output/qc_recorde.txt', header = T)
qc_res = reshape2::melt(qc_res, id='sample')
sample_info = read.csv('./placeholder_data/sample_information/sample_information.csv',row.names = 1)
qc_res$sample_id = factor(sample_info[qc_res$sample, 'sample_time'], levels=sample_levels)
qc_res$Time = factor(sample_info[qc_res$sample, 'Time'], levels=time_levels)

qc_res = qc_res%>%dplyr::group_by(Time, variable) %>% dplyr::summarise(value=sum(value))

qc_res %>%
  filter(variable=='after_doublet')

xx = qc_res %>%
  filter(variable=='after_doublet')
sum(xx$value)

xx=umap %>%
    group_by(Time) %>%
    dplyr::summarise(n=n())
min(xx$n)
max(xx$n)
mean(xx$n)

xx=umap %>%
    group_by(Time) %>%
    dplyr::summarise(n=mean(n_fragment))
min(xx$n)
max(xx$n)
mean(xx$n)

options(repr.plot.width=10, repr.plot.height=8)
a=qc_res %>%
  filter(variable=='after_doublet') %>%
  #mutate(variable=factor(variable, levels=rev(c('total', 'after_qc', 'after_doublet')))) %>%
  mutate(Time=factor(Time, rev(time_levels))) %>%
  ggplot(aes(x=value/10000, y=Time,fill=Time))+
  geom_bar(color='black',stat='identity', position = 'dodge', width=0.9, size=pt_to_mm(0.1))+
  geom_text(aes(label=value,x=2),
            size=label_size(5),
            color='#666666'
            #position = position_dodge(width =0.9), hjust=1.2
           )+
  scale_fill_manual(values = time_color2)+
  #scale_fill_manual(values = c('total'='#f0f3db', 'after_qc'='#a8ddb5', 'after_doublet'='#43a2ca'))+
  labs(x='cell number(x10,000)')+
  #scale_fill_manual(values = c('total', 'after_qc', 'after_doublet'))+
  scale_x_continuous(expand = expansion(mult = c(0,0.3)))+
  theme_classic()+
  mytheme+
  NoLegend()
a

options(repr.plot.width=10, repr.plot.height=8)
a=umap %>%
    group_by(Time) %>%
    dplyr::summarise(n=n()) %>%
  #mutate(variable=factor(variable, levels=rev(c('total', 'after_qc', 'after_doublet')))) %>%
  mutate(Time=factor(Time, rev(time_levels))) %>%
  ggplot(aes(x=n/10000, y=Time,fill=Time))+
  geom_bar(color='black',stat='identity', position = 'dodge', width=0.9, size=pt_to_mm(0.1))+
  geom_text(aes(label=n,x=2),
            size=label_size(5),
            color='#666666'
            #position = position_dodge(width =0.9), hjust=1.2
           )+
  scale_fill_manual(values = time_color2)+
  #scale_fill_manual(values = c('total'='#f0f3db', 'after_qc'='#a8ddb5', 'after_doublet'='#43a2ca'))+
  labs(x='cell number(x10,000)')+
  #scale_fill_manual(values = c('total', 'after_qc', 'after_doublet'))+
  scale_x_continuous(expand = expansion(mult = c(0,0.3)))+
  theme_classic()+
  mytheme+
  NoLegend()
a

ggsave(glue('{output_dir}/quality_filtered_cell_counts.pdf'), a,
       width=30, height=60, units='mm', dpi=600, bg='transparent')



tmp_qc = umap %>%
    group_by(Time) %>%
    dplyr::summarise(n=n()) %>%
  #mutate(variable=factor(variable, levels=rev(c('total', 'after_qc', 'after_doublet')))) %>%
  mutate(Time=factor(Time, rev(time_levels)))

qc_res[qc_res$variable=='after_doublet', 'value'] = tmp_qc$n

a=qc_res %>%
  mutate(variable=factor(variable, levels=(c('total', 'after_qc', 'after_doublet')))) %>%
  ggplot(aes(y=value, x=Time,fill=variable))+
  geom_bar(color='white',stat='identity', position = 'dodge', width=0.9)+
  geom_text(aes(label=value),
            size=3,
            position = position_dodge(width =0.9), hjust=-0.1, angle=90)+
  scale_fill_bmj()+
  labs(x='cell number')+
  #scale_fill_manual(values = c('total', 'after_qc', 'after_doublet'))+
  #scale_x_continuous(expand = expansion(mult = c(0,0.3)))+
  theme_classic()
a

ggsave(glue('{output_dir}/quality_filtered_cell_steps.pdf'), a,
       width=180, height=60, units='mm', dpi=600, bg='transparent')

# -----sample_qc_violin_plot----------------------------------------------------------
# define lable fun
give.n = function(x) {
  la.df = data.frame(y = mean(x)+ max(x) / 10,label = round(mean(x), 2));
  return(la.df)
}
# umap = read.csv('./placeholder_analysis/round_cluster02/merge/cell_meta.csv')
# sample_info = read.csv('./placeholder_data/sample_information/sample_information.csv',row.names = 1)
# umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
# umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
# rownames(umap) = umap$X

median(umap$n_fragment)

mean(umap$n_fragment)

mean(umap$tsse)

a=umap %>%
  ggplot(aes(x=Time, y=n_fragment))+
  geom_violin(aes(fill=Time), size=0.1,trim = TRUE)+
  scale_fill_manual(values = time_color2)+
  theme_classic()+
  #stat_summary(fun = mean, geom = "point", col = "black") +  # Add points to plot
  #stat_summary(fun.data = give.n,size = 3,position = 'identity',
  #             geom = "text",
  #             col = "black")+
  #theme(axis.text.x = element_text(angle=60, hjust = 1,vjust = 1))+
  NoLegend()
b=umap %>%
  ggplot(aes(x=Time, y=n_peaks))+
  geom_violin(aes(fill=Time), size=0.1,trim = TRUE)+
  scale_fill_manual(values = time_color2)+
  theme_classic()+
  #stat_summary(fun = mean, geom = "point", col = "black") +  # Add points to plot
  #stat_summary(fun.data = give.n,size = 3,position = 'identity',
  #             geom = "text",
  #             col = "black")+
  #theme(axis.text.x = element_text(angle=60, hjust = 1,vjust = 1))+
  NoLegend()
c=umap %>%
  ggplot(aes(x=Time, y=tsse))+
  geom_violin(aes(fill=Time), size=0.1,trim = TRUE)+
  scale_fill_manual(values = time_color2)+
  theme_classic()+
  #stat_summary(fun = mean, geom = "point", col = "black") +  # Add points to plot
  #stat_summary(fun.data = give.n,size = 3,position = 'identity',
  #             geom = "text",
  #             col = "black")+
  #theme(axis.text.x = element_text(angle=60, hjust = 1,vjust = 1))+
  NoLegend()
d=umap %>%
  ggplot(aes(x=Time, y=Frip))+
  geom_violin(aes(fill=Time), size=0.1,trim = TRUE)+
  scale_fill_manual(values = time_color2)+
  theme_classic()+
  #stat_summary(fun = mean, geom = "point", col = "black") +  # Add points to plot
  #stat_summary(fun.data = give.n,size = 3,position = 'identity',
  #             geom = "text",
  #             col = "black")+
  #theme(axis.text.x = element_text(angle=60, hjust = 1,vjust = 1))+
  NoLegend()

options(repr.plot.width=14, repr.plot.height=6)
a = ((a/b)|(c/d)) #+ plot_layout(align='hv')
#a=::plot_grid(a,b,c,d, ncol=2)
a

ggsave(glue('{output_dir}/quality_metric_violin.pdf'), a,
       width=240, height=80, units='mm', dpi=600, bg='transparent')

all_umap = read.csv('./placeholder_analysis/round_cluster02/merge/cell_meta.csv')
sample_info = read.csv('./placeholder_data/sample_information/sample_information.csv',row.names = 1)
all_umap$sample_id = factor(sample_info[all_umap$sample, 'sample_time'], levels=sample_levels)
all_umap$Time = factor(sample_info[all_umap$sample, 'Time'], levels=time_levels)
rownames(all_umap) = all_umap$X

base_path='./placeholder_analysis/round_cluster02/round1/Immune/'
umap = read.csv(glue('{base_path}/cell_meta_ann.csv'))
umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
rownames(umap) = umap$X

same_bc = intersect(umap$X,all_umap$X)

umap[same_bc, 'Celltype_round4'] = all_umap[same_bc, 'Celltype_round4']

umap = umap[same_bc, ]

a = myDimPlot(umap, groupby='Time', label=FALSE, group_color=time_color2, point_size=0.4)
ggsave(glue('{output_dir}/subcluster_umap_immune_time_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/subcluster_umap_immune_time_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

a = myDimPlot(umap, groupby='Celltype_round4', label=FALSE, group_color=expanded_colors, point_size=0.4)
ggsave(glue('{output_dir}/subcluster_umap_Immune_round4_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/subcluster_umap_immune_round4_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

a = myDimPlot(umap, groupby='Celltype_round2', label=FALSE, group_color=pal_igv()(10), point_size=1)
ggsave(glue('{output_dir}/subcluster_umap_immune_round2_cells.png'), a,
       width=180, height=160, units='mm', dpi=600, bg='transparent')

base_path='./placeholder_analysis/round_cluster02/round2/Tcell/'
umap = read.csv(glue('{base_path}/cell_meta_ann.csv'))
umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
rownames(umap) = umap$X

same_bc = intersect(umap$X,all_umap$X)

umap[same_bc, 'Celltype_round4'] = all_umap[same_bc, 'Celltype_round4']

umap = umap[same_bc, ]

a = myDimPlot(umap, groupby='Celltype_round3', label=FALSE, group_color=pal_igv()(10), point_size=1)
ggsave(glue('{output_dir}/subcluster_umap_immune_round3_cells.png'), a,
       width=180, height=160, units='mm', dpi=600, bg='transparent')

base_path='./placeholder_analysis/round_cluster02/round3/CD8/'
umap = read.csv(glue('{base_path}/cell_meta_ann.csv'))
umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
rownames(umap) = umap$X

same_bc = intersect(umap$X,all_umap$X)

umap[same_bc, 'Celltype_round4'] = all_umap[same_bc, 'Celltype_round4']

umap = umap[same_bc, ]

a = myDimPlot(umap, groupby='Celltype_round4', label=FALSE, group_color=pal_igv()(10), point_size=1)
ggsave(glue('{output_dir}/subcluster_umap_immune_round4_cells.png'), a,
       width=180, height=160, units='mm', dpi=600, bg='transparent')

base_path='./placeholder_analysis/round_cluster02/round1/Endothelium/'
umap = read.csv(glue('{base_path}/cell_meta_ann.csv'))
umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
rownames(umap) = umap$X

same_bc = intersect(umap$X,all_umap$X)

umap[same_bc, 'Celltype_round4'] = all_umap[same_bc, 'Celltype_round4']

umap = umap[same_bc, ]

a = myDimPlot(umap, groupby='Time', label=FALSE, group_color=time_color2, point_size=0.4)
ggsave(glue('{output_dir}/subcluster_umap_epithelium_time_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/subcluster_umap_epithelium_time_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

a = myDimPlot(umap, groupby='Celltype_round4', label=FALSE, group_color=expanded_colors, point_size=0.4)
ggsave(glue('{output_dir}/subcluster_umap_Epithelium_round4_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/subcluster_umap_Epithelium_round4_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

base_path='./placeholder_analysis/round_cluster02/round1/Epithelium/'
umap = read.csv(glue('{base_path}/cell_meta_ann.csv'))
umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
rownames(umap) = umap$X

same_bc = intersect(umap$X,all_umap$X)

umap[same_bc, 'Celltype_round4'] = all_umap[same_bc, 'Celltype_round4']

umap = umap[same_bc, ]

a = myDimPlot(umap, groupby='Time', label=FALSE, group_color=time_color2, point_size=0.4)
ggsave(glue('{output_dir}/subcluster_umap_Endothelium_Time_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/subcluster_umap_Endothelium_Time_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

a = myDimPlot(umap, groupby='Celltype_round4', label=FALSE, group_color=expanded_colors, point_size=0.4)
ggsave(glue('{output_dir}/subcluster_umap_Endothelium_round4_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/subcluster_umap_Endothelium_round4_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

base_path='./placeholder_analysis/round_cluster02/round1/Stroma//'
umap = read.csv(glue('{base_path}/cell_meta_ann.csv'))
umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
rownames(umap) = umap$X

same_bc = intersect(umap$X,all_umap$X)

umap[same_bc, 'Celltype_round4'] = all_umap[same_bc, 'Celltype_round4']

umap = umap[same_bc, ]

a = myDimPlot(umap, groupby='Time', label=FALSE, group_color=time_color2, point_size=0.4)
ggsave(glue('{output_dir}/subcluster_umap_Stroma_Time_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/subcluster_umap_Stroma_Time_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

a = myDimPlot(umap, groupby='Celltype_round4', label=FALSE, group_color=expanded_colors, point_size=0.4)
ggsave(glue('{output_dir}/subcluster_umap_Stroma_round4_legend.pdf'), get_legend(a),
       width=450, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/subcluster_umap_Stroma_round4_cells.pdf'), a+NoLegend(),
       width=160, height=160, units='mm', dpi=600, bg='transparent')

new_celltype = read.csv('./placeholder_analysis/round_cluster02/merge/celltype_modify.csv')

new_celltype$modify_name = new_celltype$Celltype_round4_new
# new_celltype$modify_name = gsub('Myeloid_', '',new_celltype$modify_name)
# new_celltype$modify_name = gsub('Stroma_', '',new_celltype$modify_name)
# new_celltype$modify_name = gsub('Fibroblast_', 'Fibro_',new_celltype$modify_name)
# new_celltype$modify_name = gsub('Monocytes_', 'Mono_',new_celltype$modify_name)
# new_celltype$modify_name = gsub('Macorphage_', 'Mac_',new_celltype$modify_name)
# 
# new_celltype$modify_name = gsub('Myofibroblast_', 'Myofibro_',new_celltype$modify_name)
# new_celltype$modify_name = gsub('EC_', '',new_celltype$modify_name)
# new_celltype$modify_name = gsub('Airway_', '',new_celltype$modify_name)
# 

head(new_celltype)

all_umap = read.csv('./placeholder_analysis/round_cluster02/merge/cell_meta.csv')
sample_info = read.csv('./placeholder_data/sample_information/sample_information.csv',row.names = 1)
all_umap$sample_id = factor(sample_info[all_umap$sample, 'sample_time'], levels=sample_levels)
all_umap$Time = factor(sample_info[all_umap$sample, 'Time'], levels=time_levels)
rownames(all_umap) = all_umap$X

# data <- data.frame(
#   level1 = "root",
#   level2 = paste0('r1:',new_celltype$MainType),
#   level3 = paste0('r2:',new_celltype$SubType),
#   level4 = paste0('r3:',new_celltype$cluster_num),
#   level5 = paste0('r4:',new_celltype$MainType,'.',new_celltype$SubType,'.',new_celltype$cluster_num)
# )
# data

data <- data.frame(
  level1 = "root",
  level2 = paste0('r1:',new_celltype$MainType),
  level3 = paste0('r2:',new_celltype$Celltype_round2),
  level4 = paste0('r3:',new_celltype$Celltype_round3),
  level5 = paste0('r4:',new_celltype$modify_name,'.',new_celltype$cluster_num)
)

# data <- data.frame(
#   level1 = "root",
#   level2 = paste0('r1:',all_umap$Celltype_round1),
#   level3 = paste0('r2:',all_umap$Celltype_round2),
#   level4 = paste0('r3:',all_umap$Celltype_round3),
#   level5 = paste0('r4:',all_umap$Celltype_round3,':',all_umap$Celltype_round4)
# )

#data = data[data$level2=='r1:Epithelium',]

# build_hierarchical_relationships:convert_data_to_edge_list
edges <- data %>%
  tidyr::pivot_longer(cols = everything(), names_to = "level", values_to = "node") %>%
  group_by(level) %>%
  mutate(parent = lag(node, default = NA)) %>%
  filter(!is.na(parent)) %>%
  distinct(parent, node) %>%
  ungroup()

edges = edges[edges$level!='level1',]

net <- graph_from_data_frame(edges[,c('parent', 'node')], directed=T)


tree <- as.phylo(net)

branch_col=data[,c('level2','level5')]
branch_col = branch_col[!duplicated(branch_col),]
xx=branch_col$level2
names(xx)=branch_col$level5

celltype_tree <- ggtree(tree, branch.length = 'none')

show_name = c()
levels = c()
for(i in celltype_tree$data$label){
    tmp_x = strsplit(i,':')[[1]]
    show_name = c(show_name, tmp_x[length(tmp_x)])
    levels = c(levels,tmp_x[1])
}

celltype_tree$data$show_name=show_name
celltype_tree$data$levels=levels
celltype_tree$data$branch_group = xx[celltype_tree$data$label]

celltype_tree$data$color_group = sapply(celltype_tree$data$show_name, function(x){
    xx = strsplit(x,'\\.')[[1]]
    as.numeric(xx[length(xx)])
})

celltype_tree$data

round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')
round4_num = readRDS('./placeholder_analysis/round_cluster02/merge/round4_cluster_num.rds')
round4_color2 = round4_color
names(round4_color2) = round4_num[names(round4_color2)]

# options(repr.plot.width=12, repr.plot.height=18)
# # use ggtree visualize_tree
# celltype_tree=celltype_tree+aes(color=branch_group)+
#   geom_tiplab(aes(label = show_name)) +
#   geom_text2(aes(label = show_name, subset = !isTip)) +
#     scale_x_continuous(expand = c(0.5,0.1))+
#     scale_color_manual(values=c('New'='black', 'Classical'='gray', 'r1:Immune'='#D7191C',
#                      'r1:Epithelium'='#1F78B4',
#                      'r1:Endothelium'='#FEE08B',
#                      'r1:Stroma'='#A6761D'))
# #celltype_tree
options(repr.plot.width=12, repr.plot.height=18)
# use ggtree visualize_tree
celltype_tree2=celltype_tree+aes(,color=levels)+
  geom_tiplab(aes(label = show_name,color=as.character(color_group)), size=label_size(6)) +
  #geom_text2(aes(label = show_name, subset = !isTip,color=show_name)) +
    scale_x_continuous(expand = c(0.5,0.1))+
    scale_color_manual(values=c('New'='black', 'Classical'='gray', 'r1:Imm'='#D7191C',
                                 'r1:Epi'='#1F78B4',
                                 'r1:EC'='#FEE08B',
                                 'r1:Str'='#A6761D', 
                                'r1'='red', 'r2'='black', 'r3' = 'blue', 'r4'='green',
                                round4_color2))+
    NoLegend()
celltype_tree2

cluster_num = readRDS('./placeholder_analysis/round_cluster02/merge/round4_cluster_num.rds')

label_odder = celltype_tree$data %>% filter(isTip) %>% arrange(y)
# name_order = sapply(label_odder$label, function(x){
#     xx=strsplit(x,':')[[1]]
#     xx[length(xx)]
# })
name_order = names(cluster_num[label_odder$color_group])
cell_prop_mat = all_umap %>%
    dplyr::group_by(Celltype_round4, Time) %>%
    dplyr::summarise(cell_num=n()) %>%
    reshape2::acast(Celltype_round4~Time, value.var='cell_num',fill = 0)

cell_prop_mat = t(t(cell_prop_mat)/colSums(cell_prop_mat))
cell_prop_mat = cell_prop_mat / rowSums(cell_prop_mat)

c=melt(cell_prop_mat) %>% 
    mutate(Celltype_round4=factor(Var1, levels = name_order)) %>%
    mutate(Time=factor(Var2, levels = rev(time_levels))) %>%
    ggplot(aes_string(x = 'Celltype_round4', y='value')) + 
    geom_bar(aes_string(fill = 'Time'), stat = "identity") + 
    scale_fill_manual(values = time_color2) + 
    theme_bw() + 
    coord_flip()+
    labs(title = "Time") + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(),
          axis.text.y = element_blank(),
         plot.title=element_text(hjust=0.5))
d=all_umap %>%
    dplyr::group_by(Celltype_round4,Celltype_round1) %>%
    dplyr::summarise(number=n()) %>%
    mutate(Celltype_round4=factor(Celltype_round4, name_order)) %>%
    ggplot(aes(x=Celltype_round4, y=log10(number),fill=Celltype_round4))+
    geom_bar(stat='identity', color='white', size=0.1)+
    scale_fill_manual(values = round4_color) + 
    theme_bw() + 
    coord_flip()+
    labs(title = "#Cell(log10)") + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          plot.title=element_text(hjust=0.5),
         legend.position='none')

xx = all_umap %>%
    dplyr::group_by(Celltype_round4,Celltype_round1) %>%
    dplyr::summarise(number=n()) 
median(xx$number)

peak_exp_rds = readRDS('./placeholder_analysis/round_cluster02/merge/cluster_gene_avg_expr.rds')

peak_order = peak_exp_rds$right_ann
peak_order$cluster = factor(peak_order$cluster, name_order)
peak_order = peak_order %>% arrange(cluster)

peak_mat = peak_exp_rds$mat[peak_order$gene,name_order]

color_func = circlize::colorRamp2(c(-2,-1.5,0,4,8),c("#2c7bb6",'#abd9e9',"white","#fdae61","#a50026"))
color_range = c(min(peak_mat),max(peak_mat))

peak_exp_plot=melt(peak_mat)%>%
    ggplot(aes(x=Var1,y=Var2, fill=value))+
    geom_tile()+
    scale_fill_gradientn(colors = color_func(color_range[1]:color_range[2]),
                        name='Peak')+
    theme_minimal()+
    labs(title = "Peak") + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(),
          plot.title=element_text(hjust=0.5),
          axis.text.x = element_blank(),
         axis.text.y = element_blank())

l1_peak_exp_plot = get_legend(peak_exp_plot)

##########
#gene_exp = read.csv('./placeholder_root/mouse_lung_project/pycode/round_cluster02_gene/cluster_gene_avg_expr.csv')

gene_exp_rds = readRDS('./placeholder_root/mouse_lung_project/pycode/round_cluster02_gene/cluster_gene_avg_expr.rds')

gene_order = gene_exp_rds$right_ann
gene_order$cluster = factor(gene_order$cluster, name_order)
gene_order = gene_order %>% arrange(cluster)

gene_mat = gene_exp_rds$mat[gene_order$gene,name_order]
gene_mat[gene_mat>2] = 2
gene_mat[gene_mat<(-2)] = -2
color_func = circlize::colorRamp2(c(-2,-1,0,1,2),c("#313695",'#74add1',"#ffffbf","#fdae61","#a50026"))
color_range = c(-2,2)

color_func = circlize::colorRamp2(c(-2,-1,0,1,2),c("#440154",'#3b528b',"#21918c","#5dc863","#a50026"))
color_range = c(min(gene_mat),max(gene_mat))
#color_func = circlize::colorRamp2(c(-2,-1,0,1,2),c("#2c7bb6",'#abd9e9',"white","#fdae61","#a50026"))
#color_range = c(min(gene_mat),max(gene_mat))
gene_exp_plot=melt(gene_mat)%>%
    ggplot(aes(x=Var1,y=Var2, fill=value))+
    geom_tile()+
    scale_fill_gradientn(colors = color_func(color_range[1]:color_range[2]),
                        name='GA')+
    theme_minimal()+
    labs(title = "Gene activity") + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(),
          plot.title=element_text(hjust=0.5),
          axis.text.x = element_blank(),
         axis.text.y = element_blank())

l1_gene_exp_plot = get_legend(gene_exp_plot)

options(repr.plot.width=16, repr.plot.height=12)
( celltype_tree2 | c | d |peak_exp_plot|gene_exp_plot) + 
    plot_layout(guides='collect', widths=c(1,1,0.5,2,2))

a=( celltype_tree2 | c | d |peak_exp_plot|gene_exp_plot) + 
    plot_layout(guides='collect', widths=c(1,1,0.5,2,2))

ggsave(glue('{output_dir}/subcluster_dendrogram_and_accessibility.pdf'), a,
       width=18, height=8, units='in', dpi=600, bg='transparent')

a=( celltype_tree2 | c | d |gene_exp_plot) + 
    plot_layout(guides='collect', widths=c(1,1,0.5,2))

ggsave(glue('{output_dir}/subcluster_dendrogram_and_accessibility2.pdf'), a,
       width=12, height=8, units='in', dpi=600, bg='transparent')

a=( celltype_tree2 | c |d) + 
    plot_layout(guides='collect', widths=c(1,0.3,0.3))

ggsave(glue('{output_dir}/subcluster_dendrogram_and_accessibility3.pdf'), a,
       width=120, height=235, units='mm', dpi=600, bg='transparent')

a=( celltype_tree2 | peak_exp_plot|gene_exp_plot) + 
    plot_layout(guides='collect', widths=c(2,2,2))

ggsave(glue('{output_dir}/subcluster_dendrogram_and_accessibility5.pdf'), a,
       width=16, height=12, units='in', dpi=600, bg='transparent')

ggsave(glue('{output_dir}/subcluster_dendrogram_and_accessibility5.png'), a,
       width=16, height=12, units='in', dpi=600, bg='transparent')

color_func = circlize::colorRamp2(c(-2,-1,0,1,2),c("#440154",'#3b528b',"#21918c","#5dc863","#a50026"))
color_range = c(min(gene_mat),max(gene_mat))
#color_func = circlize::colorRamp2(c(-2,-1,0,1,2),c("#2c7bb6",'#abd9e9',"white","#fdae61","#a50026"))
#color_range = c(min(gene_mat),max(gene_mat))
gene_exp_plot=melt(gene_mat)%>%
    ggplot(aes(x=Var1,y=Var2, fill=value))+
    geom_tile()+
    #scale_fill_gradientn(colors = color_func(color_range[1]:color_range[2]),
    #                    name='GA')+
    scale_fill_viridis()+
    theme_minimal()+
    labs(title = "Gene activity") + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(),
          plot.title=element_text(hjust=0.5),
          axis.text.x = element_blank(),
         axis.text.y = element_blank())
gene_exp_plot

au = all_umap %>%
    dplyr::group_by(Celltype_round4,Celltype_round1) %>%
    dplyr::summarise(number=n()) %>%
    as.data.frame()
rownames(au) = au$Celltype_round4
head(au)

#
cell_prop_mat = all_umap %>%
    dplyr::group_by(Celltype_round4, Time) %>%
    dplyr::summarise(cell_num=n()) %>%
    reshape2::acast(Celltype_round4~Time, value.var='cell_num',fill = 0)

cell_prop_mat = t(t(cell_prop_mat)/colSums(cell_prop_mat))
cell_prop_mat = cell_prop_mat / rowSums(cell_prop_mat)
rownames(new_celltype) = new_celltype$Celltype_round4
rownames(cell_prop_mat) = paste0(new_celltype[rownames(cell_prop_mat), 'modify_name'],
                                '.',
                                new_celltype[rownames(cell_prop_mat), 'cluster_num'])

rownames(new_celltype) = paste0(new_celltype$modify_name,
                                '.',
                                new_celltype$cluster_num)

options(repr.plot.width=8, repr.plot.height=14)
ht = Heatmap(cell_prop_mat,cluster_columns = F,
       row_split = new_celltype[rownames(cell_prop_mat), 'Celltype_round1'],
       row_dend_width= unit(20, "mm"),
        clustering_distance_rows='canberra'
       )
draw(ht)

tree_top = as.phylo(as.hclust(row_dend(ht,on_slice=T)))
tree_1 = as.phylo(as.hclust(row_dend(ht)[['Immune']]))
tree_2 = as.phylo(as.hclust(row_dend(ht)[['Epithelium']]))
tree_3 = as.phylo(as.hclust(row_dend(ht)[['Endothelium']]))
tree_4 = as.phylo(as.hclust(row_dend(ht)[['Stroma']]))

new_tree = bind.tree(tree_top, tree_1, where=which(tree_top$tip.label=='Immune'))
new_tree = bind.tree(new_tree, tree_2, where=which(new_tree$tip.label=='Epithelium'))
new_tree = bind.tree(new_tree, tree_3, where=which(new_tree$tip.label=='Endothelium'))
new_tree = bind.tree(new_tree, tree_4, where=which(new_tree$tip.label=='Stroma'))
new_tree

tree = new_tree

branch_col=data[,c('level2','level5')]
branch_col = branch_col[!duplicated(branch_col),]
xx=branch_col$level2
names(xx)=branch_col$level5

celltype_tree <- ggtree(tree, branch.length = 'none')

show_name = c()
levels = c()
for(i in celltype_tree$data$label){
    tmp_x = strsplit(i,':')[[1]]
    show_name = c(show_name, tmp_x[length(tmp_x)])
    levels = c(levels,tmp_x[1])
}

celltype_tree$data$show_name=show_name
celltype_tree$data$levels=levels
celltype_tree$data$branch_group = xx[celltype_tree$data$label]

celltype_tree$data$color_group = sapply(celltype_tree$data$show_name, function(x){
    xx = strsplit(x,'\\.')[[1]]
    as.numeric(xx[length(xx)])
})

celltype_tree$data

options(repr.plot.width=12, repr.plot.height=18)
# use ggtree visualize_tree
celltype_tree2=celltype_tree+aes(,color=levels)+
  geom_tiplab(aes(label = show_name,color=as.character(color_group)), size=label_size(6)) +
  #geom_text2(aes(label = show_name, subset = !isTip,color=show_name)) +
    scale_x_continuous(expand = c(0.5,0.1))+
    scale_color_manual(values=c('New'='black', 'Classical'='gray', 'r1:Imm'='#D7191C',
                                 'r1:Epi'='#1F78B4',
                                 'r1:EC'='#FEE08B',
                                 'r1:Str'='#A6761D', 
                                'r1'='red', 'r2'='black', 'r3' = 'blue', 'r4'='green',
                                round4_color2))+
    NoLegend()
celltype_tree2

label_odder = celltype_tree$data %>% filter(isTip) %>% arrange(y)

name_order = names(cluster_num[label_odder$color_group])
cell_prop_mat = all_umap %>%
    dplyr::group_by(Celltype_round4, Time) %>%
    dplyr::summarise(cell_num=n()) %>%
    reshape2::acast(Celltype_round4~Time, value.var='cell_num',fill = 0)

cell_prop_mat = t(t(cell_prop_mat)/colSums(cell_prop_mat))
cell_prop_mat = cell_prop_mat / rowSums(cell_prop_mat)

label_odder

name_order

output_dir

length(name_order)

write.csv(name_order, glue('{output_dir}/129cells_tree_order.txt'), row.names=F)

c=melt(cell_prop_mat) %>% 
    mutate(Celltype_round4=factor(Var1, levels = name_order)) %>%
    mutate(Time=factor(Var2, levels = rev(time_levels))) %>%
    ggplot(aes_string(x = 'Celltype_round4', y='value')) + 
    geom_bar(aes_string(fill = 'Time'), stat = "identity") + 
    scale_fill_manual(values = time_color2) + 
    theme_bw() + 
    coord_flip()+
    labs(title = "Time") + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(),
          axis.text.y = element_blank(),
         plot.title=element_text(hjust=0.5))
d=all_umap %>%
    dplyr::group_by(Celltype_round4,Celltype_round1) %>%
    dplyr::summarise(number=n()) %>%
    mutate(Celltype_round4=factor(Celltype_round4, name_order)) %>%
    ggplot(aes(x=Celltype_round4, y=log10(number),fill=Celltype_round4))+
    geom_bar(stat='identity', color='white', size=0.1)+
    scale_fill_manual(values = round4_color) + 
    theme_bw() + 
    coord_flip()+
    labs(title = "#Cell(log10)") + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          plot.title=element_text(hjust=0.5),
         legend.position='none')

options(repr.plot.width=16, repr.plot.height=12)
( celltype_tree2 | c | d ) + 
    plot_layout(guides='collect', widths=c(1,1,0.5))

a=( celltype_tree2 | c | d ) + 
    plot_layout(guides='collect', widths=c(2,1,1))

ggsave(glue('{output_dir}/subcluster_dendrogram_and_accessibility4.pdf'), a,
       width=4, height=8, units='in', dpi=600, bg='transparent')







umap = read.csv('./placeholder_analysis/round_cluster02/merge/cell_meta.csv')
sample_info = read.csv('./placeholder_data/sample_information/sample_information.csv',row.names = 1)
umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
rownames(umap) = umap$X

# usepython
# Figure1_gene_activityNMF

nmf_res = read.csv('./placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_meta_sample20wcells.csv')

head(nmf_res)


#keep_module = setdiff(1:16, c(1,2,7,10,11,13,14,15))
#keep_module

tmp_nmf_score = nmf_res[, paste0('M', 1:18)]
colnames(tmp_nmf_score) = paste0('M',1:ncol(tmp_nmf_score))

moudle_cell = apply(tmp_nmf_score,1,which.max)

moudle_cell = colnames(tmp_nmf_score)[moudle_cell]
names(moudle_cell) = nmf_res$X

nmf_res$cell_module = moudle_cell

nmf_umap = umap[umap$X%in%nmf_res$X, ]

nmf_umap$NMF_modules = moudle_cell[nmf_umap$X]

nmf_umap[, c('UMAP_1', 'UMAP_2')] = nmf_umap[, c('TSNE_1', 'TSNE_2')]

module_color = pal_igv()(19)
#module_color = module_color[c(12,16,3,4,5,6,8,9)]
names(module_color) = colnames(tmp_nmf_score)

write.csv(nmf_umap, './placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_meta_sample20wcells2.csv')
#nmf_umap = read.csv( './placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_meta_sample20wcells2.csv')

rownames(nmf_umap) = nmf_umap$X

head(nmf_umap)

table(nmf_umap$NMF_modules)

a = myDimPlot(nmf_umap, groupby='NMF_modules', group_color = module_color,label=FALSE, point_size=0.2)
ggsave(glue('{output_dir}/major_cluster_tsne_NMF_cells.png'), a+NoLegend(),
       width=60, height=60, units='mm', dpi=600, bg='transparent')

# TSNE Time
a = myDimPlot(nmf_umap, groupby='NMF_modules', group_color = module_color,label=FALSE, point_size=0.2)
ggsave(glue('{output_dir}/major_cluster_tsne_NMF_legend.pdf'), get_legend(a),
       width=150, height=150, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/major_cluster_tsne_NMF_cells.pdf'), a+NoLegend(),
       width=60, height=60, units='mm', dpi=600, bg='transparent')

nmf_umap = read.csv('./placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_meta_sample20wcells2.csv')

nmf_res  = read.csv('./placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_meta_sample20wcells.csv')

tmp_nmf_score = nmf_res[, paste0('M', 1:18)]
rownames(tmp_nmf_score) = nmf_res$X

nmf_umap[, colnames(tmp_nmf_score)] = tmp_nmf_score[nmf_umap$X,]

options(repr.plot.width=5, repr.plot.height=4)
score_list = list()
for(i in paste0('M', 1:18)){
    tmp_umap = nmf_umap[, c('UMAP_1', 'UMAP_2', i)]
    if(i=='M10'){
    q1 <- quantile(tmp_umap[,i], 0.15)
    q3 <- quantile(tmp_umap[,i], 0.85)
    }else{
    q1 <- quantile(tmp_umap[,i], 0.05)
    q3 <- quantile(tmp_umap[,i], 0.95)
    }

    tmp_umap[,i] <- pmin(pmax(tmp_umap[,i], q1), q3)
    a=tmp_umap %>%
        ggplot(aes(x=UMAP_1, y=UMAP_2))+
        ggrastr::geom_point_rast(aes_string(fill=i),shape=21,stroke=0,size=0.1, alpha=0.1)+
        scale_fill_gradientn(colours=c('#256da0','white','#bf2221'))+
        # scale_fill_gradient2(low='#256da0',mid='white', high = '#bf2221',
        #     limits = c(q1, q3),   # set_color_mapping_range
        #     oob = scales::squish # compress_out_of_range_values_to_bounds)
        #     )+
        labs(title=i)+
        theme_void()+
        theme(legend.key.width = unit(2,'mm'),
             legend.key.height = unit(2,'mm'),
             legend.position='none',
             plot.title=element_text(hjust=0.5, size=6)
             #legend.size=element_text(size=2)
             )
    score_list[[i]]=a
}

a = cowplot::plot_grid(plotlist = score_list, ncol=6)
ggsave(glue('{output_dir}/NMF_scores.pdf'), a,
       width=180, height=90, units='mm', dpi=600, bg='transparent')

options(repr.plot.width=5, repr.plot.height=4)
score_list = list()
for(i in paste0('M',c(13,14,11,8,18,15,16,9,17,5,1,6,4,10,12,3,2,7))){
    tmp_umap = nmf_umap[, c('UMAP_1', 'UMAP_2', i)]
    if(i=='M10'){
    q1 <- quantile(tmp_umap[,i], 0.15)
    q3 <- quantile(tmp_umap[,i], 0.85)
    tmp_umap[,i] <- pmin(pmax(tmp_umap[,i], q1), q3)
    }
    #
    a=tmp_umap %>%
        ggplot(aes(x=UMAP_1, y=UMAP_2))+
        ggrastr::geom_point_rast(aes_string(fill=i),shape=21,stroke=0,size=0.1, alpha=0.1)+
        scale_fill_gradientn(colours=c('#256da0','white','#bf2221'))+
        # scale_fill_gradient2(low='#256da0',mid='white', high = '#bf2221',
        #     limits = c(q1, q3),   # set_color_mapping_range
        #     oob = scales::squish # compress_out_of_range_values_to_bounds)
        #     )+
        labs(title=i)+
        theme_void()+
        theme(legend.key.width = unit(2,'mm'),
             legend.key.height = unit(2,'mm'),
             legend.position='none',
             plot.title=element_text(hjust=0.5, size=6)
             #legend.size=element_text(size=2)
             )
    score_list[[i]]=a
}

a = cowplot::plot_grid(plotlist = score_list, ncol=6)
ggsave(glue('{output_dir}/NMF_scores_raw.pdf'), a,
       width=180, height=90, units='mm', dpi=600, bg='transparent')


nmf_res_gene = read.csv('./placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_diffgene_sample20wcells.csv')

rownames(nmf_res_gene) = nmf_res_gene$Gene
head(nmf_res_gene)

table(nmf_res_gene$module)

cor_mat_z = read.csv('./placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_localcor_z_sample20wcells.csv', row.names=1, check.names=F)

gene_module = read.csv('./placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_modules_sample20wcells.csv', row.names=1)
gene_module = gene_module[gene_module$Module!=(-1),,drop=F]

table(gene_module$Module)

cor_mat_z[1:3,1:3]

cor_mat_z2 = cor_mat_z[rownames(gene_module),rownames(gene_module)]

tmp_m = paste0('M',gene_module$Module)
ra = rowAnnotation(df=data.frame(modules=tmp_m))

dim(cor_mat_z2)

tmp_rate = matrix(0.1, nrow = length(gene_module$Module), ncol=length(gene_module$Module))

for(i in unique(gene_module$Module)){
    tmp_pos = gene_module$Module==i
    tmp_rate[tmp_pos, tmp_pos] = 1
}

after_rate_mat= cor_mat_z2*tmp_rate

dim(cor_mat_z2)

module_order = c(18,8,11,15,1,5,9,17,13,14,16,10,4,6,3,12,2,7)
module_order_str = paste0('M',module_order)

options(repr.plot.width=6, repr.plot.height=8)

ht = Heatmap(cor_mat_z2,#[1:1000,1:1000],#[get_order(o,1),],
             cluster_rows=F, cluster_columns=F,
             show_row_dend = FALSE, show_column_dend = FALSE,
             show_row_names=F, show_column_names=F,
            #  col=circlize::colorRamp2(c(-10,0,10, 20), c(   "#ffffff",
            #                                                 "#d08f17",
            #                                                 "#a52f00",
            #                                                 "#620000")),
            col=circlize::colorRamp2(c(-60,-50,0,50,60), c("#91BFDB","#E0F3F8",  "#FFFFFF", "#FC8D59", "#D73027")),                                             
             right_annotation=ra,
             row_split = factor(gene_module$Module,module_order),
             column_split = factor(gene_module$Module,module_order),
             row_gap = unit(0, "mm"),           # newly_added:row_group_spacing0
             column_gap = unit(0, "mm"),          # newly_added:column_group_spacing0
             use_raster=T,
             #heatmap_legend_param=list('title'='%Cell\n(z-score)')
             #column_names_rot=0 
)

pdf(glue('{output_dir}/NMF_genecorrelation2.pdf'), width = 12,height = 12)
draw(ht)
dev.off()



options(repr.plot.width=6, repr.plot.height=8)

ht = Heatmap(cor_mat_z2,#[get_order(o,1),],
             cluster_rows=F, cluster_columns=F,
             show_row_dend = FALSE, show_column_dend = FALSE,
             show_row_names=F, show_column_names=F,
             col=circlize::colorRamp2(c(-50,-25,0,25,50), c("#91BFDB","#E0F3F8",  "#FFFFFF", "#FC8D59", "#D73027")),
             right_annotation=ra,
             row_split = factor(gene_module$Module,module_order),
             column_split = factor(gene_module$Module,module_order),
             row_gap = unit(0, "mm"),           # newly_added:row_group_spacing0
             column_gap = unit(0, "mm"),          # newly_added:column_group_spacing0
             use_raster=T,
             #heatmap_legend_param=list('title'='%Cell\n(z-score)')
             #column_names_rot=0 
)

pdf(glue('{output_dir}/NMF_genecorrelation.pdf'), width = 12,height = 12)
draw(ht)
dev.off()

nmf_gene_list = list()
for(i in 1:18){
    tmp_module = nmf_res_gene[nmf_res_gene$module==i, ]
    tmp_module = tmp_module %>% arrange(Z) 
    nmf_gene_list[[i]] = head(tmp_module$Gene, 200)
}

library(clusterProfiler)
library(org.Hs.eg.db)
mm10_hg38_homoGene = read.csv('./placeholder_root/mouse_lung_project/genome_annotation_file/mm10_hg38_homoGene_ensembl.csv')

head(mm10_hg38_homoGene)

all_module_go = c()
for(i in 1:18){
    tmp_gene = nmf_gene_list[[i]]
    tmp_gene_to_hm = mm10_hg38_homoGene[mm10_hg38_homoGene$Gene.name%in%tmp_gene,]
    tmp_gene_id = bitr(tmp_gene_to_hm$Human.gene.name, fromType = 'SYMBOL',
                      toType='ENTREZID', OrgDb='org.Hs.eg.db')
    tmp_go = enrichGO(
        gene=unique(tmp_gene_id$ENTREZID),
        keyType='ENTREZID',
        OrgDb='org.Hs.eg.db',
        ont='ALL',
        readable=TRUE,
    )
    tmp_res = tmp_go@result %>% as.data.frame()
    if(nrow(tmp_res)>0){
        tmp_res$module = i
        all_module_go = rbind(all_module_go, tmp_res)
    }
}

tmp_module = nmf_res_gene[nmf_res_gene$module==2, ]
tmp_module = tmp_module %>% arrange(Z) 
tmp_gene = head(tmp_module$Gene, 500)
tmp_gene_to_hm = mm10_hg38_homoGene[mm10_hg38_homoGene$Gene.name%in%tmp_gene,]
tmp_gene_id = bitr(tmp_gene_to_hm$Human.gene.name, fromType = 'SYMBOL',
                  toType='ENTREZID', OrgDb='org.Hs.eg.db')
tmp_go = enrichGO(
    gene=unique(tmp_gene_id$ENTREZID),
    keyType='ENTREZID',
    OrgDb='org.Hs.eg.db',
    ont='ALL',
    readable=TRUE,
)
tmp_res = tmp_go@result %>% as.data.frame()

tmp_res$module = 2
all_module_go = rbind(all_module_go, tmp_res)

write.table(all_module_go, '../share_go.txt', sep='\t', row.names = F)

all_module_go = read.table('../share_go.txt', sep='\t', header = T)

dim(all_module_go)

table(all_module_go$module)

head(all_module_go)

all_module_go_sig2

options(repr.plot.width=12, repr.plot.height=8)
all_module_go_sig = all_module_go %>% filter(p.adjust<0.05,Count>=5)

all_module_go_sig2 = all_module_go %>% filter(p.adjust<0.05,Count>=3, module==15)

all_module_go_sig = rbind(all_module_go_sig,all_module_go_sig2)

all_module_go_sig = all_module_go_sig %>%
    group_by(module) %>% 
    top_n(3, FoldEnrichment)

a=all_module_go_sig %>%
    mutate(module=factor(module, module_order)) %>%
    arrange(module) %>%
    mutate(Description=factor(Description, levels=Description)) %>%
    mutate(module = paste0('M',module)) %>%
    mutate(module = factor(module,rev(module_order_str))) %>%
    ggplot(aes(y=module, x=Description, size=-log10(p.adjust), fill=FoldEnrichment))+
    geom_point(shape=21)+
    scale_fill_gradientn(colours = c('white', '#bf2221'))+
    #scale_y_reverse()+
    theme_bw()+
    theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1))

 a   

ggsave(glue('{output_dir}/NMF_module_func.pdf'), a,
       width=230, height=150, units='mm', dpi=600, bg='transparent')

write.csv(all_module_go_sig, glue('./placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_module_GO_sig.csv'),
     row.names=F)

dim(all_module_go_sig)

all_module_go_sig = read.csv('./placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_module_GO_sig.csv')

dim(all_module_go_sig)

options(repr.plot.width=12, repr.plot.height=8)
a=all_module_go_sig %>%
    mutate(Description=abbr) %>%
    arrange(module) %>%
    mutate(Description=factor(Description, levels=Description)) %>%
    mutate(module = paste0('M',module)) %>%
    mutate(module = factor(module,paste0('M',18:1))) %>%
    ggplot(aes(y=module, x=Description, size=-log10(p.adjust), fill=FoldEnrichment))+
    geom_point(shape=21)+
    #geom_text_repel(aes(label=Description), size=3, box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"),)+
    scale_fill_gradientn(colours = c('white', '#bf2221'))+
    #scale_y_reverse()+
    theme_bw()+
    theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1))

 a   

module_accross = c()
for(i in paste0('M', 1:18)){
    tmp_umap = nmf_umap[, c('Celltype_round1', 'new_name', 'UMAP_2', i)]
    tmp_umap = tmp_umap[tmp_umap[,i]>0,]
    tmp_umap = tmp_umap %>%
        group_by(Celltype_round1,new_name) %>%
        summarise(cell_num=n())
    tmp_umap$module = i
    module_accross = rbind(module_accross, tmp_umap)
}

a=module_accross %>%
    group_by(module,Celltype_round1) %>%
    summarise(cell_num=sum(cell_num)) %>%
    mutate(module = factor(module, paste0('M',c(13,14,11,8,18,15,16,9,17,5,1,
                                                6,4,10,12,3,2,7)))) %>%
    ggplot(aes(x=module, y=cell_num, fill=Celltype_round1))+
    geom_bar(stat='identity', position='fill', color='white', size=0.1)+
    scale_fill_manual(values = celltype_round1_color)
a

 ggsave(glue('{output_dir}/NMF_module_cell_prop.pdf'), a,
       width=180, height=80, units='mm', dpi=600, bg='transparent')

tmp_umap = nmf_umap[, c('Celltype_round1', 'new_name', paste0('M', 1:18))]
tmp_umap = melt(tmp_umap, id.vars=c('Celltype_round1', 'new_name'))
tmp_umap = tmp_umap %>%
    group_by(Celltype_round1,new_name, variable) %>%
    summarise(value=mean(value))
plot_list = list()
for(i in paste0('M',c(13,14,11,8,18,15,16,9,17,5,1,6,4,10,12,3,2,7))){
    tmp_umap2 = tmp_umap[tmp_umap$variable==i,]
    a=tmp_umap2 %>% 
        as.data.frame()%>%
        mutate(new_name=fct_reorder(new_name, value)) %>%
        ggplot(aes(x=new_name,y=value, fill=Celltype_round1))+
        geom_bar(stat='identity')+
        labs(title = i, x='Round4_cluster', y='Mean module score')+
        scale_fill_manual(values = celltype_round1_color)+
        theme_bw()+
        theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.title=element_text(hjust=0.5),
            plot.background = element_rect(fill = 'white', color = NA),
            panel.grid = element_blank(),
            axis.title = element_text(siz=6)
            )
    plot_list[[i]]=a+NoLegend()
}  

ggsave(glue('{output_dir}/NMF_module_cluster.pdf'), cowplot::plot_grid(plotlist = plot_list, ncol=6),
       width=230, height=120, units='mm', dpi=600, bg='transparent')

tmp_umap = nmf_umap[, c('Celltype_round1', 'new_name', paste0('M', 1:18))]
tmp_umap = melt(tmp_umap, id.vars=c('Celltype_round1', 'new_name'))
tmp_umap = tmp_umap %>%
    group_by(Celltype_round1,new_name, variable) %>%
    summarise(value=mean(value))
plot_list = list()
for(i in paste0('M',c(13,14,11,8,18,15,16,9,17,5,1,6,4,10,12,3,2,7))){
    tmp_umap2 = tmp_umap[tmp_umap$variable==i,]
    a=tmp_umap2 %>% 
        as.data.frame()%>%
        mutate(new_name=fct_reorder(new_name, value)) %>%
        ggplot(aes(x=new_name,y=value, fill=Celltype_round1))+
        geom_bar(stat='identity', color=NA)+
        labs(title = i, x='Round4_cluster', y='Mean module score')+
        scale_fill_manual(values = celltype_round1_color)+
        theme_bw()+
        theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            plot.title=element_text(hjust=0.5),
            plot.background = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            )
    plot_list[[i]]=a+NoLegend()

}  
ggsave(glue('{output_dir}/NMF_module_cluster2.pdf'), cowplot::plot_grid(plotlist = plot_list, ncol=6),
       width=230, height=60, units='mm', dpi=600, bg='transparent')

options(repr.plot.width=6, repr.plot.height=12)
mudule_order = paste0('M',c(13,14,11,8,18,15,16,9,17,5,1,6,4,10,12,3,2,7))

xx = all_module_go_sig %>%
    mutate(Description=abbr) %>%
    mutate(module = paste0('M',module)) %>%
    mutate(module = factor(module,mudule_order)) %>%
    arrange(module, FoldEnrichment)
Description_order = rev(xx$Description)
a=all_module_go_sig %>%
    mutate(Description=as.vector(abbr)) %>%
    mutate(module = paste0('M',module)) %>%
    mutate(module = factor(module,mudule_order)) %>%
    arrange(module) %>%
    mutate(Description=factor(Description, levels=Description_order)) %>%
    ggplot(aes(x=module, y=Description, size=-log10(p.adjust), fill=FoldEnrichment))+
    geom_point(shape=21)+
    #geom_text(aes(label=Description), size=3, vjust=1,hjust=1 , angle=90)+
    scale_fill_gradientn(colours = c('white', '#bf2221'))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1),
     panel.grid = element_blank())

 a   
 ggsave(glue('{output_dir}/NMF_module_func2.pdf'), a,
       width=150, height=230, units='mm', dpi=600, bg='transparent')
 

xx = all_module_go_sig %>%
    mutate(Description=as.vector(go_abbr[Description,'Abbreviation'])) %>%
    mutate(module = paste0('M',module)) %>%
    mutate(module = factor(module,mudule_order)) %>%
    arrange(module)
rev(xx$Description)

name_order

#mudule_order = paste0('M',c(13,14,11,8,18,15,16,9,17,5,1,6,4,10,12,3,2,7))

Nek10_module_score = nmf_umap %>%
    group_by(Time) %>%
    dplyr::summarise(M15=mean(M15), M1=mean(M1), M5=mean(M5), M9=mean(M9), M17=mean(M17))# %>%
    #filter(Celltype_round4%in%c('AT1_Nek10'))
Nek10_module_score

nmf_umap2 = nmf_umap[, c(module_order_str, 'Celltype_round4')]

nmf_umap2 = nmf_umap[, c(module_order_str, 'Celltype_round4')]
nmf_umap2 = melt(nmf_umap2, id.vars='Celltype_round4')

head(nmf_umap2)

nmf_umap2= nmf_umap2 %>%
    left_join(nmf_umap2 %>%
        group_by(Celltype_round4, variable) %>%
        dplyr::summarise(mean_value=mean(value)) %>%
        mutate(mean_value_scale=scales::rescale(mean_value))
        )

head(nmf_umap2)

nmf_umap2 %>%
    filter(variable %in% c('M13','M14')) %>%
    mutate(Celltype_round4=factor(Celltype_round4, name_order)) %>%
    group_by(Celltype_round4, variable) %>%
    mutate(mean_value = scales::rescale(mean_value)) %>% head()



aa = nmf_umap2 %>%
    #filter(variable %in% c('M13','M14')) %>%
    mutate(Celltype_round4=factor(Celltype_round4, name_order)) %>%
    ggplot(aes(y=Celltype_round4, x=value,fill=mean_value_scale))+
    geom_violin(color='black',size=0.1, scale='width')+
    #scale_fill_manual(values = round4_color) + 
    scale_fill_gradientn(colours = c("#91BFDB","#E0F3F8",  "#FFFFFF", "#FC8D59", "#D73027"))+
    facet_wrap(~variable, nrow=1, scales = 'free_x')+
    theme_bw() + 
    labs(title = "#Cell(log10)") + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          plot.title=element_text(hjust=0.5),
         legend.position='none')
#aa

options(repr.plot.width=16, repr.plot.height=12)
( celltype_tree2 | aa ) + 
    plot_layout(guides='collect', widths=c(1,10))

ggsave(glue('{output_dir}/NMF_module_cluster3.pdf'), ( celltype_tree2 | aa ) + 
       plot_layout(guides='collect', widths=c(1,10)),
       width=16, height=12, units='in', dpi=600, bg='transparent')



#representativepeakplot

gene_mat2 = t(gene_mat)

xx = table(nmf_umap$Celltype_round4, nmf_umap$NMF_modules)

xx = xx[rownames(gene_mat2),]

xx = t(apply(xx, 1, function(x)x/sum(x)))

ra = rowAnnotation(modules=anno_barplot(xx, gp=gpar(fill=module_color),
                                       width=unit(20,'mm'), bar_width=1
                                       ))

options(repr.plot.width=6, repr.plot.height=8)

ht = Heatmap(gene_mat2,#[get_order(o,1),],
             cluster_rows=F, cluster_columns=F,
             show_row_names=F, show_column_names=F,
             col=circlize::colorRamp2(c(-2,-1,0,1,2), c("#91BFDB","#E0F3F8",  "#FFFFFF", "#FC8D59", "#D73027")),
             right_annotation=ra
             #heatmap_legend_param=list('title'='%Cell\n(z-score)')
             #column_names_rot=0 
)

pdf(glue('{output_dir}/NMFplots_and_accessibility.pdf'), width = 6,height = 8)
draw(ht)
dev.off()

module_gene = read.csv('./placeholder_root/mouse_lung_project/pycode//round_cluster02/merge/NMF_diffgene_sample20wcells.csv')

module_gene_filter = module_gene %>% filter(module!=(-1)) 

# new_name = paste0('M', 1:18)
#names(new_name) = as.character(new_name)

module_gene_filter$module_name = paste0('M', module_gene_filter$module)#new_name#[as.character(module_gene_filter$module)]

module_gene_filter2 = module_gene_filter %>% 
    filter(Gene%in%Annotation(pseudoObj_cells)$gene_name) %>%
    filter(FDR<0.05) %>%
    filter(!grepl('Rik?',Gene)) %>%
    filter(!grepl('^Gm',Gene)) %>%
    dplyr::group_by(module_name) %>%
    dplyr::top_n(Z, n=2)

module_gene_filter2

umap = read.csv('./placeholder_analysis/round_cluster02/merge/cell_meta.csv')
sample_info = read.csv('./placeholder_data/sample_information/sample_information.csv',row.names = 1)
umap$sample_id = factor(sample_info[umap$sample, 'sample_time'], levels=sample_levels)
umap$Time = factor(sample_info[umap$sample, 'Time'], levels=time_levels)
rownames(umap) = umap$X

# new_celltype = read.csv('./placeholder_analysis/round_cluster02/merge/celltype_modify.csv')
# rownames(new_celltype) = new_celltype$Celltype_round4
# new_celltype$new_name = paste0(new_celltype$MainType,'_', new_celltype$SubType, '.', new_celltype$cluster_num)

# new_celltype$new_name2 = paste0(new_celltype$MainType,'_', new_celltype$SubType)

new_celltype = read.csv('./placeholder_analysis/round_cluster02/merge/celltype_modify.csv')
rownames(new_celltype) = new_celltype$Celltype_round4

new_celltype$modify_name = new_celltype$Celltype_round4_new
# new_celltype$modify_name = gsub('^Myeloid_', '',new_celltype$modify_name)
# new_celltype$modify_name = gsub('^Stroma_', '',new_celltype$modify_name)
# new_celltype$modify_name = gsub('^Fibroblast_', 'Fibro_',new_celltype$modify_name)
# new_celltype$modify_name = gsub('^Monocytes_', 'Mono_',new_celltype$modify_name)
# new_celltype$modify_name = gsub('^Macorphage_', 'Mac_',new_celltype$modify_name)
# 
# new_celltype$modify_name = gsub('^Myofibroblast_', 'Myofibro_',new_celltype$modify_name)
# new_celltype$modify_name = gsub('^EC_', '',new_celltype$modify_name)
# new_celltype$modify_name = gsub('^Airway_', '',new_celltype$modify_name)

new_celltype$new_name = paste0(new_celltype$modify_name,'.', new_celltype$cluster_num)

500/775

round4_color_new = round4_color

names(round4_color_new) = new_celltype[names(round4_color), 'new_name']

round4_color

all_cell_prop = c()
for(i in unique(umap$Celltype_round1)){
    tmp_umap = umap[umap$Celltype_round1==i,]
    tmp_count_mat = tmp_umap %>% 
        dplyr::group_by(Celltype_round4, Time) %>%
        dplyr::summarise(n=n()) %>%
        reshape2::acast(Celltype_round4~Time, value.var = 'n', fill=0)
    tmp_prop = t(t(tmp_count_mat)/colSums(tmp_count_mat))
    tmp_prop = melt(tmp_prop)
    colnames(tmp_prop) = c('Celltype_round4', 'Time', 'prop')
    tmp_prop$Celltype_round1 = i
    all_cell_prop = rbind(all_cell_prop, tmp_prop)
}

all_cell_prop$new_name = new_celltype[as.vector(all_cell_prop$Celltype_round4), 'new_name']

umap$new_name = new_celltype[as.vector(umap$Celltype_round4), 'new_name']

umap$new_name2 = new_celltype[as.vector(umap$Celltype_round4), 'new_name2']

library(seriation)

register_DendSer()

cell_prop_mat = umap %>%
  dplyr::group_by(Time, Celltype_round4) %>%
  dplyr::summarise(cell_num=n()) %>%
  reshape2::acast(Time~Celltype_round4, fill=0)
cell_prop_mat = cell_prop_mat / rowSums(cell_prop_mat)
tmp_mat = t(apply(t(cell_prop_mat),1,scale))

cell_prop_mat = reshape2::acast(all_cell_prop,Time~new_name, value.var = 'prop')

cell_prop_mat = cell_prop_mat / rowSums(cell_prop_mat)
tmp_mat = t(apply(t(cell_prop_mat),1,scale))
#tmp_mat = t(apply(t(cell_prop_mat),1,function(x)(x-min(x))/(max(x)-min(x))))

min(cell_prop_mat)

#tmp_mat = t(apply(cell_prop_mat,2, function(x)(x-min(x))/(max(x)-min(x))))

predict_seq = seq(1,15,1)
tmp_mat_smmoth = t(apply(tmp_mat, 1, function(x){
    tmp_lm = loess(formula = y~x, data=data.frame(x=1:15,y=x))
    predict(tmp_lm, predict_seq)
}))
colnames(tmp_mat_smmoth) = sapply(predict_seq, function(x)ifelse(x%in%1:15, x, ''))
o = seriate(tmp_mat_smmoth,method='Heatmap', 
            #margin=1,
            seriation_method = 'DendSer_BAR',dist_fun=dist)

colnames(tmp_mat_smmoth) = paste0('P', as.numeric(colnames(tmp_mat_smmoth))-1)

tmp_meta = umap[,c('Celltype_round1', 'new_name')] %>% unique()
rownames(tmp_meta) = tmp_meta$new_name

leftanno = rowAnnotation(
    df=data.frame(
        Celltype_round1=tmp_meta[rownames(tmp_mat_smmoth),'Celltype_round1']
    ),
    col=list(
        Celltype_round1=celltype_round1_color
    ),
    show_legend = TRUE,
    show_annotation_name = FALSE
)

options(repr.plot.width=6, repr.plot.height=8)

ht = Heatmap(tmp_mat_smmoth,#[get_order(o,1),],
             cluster_rows=stats::as.dendrogram(o[[1]]),
             left_annotation=leftanno, 
             cluster_columns=F,
             show_row_names=F,
             show_column_dend=F,
             show_row_dend = F,
             col=circlize::colorRamp2(c(-2,-1,0,1,2), c("#E0F3F8", "#91BFDB", "#FFFFFF", "#FC8D59", "#D73027")),
             heatmap_legend_param=list('title'='%Cell\n(z-score)')
             #column_names_rot=0 
)
ht

pdf(glue('{output_dir}/cell_proportion_heatmap15.pdf'), width = 6,height = 8)
draw(ht)
dev.off()

saveRDS(list('smooth'=tmp_mat_smmoth,
             'raw'=t(cell_prop_mat),
             'cell_order'=rownames(tmp_mat_smmoth)[row_order(ht)]),
       '../out_data/cell_prop_time_mat.rds')

predict_seq = seq(1,15,0.01)
tmp_mat_smmoth = t(apply(tmp_mat, 1, function(x){
    tmp_lm = loess(formula = y~x, data=data.frame(x=1:15,y=x))
    predict(tmp_lm, predict_seq)
}))
colnames(tmp_mat_smmoth) = sapply(predict_seq, function(x)ifelse(x%in%1:15, x, ''))
#o = seriate(tmp_mat_smmoth,method='Heatmap', 
#            margin=1,
#            seriation_method = 'DendSer_BAR',dist_fun=dist)

options(repr.plot.width=6, repr.plot.height=8)

ht = Heatmap(tmp_mat_smmoth,#[get_order(o,1),],
             cluster_rows=stats::as.dendrogram(o[[1]]), 
             cluster_columns=F,
             show_row_names=F,
             show_column_dend=F,
             col=circlize::colorRamp2(c(-2,-1,0,1,2), c("#E0F3F8", "#91BFDB", "#FFFFFF", "#FC8D59", "#D73027")),
             #column_names_rot=0 
)
ht

pdf(glue('{output_dir}/cell_proportion_heatmap150.pdf'), width = 6,height = 8)
draw(ht)
dev.off()

la = data.frame(stage = cutree(o[[1]], 3))
la$stage = paste0('C',la$stage)
la_col = pal_igv()(20)[1:length(unique(la$stage))]
names(la_col) = unique(la$stage)
la = rowAnnotation(df=la,
                  col=list('stage'=la_col))


options(repr.plot.width=6, repr.plot.height=8)

ht = Heatmap(tmp_mat_smmoth,#[get_order(o,1),],
             cluster_rows=stats::as.dendrogram(o[[1]]), 
             cluster_columns=F,
             show_row_names=F,
             show_column_dend=F,
             show_row_dend=F,
             col=circlize::colorRamp2(c(-2,-1,0,1,2), c("#E0F3F8", "#91BFDB", "#FFFFFF", "#FC8D59", "#D73027")),
             left_annotation=la
)
draw(ht)

draw(ht)

pdf(glue('{output_dir}/cell_proportion_heatmap150(cluster).pdf'), width = 6,height = 8)
draw(ht)
dev.off()

la = data.frame(stage = cutree(o[[1]], 3))
la$stage = paste0('C',la$stage)

all_prop = melt(tmp_mat_smmoth) 
colnames(all_prop) = c('Celltype', 'Time', 'prop')
all_prop$cluster = la[all_prop$Celltype, 'stage']
#all_prop$Time = all_prop$Time - 1

data_summary2 <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  #print(data_sum)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

all_prop_error <- data_summary2(all_prop, varname="prop", 
                    groupnames=c("Time", 'cluster'))


a=all_prop_error %>%
    mutate(cluster=factor(cluster, c('C1', 'C3', 'C2'))) %>%
  ggplot(aes(x=Time, y=prop))+
  geom_line(aes(group=1), size=0.1)+
  geom_errorbar(aes(x=Time,ymin=prop-sd,ymax=prop+sd),data=all_prop_error, width=0.3,size=0.2,
                inherit.aes = F)+
      facet_wrap(~cluster, ncol=1,strip.position = 'right')+
  geom_point(aes(fill=Time), shape=21, size=0.8, strokes=0.1)+
  scale_x_discrete(breaks =c('P0', 'P4', 'P12', 'P14'))+
  scale_fill_manual(values = time_color2)+
   mytheme+
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_line(color='black', size=0.1, linetype = 'dashed'))
a

ggsave(glue('{output_dir}/cell_proportion_heatmapclusterline_plot.pdf'), a,
       width=45, height=85, units='mm', dpi=600, bg='transparent')

new_celltype

tmp_meta = umap[,c('Celltype_round1', 'new_name')] %>% unique()
rownames(tmp_meta) = tmp_meta$new_name

all_prop$Celltype_round1 = tmp_meta[as.vector(all_prop$Celltype), 'Celltype_round1']
head(all_prop)

options(repr.plot.width=12, repr.plot.height=4)
a=all_prop %>%
    mutate(cluster=factor(cluster, c('C1', 'C3', 'C2'))) %>%
    ggplot(aes(x=Time, y=prop))+
    geom_line(aes(group=Celltype, color=Celltype_round1), linewidth=0.5, alpha=0.53)+
    geom_smooth(method='loess', color='black', se=F)+
    facet_wrap(~cluster, ncol=3)+
    scale_color_manual(values = celltype_round1_color)+
    #scale_color_manual(values = c('C2'='#313695', 'C1'='#74ADD1', 'C3'='#FEE090', 'C4'='#F46D43'))+
    #scale_x_discrete(breaks = c('P4', 'P2'), labels = c('P0', 'P3', 'P7', 'P13'))+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          #panel.grid.minor.y = element_blank(),
          #panel.grid.minor.x=element_blank(),
          #panel.grid.major.x=element_line(color='black', size=0.1),
          #panel.grid.major.y
         strip.background = element_blank(),
         panel.border=element_rect(size=0.1, fill=NA))
a

ggsave(glue('{output_dir}/cell_proportion_heatmapclusterline_plot(2).pdf'), a,
       width=240, height=70, units='mm', dpi=600, bg='transparent')

umap_meta = umap[, c('Celltype_round1','new_name')] %>% unique()
rownames(umap_meta) = umap_meta$new_name

la$Celltype_round1 = umap_meta[rownames(la), 'Celltype_round1']

xx = table(la$stage, la$Celltype_round1)
xx = xx / rowSums(xx)
xx_ypos = melt(apply(xx,1, cumsum))
colnames(xx_ypos) = c('Celltype_round1', 'stage', 'ypos')

xx = melt(xx)


colnames(xx) = c('stage', 'Celltype_round1', 'prop')
xx$label = paste0(round(xx$prop,4)*100, '%')

xx = left_join(xx, xx_ypos)

xx

options(repr.plot.width=8, repr.plot.height=4)
a=xx %>% 
    mutate(stage=factor(stage, c('C1', 'C3', 'C2'))) %>%
    ggplot(aes(x="", y=prop, fill=Celltype_round1))+
    geom_bar(width = 1, stat = "identity")+
    geom_text(aes(y=1-ypos+0.1, label=label))+
    coord_polar("y", start=0)+
    facet_wrap("~ stage",ncol = 1, scale = 'free') +
    theme_void()+
    scale_fill_manual(name='No.clusters', values=celltype_round1_color)+
    #labs(title = 'Non-overlap cCREs') +
    theme( plot.title = element_text(hjust = 0.5))
a

ggsave(glue('{output_dir}/cell_proportion_heatmapclusterpie_chart.pdf'), a+NoLegend(),
       width=60, height=100, units='mm', dpi=600, bg='transparent')

options(repr.plot.width=8, repr.plot.height=8)
sample_cor = cor(tmp_mat)
sample_cor[lower.tri(sample_cor,diag = F)] = NA

ht = Heatmap(sample_cor, 
        #col = circlize::colorRamp2(c(0,1),c("white","#a50f15")),
        cluster_rows = F, cluster_columns = F,show_column_names = F,
        na_col = 'white',
        width = ncol(sample_cor)*unit(2, "mm"), # width
        height = nrow(sample_cor)*unit(2, "mm"), # height
        #row_names_gp = gpar(col=time_color,fontsize=12, font='bold'),
        heatmap_legend_param = list(title='Cor')
        )
ht

pdf(glue('{output_dir}/cell_proportion_time_similarity.pdf'), width = 3,height = 3)
draw(ht)
dev.off()

la

# between_subtypespropsimilarity
sample_cor = cor(t(tmp_mat), method='spearman')


options(repr.plot.width=14, repr.plot.height=14)
ht=Heatmap(sample_cor, 
        col = circlize::colorRamp2(c(-1,0,1),c("#045a8d","white","#a50f15")),
        cluster_rows = T, cluster_columns =T,show_column_names = T,
        #row_split=r4_relation[rownames(sample_cor), 'Celltype_round1'],
        #column_split=r4_relation[rownames(sample_cor), 'Celltype_round1'],
        na_col = 'white',
        #width = ncol(sample_cor)*unit(2, "mm"), # width
        #height = nrow(sample_cor)*unit(2, "mm"), # height
        row_names_gp = gpar(fontsize=6, font='bold'),
        column_names_gp = gpar(fontsize=6, font='bold'),
        heatmap_legend_param = list(title='Cor'),
        right_ann=rowAnnotation(df=la,
                               col=list('stage'=c('C1'='#637abb','C2'='#b2d688','C3'='#ea5142'), 
                                        'Celltype_round1'=celltype_round1_color))
)
ht

tmp_mat

draw(ht)

pdf(glue('{output_dir}/cell_proportion_similarity_heatmap(clustering).pdf'), width = 12,height = 10)
draw(ht)
dev.off()

# between_subtypespropsimilarity
sample_cor = cor(t(tmp_mat), method='spearman')


options(repr.plot.width=14, repr.plot.height=14)
ht=Heatmap(sample_cor, 
        col = circlize::colorRamp2(c(-1,0,1),c("#045a8d","white","#a50f15")),
        cluster_rows = T, cluster_columns =T,show_column_names = T,
        row_split=la$stage,
        column_split=la$stage,
        na_col = 'white',
        #width = ncol(sample_cor)*unit(2, "mm"), # width
        #height = nrow(sample_cor)*unit(2, "mm"), # height
        row_names_gp = gpar(fontsize=6, font='bold'),
        column_names_gp = gpar(fontsize=6, font='bold'),
        heatmap_legend_param = list(title='Cor'),
        right_ann=rowAnnotation(df=la,
                               col=list('stage'=c('C1'='#637abb','C2'='#b2d688','C3'='#ea5142'), 
                                        'Celltype_round1'=celltype_round1_color))
)
ht

pdf(glue('{output_dir}/cell_proportion_similarity_heatmap(stagegrouping).pdf'), width = 12,height = 10)
draw(ht)
dev.off()

orig_df = as.data.frame.matrix(table(umap[,'Time'], umap[,'Celltype_round1']))
orig_df = orig_df / rowSums(orig_df)
orig_df$x = 1:15
orig_df$Time=rownames(orig_df)
orig_df = melt(orig_df, id.vars=c('Time','x'))

options(repr.plot.width=6, repr.plot.height=4)
a=orig_df %>%
    mutate(Time=factor(Time, levels = time_levels)) %>%
    ggplot()+
    #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
    geom_point(mapping = aes(x=Time, y=value,color=variable),  size=0.5,fill=NA,shape=21, stroke=0.1)+
    stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
    geom_smooth(mapping = aes(x=x, y=value,color=variable),
                method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=0.5)+
    labs(y='Cell / Total')+
    scale_color_manual(values = c('Immune'='#D7191C',
                 'Epithelium'='#1F78B4',
                 'Endothelium'='#FEE08B',
                 'Stroma'='#A6761D'))+
    mytheme+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
a

ggsave(glue('{output_dir}/cell_proportion_changes_major_subclusters.pdf'), a+NoLegend(),
       width=40, height=40, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/cell_proportion_changes_major_subclusters_legend.pdf'), get_legend(a),
       width=40, height=40, units='mm', dpi=600, bg='transparent')

Immune_prop = orig_df %>% filter(variable=='Immune')

options(repr.plot.width=16, repr.plot.height=4)
a1=distribution_func(umap, 'Time', 'Celltype_round1',c('Immune'='#D7191C',
                 'Epithelium'='#1F78B4',
                 'Endothelium'='#FEE08B',
                 'Stroma'='#A6761D'))
a2=distribution_func(umap, 'Time', 'Celltype_round1',c('Immune'='#D7191C',
                 'Epithelium'='#1F78B4',
                 'Endothelium'='#FEE08B',
                 'Stroma'='#A6761D'), 'dodge')
a1+a2

ggsave(glue('{output_dir}/cell_proportion_changes_major_subclustersbar.pdf'), a1+a2,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

options(repr.plot.width=16, repr.plot.height=4)
res=distribution_Roe_func(umap, 'Celltype_round1','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)

dev.off()
pdf(glue('{output_dir}/cell_proportion_changes_major_subclusters_Roe.pdf'), width = 8,height = 2)
res=distribution_Roe_func(umap, 'Celltype_round1','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()

round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')
round4_num = readRDS('./placeholder_analysis/round_cluster02/merge/round4_cluster_num.rds')
round4_color2 = round4_color
names(round4_color2) = round4_num[names(round4_color2)]

round4_color_new = round4_color

names(round4_color_new) = new_celltype[names(round4_color), 'new_name']

tmp_umap = umap[umap$Celltype_round1=='Epithelium',]
tmp_umap = tmp_umap%>% filter(!Celltype_round3%in%c('Tip', 'SMG_Oit1', 'SMG_Rnase1'))
tmp_round='Celltype_round2'

options(repr.plot.width=16, repr.plot.height=4)
#a=distribution_func(tmp_umap, 'Time', tmp_round)
#a

a1=distribution_func(tmp_umap, 'Time', tmp_round)+labs(y='Cell proportion')
a2=distribution_func(tmp_umap, 'Time', tmp_round, , position = 'dodge')+labs(y='Cell number')
a1+a2
ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}.pdf'), a1+a2,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

#ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}.pdf'), a,
#       width=8, height=6, units='in', dpi=600, bg='transparent')

options(jupyter.plot_mimetypes = "image/png")

tmp_umap2 = tmp_umap%>% filter(!Celltype_round3%in%c('Tip', 'SMG_Oit1', 'SMG_Rnase1'))
res=distribution_Roe_func(tmp_umap2, tmp_round,'Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)

dev.off()
pdf(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}_Roe.pdf'), width = 8,height = 2)
res=distribution_Roe_func(tmp_umap2, tmp_round,'Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()

orig_df = as.data.frame.matrix(table(tmp_umap[,'Time'], tmp_umap[,tmp_round]))
orig_df = orig_df / rowSums(orig_df)
orig_df$x = 1:15
orig_df$Time=rownames(orig_df)
orig_df = melt(orig_df, id.vars=c('Time','x'))


options(repr.plot.width=10, repr.plot.height=8)
# a=orig_df %>%
#     mutate(Time=factor(Time, levels = time_levels)) %>%
#     ggplot()+
#     #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
#     geom_point(mapping = aes(x=Time, y=value,color=variable),  size=2)+
#     stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
#     geom_smooth(mapping = aes(x=x, y=value,color=variable),
#                 method = 'glm' , formula = 'y ~ poly(x,3)', se=F,  size=1)+
#     scale_color_manual(values = expanded_colors2)+
#     labs(y='Cell / Epithelial')+
#     theme_classic()+
#     mytheme
# a
a=orig_df %>%
    mutate(Time=factor(Time, levels = time_levels)) %>%
    ggplot()+
    #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
    geom_point(mapping = aes(x=Time, y=value,color=variable),  size=0.5,fill=NA,shape=21, stroke=0.1)+
    stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
    geom_smooth(mapping = aes(x=x, y=value,color=variable),
                method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=0.5)+
    labs(y='Cell / Epithelial')+
    scale_color_manual(values = expanded_colors2)+
    mytheme+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
a
ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}_line_plot.pdf'), a+NoLegend(),
       width=40, height=40, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}_line_plot_legend.pdf'), get_legend(a),
       width=40, height=40, units='mm', dpi=600, bg='transparent')

# ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}_line_plot.pdf'), a,
#        width=70, height=50, units='mm', dpi=600, bg='transparent')

tmp_umap = umap[umap$Celltype_round1=='Epithelium',]
tmp_umap = tmp_umap%>% filter(!Celltype_round3%in%c('Tip', 'SMG_Oit1', 'SMG_Rnase1'))
tmp_round='new_name'

options(repr.plot.width=8, repr.plot.height=6)
a=distribution_func(tmp_umap, 'Time', tmp_round, color=round4_color_new)
a

ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}.pdf'), a,
       width=8, height=6, units='in', dpi=600, bg='transparent')

tmp_umap2 = tmp_umap%>% filter(!Celltype_round3%in%c('Tip', 'SMG_Oit1', 'SMG_Rnase1'))
res=distribution_Roe_func(tmp_umap2, tmp_round,'Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)

dev.off()
pdf(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}_Roe.pdf'), width = 8,height = 6)
res=distribution_Roe_func(tmp_umap2, tmp_round,'Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()

tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='AT1')
orig_df = as.data.frame.matrix(table(tmp_umap1[,'Time'], tmp_umap1[,'new_name']))
orig_df = orig_df / rowSums(orig_df)
orig_df$x = 1:15
orig_df$Time=rownames(orig_df)
orig_df = melt(orig_df, id.vars=c('Time','x'))


options(repr.plot.width=10, repr.plot.height=8)
# a=orig_df %>%
#     mutate(Time=factor(Time, levels = time_levels)) %>%
#     ggplot()+
#     #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
#     geom_point(mapping = aes(x=Time, y=value,color=variable),  size=2)+
#     stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
#     geom_smooth(mapping = aes(x=x, y=value,color=variable),
#                 method = 'glm' , formula = 'y ~ poly(x,3)', se=F,  size=1)+
#     scale_color_manual(values = expanded_colors2)+
#     labs(y='Cell / Epithelial')+
#     theme_classic()+
#     mytheme
# a
a=orig_df %>%
    mutate(Time=factor(Time, levels = time_levels)) %>%
    ggplot()+
    #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
    geom_point(mapping = aes(x=Time, y=value,color=variable),  size=0.5,fill=NA,shape=21, stroke=0.1)+
    stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
    geom_smooth(mapping = aes(x=x, y=value,color=variable),
                method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=0.5)+
    labs(y='Cell / Epithelial')+
    scale_color_manual(values = expanded_colors2)+
    mytheme+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
a
# ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}_line_plot.pdf'), a+NoLegend(),
#        width=40, height=40, units='mm', dpi=600, bg='transparent')
# ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_{tmp_round}_line_plot_legend.pdf'), get_legend(a),
#        width=40, height=40, units='mm', dpi=600, bg='transparent')

# tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='AT1')
# orig_df = melt(cell_prop_mat[,unique(tmp_umap1$new_name)])
# colnames(orig_df) = c('Time', 'Celltype', 'value')
# orig_df$x = as.numeric(gsub('P','',orig_df$Time))+1

tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='AT1')
orig_df = as.data.frame.matrix(table(tmp_umap1[,'Time'], tmp_umap1[,'new_name']))
orig_df = orig_df / rowSums(orig_df)
orig_df$x = 1:15
orig_df$Time=rownames(orig_df)
orig_df = melt(orig_df, id.vars=c('Time','x'))

options(repr.plot.width=10, repr.plot.height=8)
# a=orig_df %>%
#     mutate(Time=factor(Time, levels = time_levels)) %>%
#     ggplot()+
#     #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
#     geom_point(mapping = aes(x=Time, y=value,color=Celltype),  size=1)+
#     stat_cor(mapping = aes(x=x, y=value,color=Celltype), size=label_size(6))+
#     geom_smooth(mapping = aes(x=x, y=value,color=Celltype),
#                 method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=1)+
#     scale_color_manual(values = round4_color_new)+
#     labs(y='Cell / AT1')+
#     theme_classic()+
#     mytheme
# a

a=orig_df %>%
    mutate(Time=factor(Time, levels = time_levels)) %>%
    ggplot()+
    #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
    geom_point(mapping = aes(x=Time, y=value,color=variable),  size=0.5,fill=NA,shape=21, stroke=0.1)+
    stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
    geom_smooth(mapping = aes(x=x, y=value,color=variable),
                method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=0.5)+
    labs(y='Cell / AT1')+
    scale_color_manual(values = round4_color_new)+
    mytheme+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
a
ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_AT1_major_subclusters.pdf'), a+NoLegend(),
       width=40, height=40, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_AT1_major_subclusters_legend.pdf'), get_legend(a),
       width=40, height=40, units='mm', dpi=600, bg='transparent')

nek10_prop = orig_df %>% filter(variable%in%c('AT1_Nek10.112'))
mfrp_prop = orig_df %>% filter(variable%in%c('AT1_Mfrp.103'))

Nek10_module_score=Nek10_module_score %>% 
    mutate(Time=factor(Time, time_levels)) %>%
    arrange(Time)
Nek10_module_score

all_merge_df = cbind(nek10_prop,mfrp_prop$value, Immune_prop$value)
all_merge_df = cbind(all_merge_df,Nek10_module_score[, c('M15', 'M1', 'M5', 'M9', 'M17')])
colnames(all_merge_df) = c('Time', 'x', 'Nek10', 'Nek10_prop', 'Mfrp_prop', 'Immune_prop', 'M15', 'M1', 'M5', 'M9', 'M17')
all_merge_df = melt(all_merge_df, id.vars = c('Time', 'x', 'Nek10', 'Nek10_prop'))
all_merge_df %>% head()

options(repr.plot.width=16, repr.plot.height=3)
a=all_merge_df %>%
    ggplot(aes(x=Nek10_prop, y=value))+
    geom_point(aes(, color=Time))+
    geom_smooth(se=F, method='lm', color='black')+
    stat_cor()+
    facet_wrap(~variable, scales = 'free_y', nrow=1)+
    scale_color_manual(values = time_color2)+
    theme_bw()+
    theme(panel.grid = element_blank(), strip.background = element_blank())
a
ggsave(glue('{output_dir}/cell_proportion_changes_Nek10_Immune_cor.pdf'), a,
       width=280, height=50, units='mm', dpi=600, bg='transparent')

ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_AT1_major_subclusters.pdf'), a,
       width=80, height=50, units='mm', dpi=600, bg='transparent')

res=distribution_Roe_func(tmp_umap1, 'new_name','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)

dev.off()
pdf(glue('{output_dir}/cell_proportion_changes_Epithelium_AT1_Roe.pdf'), width = 8,height = 3)
res=distribution_Roe_func(tmp_umap1, 'new_name','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()

options(repr.plot.width=16, repr.plot.height=4)
#a=distribution_func(tmp_umap, 'Time', tmp_round)
#a

a1=distribution_func(tmp_umap1, 'Time', 'new_name')+labs(y='Cell proportion')
a2=distribution_func(tmp_umap1, 'Time', 'new_name', , position = 'dodge')+labs(y='Cell number')
a1+a2
ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_AT1.pdf'), a1+a2,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='AT2')

table(tmp_umap1$new_name)
xx=table(tmp_umap1$new_name)/nrow(tmp_umap1)
xx
xx[1]+xx[3]

late_tmp_umap1 = tmp_umap1 %>% filter(Time%in%c('P12','P13','P14'))
table(late_tmp_umap1$new_name)
xx=table(late_tmp_umap1$new_name)/nrow(late_tmp_umap1)
xx

late_tmp_umap1 = tmp_umap1 %>% filter(Time%in%c('P14'))
table(late_tmp_umap1$new_name)
xx=table(late_tmp_umap1$new_name)/nrow(late_tmp_umap1)
xx

# tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='AT2')
# orig_df = melt(cell_prop_mat[,unique(tmp_umap1$new_name)])
# colnames(orig_df) = c('Time', 'Celltype', 'value')
# orig_df$x = as.numeric(gsub('P','',orig_df$Time))+1
tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='AT2')
orig_df = as.data.frame.matrix(table(tmp_umap1[,'Time'], tmp_umap1[,'new_name']))
orig_df = orig_df / rowSums(orig_df)
orig_df$x = 1:15
orig_df$Time=rownames(orig_df)
orig_df = melt(orig_df, id.vars=c('Time','x'))

options(repr.plot.width=10, repr.plot.height=8)
# a=orig_df %>%
#     mutate(Time=factor(Time, levels = time_levels)) %>%
#     ggplot()+
#     #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
#     geom_point(mapping = aes(x=Time, y=value,color=Celltype),  size=1)+
#     stat_cor(mapping = aes(x=x, y=value,color=Celltype), size=label_size(6))+
#     geom_smooth(mapping = aes(x=x, y=value,color=Celltype),
#                 method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=1)+
#     scale_color_manual(values = round4_color_new)+
#     labs(y='Cell / AT2')+
#     theme_classic()+
#     mytheme
# a
a=orig_df %>%
    mutate(Time=factor(Time, levels = time_levels)) %>%
    ggplot()+
    #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
    geom_point(mapping = aes(x=Time, y=value,color=variable),  size=0.5,fill=NA,shape=21, stroke=0.1)+
    stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
    geom_smooth(mapping = aes(x=x, y=value,color=variable),
                method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=0.5)+
    labs(y='Cell / AT2')+
    scale_color_manual(values = round4_color_new)+
    mytheme+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
a
ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_AT2_major_subclusters.pdf'), a+NoLegend(),
       width=40, height=40, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_AT2_major_subclusters_legend.pdf'), get_legend(a),
       width=40, height=40, units='mm', dpi=600, bg='transparent')

res=distribution_Roe_func(tmp_umap1, 'new_name','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)

dev.off()
pdf(glue('{output_dir}/cell_proportion_changes_Epithelium_AT2_Roe.pdf'), width = 8,height = 2)
res=distribution_Roe_func(tmp_umap1, 'new_name','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()

options(repr.plot.width=16, repr.plot.height=4)
#a=distribution_func(tmp_umap, 'Time', tmp_round)
#a

a1=distribution_func(tmp_umap1, 'Time', 'new_name')+labs(y='Cell proportion')
a2=distribution_func(tmp_umap1, 'Time', 'new_name', , position = 'dodge')+labs(y='Cell number')
a1+a2
ggsave(glue('{output_dir}/cell_proportion_changes_Epithelium_AT2.pdf'), a1+a2,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

tmp_umap = umap[umap$Celltype_round1=='Immune',]
#tmp_umap[tmp_umap$]
tmp_round='Celltype_round2'
tmp_umap[tmp_umap$Celltype_round3%in%c('NK'), tmp_round]='NK'
#tmp_umap[tmp_umap$Celltype_round3%in%c('NKT'), tmp_round]='NKT'

options(repr.plot.width=12, repr.plot.height=6)
#a=distribution_func(tmp_umap, 'Time', tmp_round)
#a

options(repr.plot.width=16, repr.plot.height=4)
#a=distribution_func(tmp_umap, 'Time', tmp_round)
#a

a1=distribution_func(tmp_umap, 'Time', tmp_round)+labs(y='Cell proportion')
a2=distribution_func(tmp_umap, 'Time', tmp_round, , position = 'dodge')+labs(y='Cell number')
a1+a2
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_{tmp_round}.pdf'), a1+a2,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

tmp_umap2 = tmp_umap%>% filter(!Celltype_round3%in%c('Tip', 'SMG_Oit1', 'SMG_Rnase1'))
res=distribution_Roe_func(tmp_umap2, tmp_round,'Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)

dev.off()
pdf(glue('{output_dir}/cell_proportion_changes_Immune_{tmp_round}_Roe.pdf'), width = 8,height = 2)
res=distribution_Roe_func(tmp_umap2, tmp_round,'Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()

orig_df = as.data.frame.matrix(table(tmp_umap[,'Time'], tmp_umap[,tmp_round]))
orig_df = orig_df / rowSums(orig_df)
orig_df$x = 1:15
orig_df$Time=rownames(orig_df)
orig_df = melt(orig_df, id.vars=c('Time','x'))

options(repr.plot.width=10, repr.plot.height=8)
# a=orig_df %>%
#     mutate(Time=factor(Time, levels = time_levels)) %>%
#     ggplot()+
#     #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
#     geom_point(mapping = aes(x=Time, y=value,color=variable),  size=2)+
#     stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
#     geom_smooth(mapping = aes(x=x, y=value,color=variable),
#                 method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=1)+
#     scale_color_npg()+
#     labs(y='Cell / Immune')+
#     theme_classic()+
#     mytheme
# a

a=orig_df %>%
    mutate(Time=factor(Time, levels = time_levels)) %>%
    ggplot()+
    #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
    geom_point(mapping = aes(x=Time, y=value,color=variable),  size=0.5,fill=NA,shape=21, stroke=0.1)+
    stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
    geom_smooth(mapping = aes(x=x, y=value,color=variable),
                method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=0.5)+
    labs(y='Cell / Immune')+
    scale_color_npg()+
    mytheme+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
a
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_{tmp_round}_line_plot.pdf'), a+NoLegend(),
       width=40, height=40, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_{tmp_round}_line_plot_legend.pdf'), get_legend(a),
       width=40, height=40, units='mm', dpi=600, bg='transparent')

tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='Myeloid') 
orig_df = as.data.frame.matrix(table(tmp_umap1[,'Time'], tmp_umap1[,'Celltype_round3']))
orig_df = orig_df / rowSums(orig_df)
orig_df$x = 1:15
orig_df$Time=rownames(orig_df)
orig_df = melt(orig_df, id.vars=c('Time','x'))

options(repr.plot.width=10, repr.plot.height=8)
# a=orig_df %>%
#     filter(variable%in%c('Macorphage', 'Monocytes', 'Neutrophils'))%>%
#     mutate(Time=factor(Time, levels = time_levels)) %>%
#     ggplot()+
#     #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
#     geom_point(mapping = aes(x=Time, y=value,color=variable),  size=2)+
#     stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
#     geom_smooth(mapping = aes(x=x, y=value,color=variable),
#                 method = 'glm' , formula = 'y ~ poly(x,3)', se=F,  size=1)+
#     scale_color_manual(values = expanded_colors2)+
#     labs(y='Cell / Myeloid')+
#     theme_classic()+
#     mytheme
# a

a=orig_df %>%
    filter(variable%in%c('Macorphage', 'Monocytes', 'Neutrophils'))%>%
    mutate(Time=factor(Time, levels = time_levels)) %>%
    ggplot()+
    #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
    geom_point(mapping = aes(x=Time, y=value,color=variable),  size=0.5,fill=NA,shape=21, stroke=0.1)+
    stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
    geom_smooth(mapping = aes(x=x, y=value,color=variable),
                method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=0.5)+
    labs(y='Cell / Myeloid')+
    scale_color_manual(values = expanded_colors2)+
    mytheme+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
a
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_round3_myeloid_line_plot.pdf'), a+NoLegend(),
       width=40, height=40, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_round3_myeloid_line_plot_legend.pdf'), get_legend(a),
       width=40, height=40, units='mm', dpi=600, bg='transparent')

a1=distribution_func(tmp_umap1%>%filter(Celltype_round3%in%c('Macorphage', 'Monocytes', 'Neutrophils')), 'Time', 'Celltype_round3')+labs(y='Cell proportion')
a2=distribution_func(tmp_umap1%>%filter(Celltype_round3%in%c('Macorphage', 'Monocytes', 'Neutrophils')), 'Time', 'Celltype_round3', , position = 'dodge')+labs(y='Cell number')
a1+a2
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_Myeloid_{tmp_round}.pdf'), a1+a2,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

tmp_umap2 = tmp_umap%>% filter(Celltype_round3%in%c('Macorphage', 'Monocytes', 'Neutrophils'))
res=distribution_Roe_func(tmp_umap2, 'Celltype_round3','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)

dev.off()
pdf(glue('{output_dir}/cell_proportion_changes_Immune_Myeloid_{tmp_round}_Roe.pdf'), width = 8,height = 2)
res=distribution_Roe_func(tmp_umap2, 'Celltype_round3','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()

tmp_umap1 = tmp_umap %>% 
    filter(Celltype_round2=='Tcell') %>%
    filter(!Celltype_round4%in%c('ILC1','ILC2', 'ILC3') )
tmp_umap1[tmp_umap1$Celltype_round3=='Treg', 'Celltype_round3'] = 'CD4'
tmp_umap1[tmp_umap1$Celltype_round3%in%c('NKT'), 'Celltype_round3']='NKT'

#tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='Tcell') 
orig_df = as.data.frame.matrix(table(tmp_umap1[,'Time'], tmp_umap1[,'Celltype_round3']))
orig_df = orig_df / rowSums(orig_df)
orig_df$x = 1:15
orig_df$Time=rownames(orig_df)
orig_df = melt(orig_df, id.vars=c('Time','x'))

options(repr.plot.width=10, repr.plot.height=8)
# a=orig_df %>%
#     filter(!variable%in%c('ILCs'))%>%
#     mutate(Time=factor(Time, levels = time_levels)) %>%
#     ggplot()+
#     #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
#     geom_point(mapping = aes(x=Time, y=value,color=variable),  size=2)+
#     stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
#     geom_smooth(mapping = aes(x=x, y=value,color=variable),
#                 method = 'glm' , formula = 'y ~ poly(x,3)', se=F,  size=1)+
#     scale_color_manual(values = expanded_colors2)+
#     labs(y='Cell / Tcell')+
#     theme_classic()+
#     mytheme
# a
a=orig_df %>%
    mutate(Time=factor(Time, levels = time_levels)) %>%
    ggplot()+
    #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
    geom_point(mapping = aes(x=Time, y=value,color=variable),  size=0.5,fill=NA,shape=21, stroke=0.1)+
    stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
    geom_smooth(mapping = aes(x=x, y=value,color=variable),
                method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=0.5)+
    labs(y='Cell / Tcell')+
    scale_color_manual(values = expanded_colors2)+
    mytheme+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
a
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_round3_Tcell_line_plot.pdf'), a+NoLegend(),
       width=40, height=40, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_round3_Tcell_line_plot_legend.pdf'), get_legend(a),
       width=40, height=40, units='mm', dpi=600, bg='transparent')

a1=distribution_func(tmp_umap1, 'Time', 'Celltype_round3')+labs(y='Cell proportion')
a2=distribution_func(tmp_umap1, 'Time', 'Celltype_round3', , position = 'dodge')+labs(y='Cell number')
a1+a2
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_Tcell_{tmp_round}.pdf'), a1+a2,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

res=distribution_Roe_func(tmp_umap1, 'Celltype_round3','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()
pdf(glue('{output_dir}/cell_proportion_changes_Immune_Tcell_{tmp_round}_Roe.pdf'), width = 8,height = 2)
res=distribution_Roe_func(tmp_umap1, 'Celltype_round3','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()

tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='Bcell') 
tmp_umap1[tmp_umap1$Celltype_round3%in%c('MatureB_(Cd5+)Ccl22','MatureB_(Cd5+)Dusp2','MatureB(Cd5-)'), 
          'Celltype_round3'] = 'MatureB'

table(tmp_umap1$Celltype_round3)

#tmp_umap1 = tmp_umap %>% filter(Celltype_round2=='Bcell') 
orig_df = as.data.frame.matrix(table(tmp_umap1[,'Time'], tmp_umap1[,'Celltype_round3']))
orig_df = orig_df / rowSums(orig_df)
orig_df$x = 1:15
orig_df$Time=rownames(orig_df)
orig_df = melt(orig_df, id.vars=c('Time','x'))

options(repr.plot.width=10, repr.plot.height=8)
# a=orig_df %>%
#     filter(!variable%in%c('ILCs'))%>%
#     mutate(Time=factor(Time, levels = time_levels)) %>%
#     ggplot()+
#     #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
#     geom_point(mapping = aes(x=Time, y=value,color=variable),  size=2)+
#     stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
#     geom_smooth(mapping = aes(x=x, y=value,color=variable),
#                 method = 'glm' , formula = 'y ~ poly(x,3)', se=F,  size=1)+
#     scale_color_manual(values = expanded_colors2)+
#     labs(y='Cell / Bcell')+
#     theme_classic()+
#     mytheme
# a
a=orig_df %>%
    mutate(Time=factor(Time, levels = time_levels)) %>%
    ggplot()+
    #geom_line(mapping = aes(x=Time, y=value,color=variable,group=variable), size=2)+
    geom_point(mapping = aes(x=Time, y=value,color=variable),  size=0.5,fill=NA,shape=21, stroke=0.1)+
    stat_cor(mapping = aes(x=x, y=value,color=variable), size=label_size(6))+
    geom_smooth(mapping = aes(x=x, y=value,color=variable),
                method = 'glm' , formula = 'y ~ poly(x,3)', se=F, size=0.5)+
    labs(y='Cell / Tcell')+
    scale_color_manual(values = expanded_colors2)+
    mytheme+
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
a
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_round3_Bcell_line_plot.pdf'), a+NoLegend(),
       width=40, height=40, units='mm', dpi=600, bg='transparent')
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_round3_Bcell_line_plot_legend.pdf'), get_legend(a),
       width=40, height=40, units='mm', dpi=600, bg='transparent')

a1=distribution_func(tmp_umap1, 'Time', 'Celltype_round3')+labs(y='Cell proportion')
a2=distribution_func(tmp_umap1, 'Time', 'Celltype_round3', , position = 'dodge')+labs(y='Cell number')
a1+a2
ggsave(glue('{output_dir}/cell_proportion_changes_Immune_Bcell_{tmp_round}.pdf'), a1+a2,
       width=210, height=40, units='mm', dpi=600, bg='transparent')

res=distribution_Roe_func(tmp_umap1, 'Celltype_round3','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()
pdf(glue('{output_dir}/cell_proportion_changes_Immune_Bcell_{tmp_round}_Roe.pdf'), width = 8,height = 2)
res=distribution_Roe_func(tmp_umap1, 'Celltype_round3','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)
dev.off()

length(umap$new_name %>% unique())

tmp_umap = umap[umap$Celltype_round1=='Immune',]



options(repr.plot.width=12, repr.plot.height=4)
tmp_umap2 = tmp_umap%>% filter(Celltype_round3%in%c('Macorphage'))
tmp_umap2$Celltype_round4 = gsub('Macorphage_','Macrophage_', tmp_umap2$Celltype_round4)
res=distribution_Roe_func(tmp_umap2, 'Celltype_round4','Time', obj_is_meta=TRUE, 
                      show.sign = 'value', angle_col = 90,size=12)

