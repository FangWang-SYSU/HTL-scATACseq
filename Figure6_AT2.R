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

load_peak<-function(name){
    # Load differential peak matrix
    mat <- Matrix::readMM(glue("{base_path}/{name}_aggr_sparse.mtx"))
    row_names <- read.csv(glue("{base_path}/{name}_aggr_sparse_index.csv"), header = T)[,1]
    col_names <- read.csv(glue("{base_path}/{name}_aggr_sparse_columns.csv"), header = T)[,1]
    row_names = gsub(':','-',row_names)
    rownames(mat) <- row_names
    colnames(mat) <- col_names
    return(mat)
}

output_dir='./placeholder_output/raw_figure_output/Figure3_AT2/'

celltype = 'AT2'
act_assay='CRE' # 
base_path='./placeholder_analysis/round_cluster02/cCRE/'

srt = readRDS(glue('{base_path}/celltype_rds/Epi_sub.rds'))

srt = subset(srt, Celltype_round2=='AT2')
Idents(srt) = srt$Time
srt

srt = my_process_srt(srt, assay= act_assay)

srt_marker = read.csv(glue("{base_path}/celltype_rds/Epi_{celltype}_time_peaks.csv"))
colnames(srt_marker)[c(1,2,4,5,6)] = c('cluster','gene', 'avg_log2FC', 'p_val','p_val_adj')
srt_marker$gene = gsub(':','-',srt_marker$gene)

merge.peaks = read.csv(glue('./placeholder_analysis/round_cluster02/cCRE/all_merge_peaks_HMMann.csv'), check.names = F)

merge.peaks = merge.peaks[merge.peaks$Tissue=='lung', ]
merge.peaks$cre = gsub(':','-',merge.peaks$cre)
merge.peaks = merge.peaks[, c('cre', 'hmm_state', 'isin_encode')]

unique_time_cre = srt_marker %>% 
  dplyr::filter(avg_log2FC >=0.1&p_val_adj <=0.05) %>%
  pull(gene) %>% unique() 
unique_time_cre%>% length()

unique_time_cre_df = srt_marker %>% 
  dplyr::filter(avg_log2FC >=0.1&p_val_adj <=0.05)

unique_time_cre_df2 = merge(unique_time_cre_df,merge.peaks, by.x='gene', by.y='cre', all.x=TRUE)

head(unique_time_cre_df2)
dim(unique_time_cre_df2)

chromHMM_colors_simple <- c(
  "Pr-A" = "#006400",
  "Pr-W" = "#90EE90",
  "Pr-B" = "#D3D3D3",
  "Pr-F" = "#008000",
  "En-Sd" = "#FFFF00",
  "En-Sp" = "#FFFF00",
  "En-W" = "#FFFFE0",
  "En-Pd" = "#696969",
  "En-Pp" = "#696969",
  "Tr-S" = "#00008B",
  "Tr-P" = "#4169E1",
  "Tr-I" = "#ADD8E6",
  "Hc-P" = "#FA8072",
  "Hc-H" = "#FFC0CB",
  "NS"   = "#FFFFFF"
)
chromHMM_colors_simple2 = data.frame(chromHMM_colors_simple)
chromHMM_colors_simple2$state=rownames(chromHMM_colors_simple2)
chromHMM_colors_simple2$full_name = c(
    "Promoter-Active",
    "Promoter-Weak/Inactive",
    "Promoter-Bivalent",
    "Promoter-Flanking",
    "Enhancer-Strong TSS-distal",
    "Enhancer-Strong TSS-proximal",
    "Enhancer-Weak TSS-distal",
    "Enhancer-Poised TSS-distal",
    "Enhancer-Poised TSS-proximal",
    "Transcription-Strong",
    "Transcription-Permissive",
    "Transcription-Initiation",
    "Heterochromatin-Polycomb",
    "Heterochromatin-H3K9me3",
    "No Chromatin Signal"
)

chromHMM_colors_simple_color2 = chromHMM_colors_simple
names(chromHMM_colors_simple_color2) = chromHMM_colors_simple2$full_name

unique_time_cre_df2$hmm_state_name = chromHMM_colors_simple2[unique_time_cre_df2$hmm_state, 'full_name']

unique_time_cre_df2 %>%
    mutate(activated_state = ifelse(!hmm_state_name %in% c('No Chromatin Signal', 'Heterochromatin-H3K9me3', 'Heterochromatin-Polycomb'),
        'Activated','Other')) %>%
    dplyr::group_by(activated_state) %>%
    dplyr::summarise(n = n()) %>%
    mutate(freq = n/sum(n))

a=unique_time_cre_df2 %>%
    mutate(activated_state = ifelse(!hmm_state_name %in% c('No Chromatin Signal', 'Heterochromatin-H3K9me3', 'Heterochromatin-Polycomb'),
        'Activated','Other')) %>%
    dplyr::group_by(activated_state) %>%
    dplyr::summarise(n = n()) %>%
    ggplot(aes(x="", y=n, fill=activated_state))+
    geom_bar(stat='identity',width = 1, position = 'fill')+
    scale_fill_manual(values=c('Activated'='#FF7F0E','Other'='#1F77B4'))+
    coord_polar("y", start=0)+
    theme_void()
a
ggsave(glue('{output_dir}/AT2_time_cre_HMM_pie.pdf'), a,  
       width=210, height = 120, units='mm', dpi=600, bg='transparent')

a=unique_time_cre_df2 %>%
    mutate(cluster=factor(cluster, time_levels)) %>%
    dplyr::group_by(cluster, hmm_state_name) %>%
    dplyr::summarise(Num=n()) %>%
    ggplot(aes(x=cluster, y=Num, fill=hmm_state_name)) +
    geom_bar(stat='identity', position = 'fill', width=0.6, color='black', size=0.1)+
    scale_fill_manual(values = chromHMM_colors_simple_color2)+
    theme_bw()+
    theme(panel.grid = element_blank())
a
ggsave(glue('{output_dir}/AT2_time_cre_HMM.pdf'), a,  
       width=210, height = 120, units='mm', dpi=600, bg='transparent')

top_peaks = srt_marker %>% 
  dplyr::filter(avg_log2FC >=0.1&p_val_adj <=0.05) %>%
  group_by(cluster) %>% 
  top_n(2000,avg_log2FC) %>%
  arrange(cluster, -avg_log2FC)

top_peaks$cluster = factor(top_peaks$cluster, levels = time_levels)
top_peaks = top_peaks %>% arrange(cluster)
meanM = myAveragePeak(srt, 'Time', assay=act_assay, slot='data')
meanM = log2(meanM+1)
meanM = as.matrix(meanM)
meanM = myRowScale(meanM, max = 2, min = -2,limit=TRUE)


options(repr.plot.width=6, repr.plot.height=8)
ht = Heatmap(meanM[top_peaks$gene, time_levels],
        col = paletteContinuous(set = "solarExtra", n = 100),
        cluster_columns = F,cluster_rows = F,
        show_column_names = T, show_row_names = F,
        use_raster = T)
ht

pdf(glue('{output_dir}/AT2_differentialCREheatmap.pdf'), width = 6,height = 8)
draw(ht)
dev.off()

library(rGREAT)
Species = 'mm10'

tmp_group = list('stage1'=c('P0','P1','P2', 'P3'),
              'stage2'=c('P4','P5','P6', 'P7', 'P8', 'P9', 'P10','P11'),
              'stage3'=c('P12', 'P13', 'P14'))

Go.dataList = list()
for (i in names(tmp_group)) {
    peaks = srt_marker %>% 
                  dplyr::filter(avg_log2FC >=0.1&p_val_adj <=0.05) %>%
                  dplyr::filter(cluster%in%tmp_group[[i]])
    peaks = peaks$gene
    peaks = strsplit(peaks,split = '-')
    peak_dat = data.frame(chr = unlist(lapply(peaks, '[',1)),
                          start = unlist(lapply(peaks, '[',2)),
                          end =unlist(lapply(peaks, '[',3)))

    peak_dat = peak_dat[peak_dat$chr%in%paste0('chr',c(1:22,'X','Y')),]
    peak_gr = GRanges(seqnames = peak_dat$chr,
                      ranges = IRanges(start =as.numeric(peak_dat$start),
                                       end = as.numeric(peak_dat$end)))
    
    peak_gr = GenomicRanges::reduce(peak_gr)
    job = submitGreatJob(peak_gr,species =  Species )
    tbl = getEnrichmentTables(job)
    GO_BP = tbl$`GO Biological Process`
    GO_BP$Cluster = 'Biological Process'
    
    
    GO_MF = tbl$`GO Molecular Function`
    GO_MF$Cluster = 'Molecular Function'

    GO_CC =  tbl$`GO Cellular Component`
    GO_CC$Cluster = 'Cellular Component'

    Go.data = rbind(GO_BP,GO_CC,GO_MF)
    Go.data$Group = i 
    Go.dataList[[i]] = Go.data
}

saveRDS(Go.dataList, '../pycode/round_cluster02/cCRE/AT2_diff_cre_GO.rds')

Go.dataList = readRDS('../pycode/round_cluster02/cCRE/AT2_diff_cre_GO.rds')

positional_encoding <- function(time_vec, d_model = 16) {
  time_vec <- as.matrix(time_vec)
  N <- length(time_vec)
  i <- matrix(0:(d_model - 1), nrow = 1)

  angle_rates <- 1 / (10000 ^ ((2 * floor(i / 2)) / d_model))
  angle_matrix <- time_vec %*% angle_rates  # (N x d_model)

  pos_encoding <- matrix(0, nrow = N, ncol = d_model)
  even_idx <- seq(1, d_model, by = 2)
  odd_idx  <- seq(2, d_model, by = 2)

  pos_encoding[, even_idx] <- sin(angle_matrix[, even_idx])
  pos_encoding[, odd_idx]  <- cos(angle_matrix[, odd_idx])

  return(pos_encoding)
}
time_umap <- function(obj, time_vec,n.neighbors=30,dims=NULL,
                      reduction='lsi', time_weight=0.5, key='TimeLSI', assay=NULL,
                     time_coding=FALSE, time_coding_models = 30){
    if(is.null(assay)){
        assay = DefaultAssay(obj)
    }
    key_lower = tolower(key)
    X = Embeddings(obj, reduction)
    if(!is.null(dims)){
        dims = min(ncol(X), dims)
    }else{
        dims = 30
    }
    X = X[, 2:dims]
    
    if(time_coding){
        time_vec = positional_encoding(time_vec,time_coding_models)
    }
    
    New_X = cbind(X,time_vec*time_weight)
    colnames(New_X) = paste0(key,'_', 1:ncol(New_X))
    
    if(is.vector(time_vec))
    
    obj[[key_lower]] = CreateDimReducObject(New_X, key=paste0(key,'_'), assay = assay)
    obj <- RunUMAP(object = obj, 
                   reduction = key_lower,
                   dims = 1:ncol(New_X),
                   n.neighbors=n.neighbors,
                   metric = "euclidean",min.dist=1)
    return(obj)
}

time_numeric = as.numeric(gsub('P','', srt$Time))
srt[['Time_num']] = time_numeric

srt=time_umap(srt, srt$Time_num, dims=50, n.neighbors=30,
                reduction='lsi',time_weight=0.5)



srt <- FindNeighbors(object = srt, reduction = 'timelsi', dims = c(1:30,50))

srt <- FindClusters(object = srt, resolution = 0.2)

a = DimPlot(srt, group.by='seurat_clusters', label=TRUE,label.size=12)+scale_color_igv()
b = DimPlot(srt, group.by='Time', label=FALSE, cols=time_color2, raster = F)
c = DimPlot(srt, group.by='Celltype_round2', label=FALSE)+scale_color_npg()
d = DimPlot(srt, group.by='Celltype_round4', label=TRUE, raster = F)+scale_color_igv()
options(repr.plot.width=18, repr.plot.height=16)
a+b+c+d

round4_color = readRDS('./placeholder_analysis/round_cluster02/merge/round4_color.rds')

options(repr.plot.width=18, repr.plot.height=8)
b = DimPlot(srt, group.by='Time', label=FALSE, cols=time_color2, raster = F)
d = DimPlot(srt, group.by='Celltype_round4', label=TRUE, raster = F,cols=round4_color)
b+d

ggsave(glue('{output_dir}/AT2_UMAP.pdf'), b/d,  
       width=120, height = 180, units='mm', dpi=600, bg='transparent')

library(BSgenome.Mmusculus.UCSC.mm10)
jaspar <- JASPAR2024::JASPAR2024()
pfm <- getMatrixSet(
  x = jaspar@db,opts = list(collection = "CORE",tax_group='vertebrates' , all_versions = FALSE))

hg38_mm10_TF_mapping = read.table('../hg38_mm10_TF_mapping.txt', sep='\t', header = TRUE)
hg38_mm10_TF_mapping = hg38_mm10_TF_mapping[hg38_mm10_TF_mapping$Mouse_Symbol%in%Annotation(srt)$gene_name,]

keep_TF = unique(hg38_mm10_TF_mapping$tf_name)

keep_TF_upper = toupper(keep_TF)

dup_tf = keep_TF_upper[duplicated(keep_TF_upper)]

keep_TF = setdiff(keep_TF,dup_tf)

keep_TF_motifname = c()
for(i in pfm@listData){
    if(i@name %in% keep_TF){
        keep_TF_motifname = c(keep_TF_motifname,i@ID)
    }
}
pfm_filter = pfm[names(pfm)%in%keep_TF_motifname]

jaspar = JASPAR2024::JASPAR2024()
pfm = TFBSTools::getMatrixSet(x=jaspar@db, opts=list(all_version=T, species=9606, collection='CORE'))
tf_meta = c()
for(i in pfm@listData){
  tf_meta = rbind(tf_meta, c(i@ID,i@name,i@tags$family,i@tags$symbol))
}


tf_meta = as.data.frame(tf_meta)

if_names = tf_meta$V2

if_names2 = sapply(if_names, function(x)strsplit(x,' ')[[1]][1])
if_names2 = gsub('_extended', '',if_names2)

tf_meta = tf_meta[tf_meta$V2%in%if_names2,]

motif_meta = c()
for(m in pfm_filter@listData){
    motif_meta = rbind(motif_meta, c(m@ID, paste0(m@tags$family,collapse = '::'), m@name))
}
motif_meta = as.data.frame(motif_meta)
colnames(motif_meta) = c('motif', 'family', 'TF')
rownames(motif_meta) = motif_meta$motif

srt <- AddMotifs(object = srt,genome = BSgenome.Mmusculus.UCSC.mm10,pfm = pfm_filter,verbose = F)

srt = AddChromVARSelf(obj = srt, Species = 'mm10')

saveRDS(srt, './placeholder_analysis/round_cluster02/cCRE//celltype_rds/Epi_AT2_sub.rds')

srt = readRDS('./placeholder_analysis/round_cluster02/cCRE//celltype_rds/Epi_AT2_sub.rds')

top_peaks = srt_marker %>% 
  dplyr::filter(avg_log2FC >=0.1&p_val_adj <=0.05) %>%
  group_by(cluster) %>% 
  top_n(2000,avg_log2FC) %>%
  arrange(cluster, -avg_log2FC)

meta.feature <- GetAssayData(srt, assay = act_assay, layer = "meta.features")

Cluster = unique(top_peaks$cluster)
enrich.list <- lapply(Cluster,function(i){
    message(i)
     foreground.peaks = unique(subset(top_peaks,cluster == i)$gene)
     background.peaks = AccessiblePeaks(srt,assay = act_assay)
     set.seed(123)
     peaks.matched <- MatchRegionStats(meta.feature = meta.feature[background.peaks, ],
                                      query.feature = meta.feature[foreground.peaks, ],
                                      n = 50000)
     # background.peaks = AccessiblePeaks(object,assay = assay,idents = c(Cluster,i))
     enrich.motif <- FindMotifs(object = srt, 
                              features = foreground.peaks,
                              background = peaks.matched)
     enrich.motif$Cluster = i
     return(enrich.motif)
}) 
enrich.mat = do.call(rbind,enrich.list)

write.table(enrich.mat, glue('{output_dir}/AT2_Motif_Enrichment.txt'), quote = F, sep='\t', row.names = F)

enrich.mat = read.table(glue('{output_dir}/AT2_Motif_Enrichment.txt'), header = T, sep='\t')

enrich.mat_AT1 = read.table('./placeholder_output/raw_figure_output/Figure3//AT1_Motif_Enrichment.txt', header = T)
enrich.mat_AT1$Celltype='AT1'

enrich.mat_AT2 = enrich.mat
enrich.mat_AT2$Celltype='AT2'

enrich.mat_merge = rbind(enrich.mat_AT1, enrich.mat_AT2)

write.table(enrich.mat_merge, glue('{output_dir}/AT1_AT2_Motif_Enrichment.txt'), quote = F, sep='\t', row.names = F)

enrich.mat %>%
    filter(p.adjust<0.05) %>% pull(motif) %>% unique() %>% length()

plot.list=list()
for(i in time_levels){
    temp = subset(enrich.mat, Cluster == i )
    temp <- temp[order(temp$mlog10Padj, decreasing = TRUE),]
    temp$rank <- seq_len(nrow(temp))
    
    p = ggplot(temp, aes(rank, mlog10Padj, color = mlog10Padj)) + 
        geom_point(size = 2) +
        ggrepel::geom_text_repel(
            data = temp[rev(seq_len(10)), ], 
            aes(x = rank, y = mlog10Padj, label = motif.name), 
            size = label_size(6),
            nudge_x = 2,
            max.overlaps = 100,
            color = "black") +
        labs(title = i,
             color = "- log10(p.adjust)")+
        theme_few() + 
        ylab("-log10(P-adj) Motif Enrichment") + 
        xlab("Rank Sorted Motif Enriched") +
        scale_color_viridis(option = "viridis")+
        theme(plot.title = element_text(hjust = 0.5, size = 12))+NoLegend()
    plot.list[[i]]=p
}

ggsave(glue('{output_dir}/AT2_Motif_Enrichment_point_rank.pdf'),
       cowplot::plot_grid(plotlist = plot.list, ncol=5),
       width=320, height = 220, units='mm', dpi=600, bg='transparent')



plotlist = lapply(time_levels, function(i){
    motif_Clus = subset(enrich.mat,Cluster == i)
    motif_Clus <- motif_Clus[order(motif_Clus$mlog10Padj, decreasing = TRUE),]
    MotifPlot(object = srt,motifs = head(motif_Clus$motif,2)) + ggtitle(label = i)+
    theme_void()+
        theme(plot.title = element_text(hjust = 0.5, size = 6,colour = 'black'), 
        strip.text = element_text(hjust = 0.5, size = 6,colour = 'black'),
        plot.background = element_rect(fill=NA, color='black', size=0.1))
        
})

options(repr.plot.width=16, repr.plot.height=8)
cowplot::plot_grid(plotlist = plotlist, ncol=5)

ggsave(glue('{output_dir}/AT2_Motif_Enrichment_LOGO.pdf'),
       cowplot::plot_grid(plotlist = plotlist, ncol=3),
       width=58, height = 80, units='mm', dpi=600, bg='transparent')

meanVar = myAveragePeak(srt, 'Time', assay='chromvar', slot='data')
meanVar = as.matrix(meanVar)
meanVar = myRowScale(meanVar, max = 2, min = -2,limit=TRUE)

sig_motif = enrich.mat %>%
    filter(p.adjust<0.05) %>%
    arrange(-fold.enrichment)
sig_motif$Cluster = factor(sig_motif$Cluster, levels = time_levels)
sig_motif = sig_motif %>% arrange(Cluster)

#sig_gene = ConvertMotifID(srt,id = sig_motif$motif,assay = act_assay)
sig_motif$sig_gene = motif_meta[sig_motif$motif, 'TF']

sig_motif = sig_motif[!duplicated(sig_motif$motif),]

sig_motif = sig_motif %>% arrange(Cluster,-mlog10p)

ht_data = meanVar[sig_motif$motif, ]
#rownames(ht_data) = sig_motif$sig_gene

ht_data = ht_data[,order(as.numeric(gsub('P','',colnames(ht_data))))]
ht_data = t(apply(ht_data, 1, function(x){
    tmp_lm = loess(formula = y~x, data=data.frame(x=1:ncol(ht_data),y=x))
    #predict(tmp_lm, x=seq(1,15,0.1))
    predict(tmp_lm, newdata=data.frame(x=seq(1,ncol(ht_data),1)))
}))

options(repr.plot.width=6, repr.plot.height=8)
Heatmap(ht_data,
        col = paletteContinuous(set = "solarExtra", n = 100),
        cluster_columns = F,cluster_rows = T,
        show_column_names = T, show_row_names = F,
        use_raster = F)

meanVar_pseudo = myAveragePeak(srt, 'Time', assay='chromvar', slot='data')
meanVar_pseudo = as.matrix(meanVar_pseudo)
meanVar_pseudo = myRowScale(meanVar_pseudo, max = 2, min = -2,limit=TRUE)

meanVar_pseudo_na = na.omit(meanVar_pseudo)

ht_data = meanVar_pseudo[sig_motif$motif, ]
ht_data = ht_data[,order(as.numeric(gsub('P','',colnames(ht_data))))]

rownames(ht_data) = sig_motif$sig_gene

ht_data = t(apply(ht_data, 1, function(x){
    tmp_lm = loess(formula = y~x, data=data.frame(x=1:ncol(ht_data),y=x))
    #predict(tmp_lm, x=seq(1,15,0.1))
    predict(tmp_lm, newdata=data.frame(x=seq(1,ncol(ht_data),1)))
}))

options(repr.plot.width=8, repr.plot.height=8)
ht=Heatmap(ht_data,
        col = paletteContinuous(set = "solarExtra", n = 100),
        cluster_columns = F,cluster_rows = T,
        show_column_names = F, show_row_names = F,
           #row_split=sig_motif$Cluster,
        use_raster = F)
ht

cre_acces_pseudo = myAveragePeak(srt, 'Time', assay='CRE', slot='data')
cre_acces_pseudo = log2(cre_acces_pseudo+1)
cre_acces_pseudo = as.matrix(cre_acces_pseudo)
cre_acces_pseudo = myRowScale(cre_acces_pseudo, max = 2, min = -2,limit=TRUE)

pesutotime_levels = colnames(cre_acces_pseudo)[order(as.numeric(gsub('P','',colnames(cre_acces_pseudo))))]
col_name_levels = pesutotime_levels

tf_ht_order = row_order(ht)
tf_ht_order = sig_motif$sig_gene[tf_ht_order]

tf_ht_order_data = ht_data[row_order(ht),]
rownames(tf_ht_order_data) = tf_ht_order
colnames(tf_ht_order_data) = col_name_levels

options(repr.plot.width=8, repr.plot.height=8)
ht_tf=Heatmap(tf_ht_order_data,
        col = paletteContinuous(set = "solarExtra", n = 100),
        cluster_columns = F,cluster_rows = F,
        show_column_names = F, show_row_names = F,
        use_raster = F)
ht_tf

dim(tf_ht_order_data)

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

gb = get_genebody(object = EnsDb.Mmusculus.v79)
seqlevels(gb) = paste0('chr',seqlevels(gb))

keep_gb = gb[gb$gene_name%in%hg38_mm10_TF_mapping$Mouse_Symbol,]

gn = keep_gb$gene_name
names(gn) = GRangesToString(keep_gb)

gene_body_matrix=FeatureMatrix(fragments = myfragments, features = keep_gb)

saveRDS(gene_body_matrix, './placeholder_analysis//round_cluster02/merge/AT2_gene_matrix.rds')

head(gn)


gene_body_matrix = readRDS('./placeholder_analysis//round_cluster02/merge/AT2_gene_matrix.rds')

gene_body_matrix = gene_body_matrix[names(gn),]

rownames(gene_body_matrix) = gn

srt[["GeneACC"]]=NULL

srt[["GeneACC"]] = CreateAssayObject(gene_body_matrix,key = 'GeneACC')

meanACC_pseudo = myAveragePeak(srt, 'Time', assay="GeneACC", slot='data')
meanACC_pseudo = log2(meanACC_pseudo+1)
meanACC_pseudo = as.matrix(meanACC_pseudo)
meanACC_pseudo = myRowScale(meanACC_pseudo, max = 4, min = -4,limit=TRUE)

tmp_tf_acc = meanACC_pseudo[, col_name_levels]

tmp_tf_acc = t(apply(tmp_tf_acc, 1, function(x){
    tmp_lm = loess(formula = y~x, data=data.frame(x=1:ncol(tmp_tf_acc),y=x))
    #predict(tmp_lm, x=seq(1,15,0.1))
    predict(tmp_lm, newdata=data.frame(x=seq(1,ncol(tmp_tf_acc),1)))
}))

options(repr.plot.width=8, repr.plot.height=8)
#pheatmap::pheatmap(tmp_tf_acc,cluster_cols = F)
Heatmap(tmp_tf_acc,
          #col=circlize::colorRamp2(c(-1,-0.5,0,0.5,1),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          #col=circlize::colorRamp2(c(0,1,2,3,4),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          cluster_rows = T, cluster_columns = F,
          show_row_names = T, show_column_names = T,
          use_raster=FALSE,
          border = TRUE,
          #row_split=factor(time_levels, levels=time_levels),
          #column_split = rownames(dat_fc2),
          border_gp = gpar(col = "black",lwd = 0.2) ,
          column_names_gp = gpar(fontsize=12),
         na_col="#2c7bb6"
)

meanVar_pseudo2 = meanVar_pseudo[,order(as.numeric(gsub('P','',colnames(meanVar_pseudo))))]
meanVar_pseudo2 = t(apply(na.omit(meanVar_pseudo2), 1, function(x){
    tmp_lm = loess(formula = y~x, data=data.frame(x=1:ncol(meanVar_pseudo2),y=x))
    #predict(tmp_lm, x=seq(1,15,0.1))
    predict(tmp_lm, newdata=data.frame(x=seq(1,ncol(meanVar_pseudo2),1)))
}))
meanVar_pseudo2[1:3,1:3]

tmp_tf_acc[1:3,1:3]

all_var_df = meanVar_pseudo2[sig_motif$motif,]


all_acc_df = c()
for(select_tf in sig_motif$sig_gene){
    tmp_gene = hg38_mm10_TF_mapping[hg38_mm10_TF_mapping$tf_name==select_tf,'Mouse_Symbol']
    if(length(tmp_gene)==1){
        x = tmp_tf_acc[tmp_gene, ]
    }else{
        x = colMeans(tmp_tf_acc[tmp_gene, ,drop=F])
    }
    all_acc_df = rbind(all_acc_df, x)
}
rownames(all_acc_df) = sig_motif$sig_gene
colnames(all_acc_df) = colnames(all_var_df)

all_causality = t(sapply(1:nrow(all_acc_df), function(x){
    xx = cor.test(all_acc_df[x,],all_var_df[x,])
    return(c(xx$estimate,xx$p.value))
})) %>% as.data.frame()
colnames(all_causality) = c('cor', 'pvalue')
all_causality$tf = sig_motif$sig_gene

all_causality = c()
all_acc_df = c()
all_var_df = c()
for(select_tf in unique(sig_motif$sig_gene)){
    tmp_gene = hg38_mm10_TF_mapping[hg38_mm10_TF_mapping$tf_name==select_tf,'Mouse_Symbol']
    if(length(tmp_gene)==1){
        x = tmp_tf_acc[tmp_gene, ,drop=T]
    }else{
        x = colMeans(tmp_tf_acc[tmp_gene, ,drop=F])
    }
    tmp_motif = sig_motif[sig_motif$sig_gene==select_tf, 'motif']
    if(length(tmp_motif)==1){
        y = meanVar_pseudo2[tmp_motif,,drop=T]
    }else{
        y = colMeans(meanVar_pseudo2[tmp_motif,,drop=F])
    }
    all_acc_df = rbind(all_acc_df, x)
    all_var_df = rbind(all_var_df,y)
    xx = cor.test(x,y)
    #causal_edges$cor = xx$estimate
    #causal_edges$cor_pvalue = xx$p.value
    all_causality = rbind(all_causality, c(xx$estimate,xx$p.value,select_tf))
}
all_causality = as.data.frame(all_causality)
colnames(all_causality) = c('cor', 'pvalue', 'tf')
all_causality$cor = as.numeric(all_causality$cor)
all_causality$pvalue = as.numeric(all_causality$pvalue)

rownames(all_acc_df) = all_causality$tf
rownames(all_var_df) = all_causality$tf

all_causality2 = all_causality %>% filter(pvalue<0.05, abs(cor)>0.5)
acc_ht = all_acc_df[all_causality2$tf,]
var_ht = all_var_df[all_causality2$tf,]

options(repr.plot.width=8, repr.plot.height=16)
ht = Heatmap(cbind(acc_ht,var_ht),
          #col=circlize::colorRamp2(c(-1,-0.5,0,0.5,1),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          #col=circlize::colorRamp2(c(0,1,2,3,4),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          cluster_rows = T, cluster_columns = F,
          show_row_names = T, show_column_names = T,
          use_raster=FALSE,
          border = TRUE,
          #row_split=factor(time_levels, levels=time_levels),
          #column_split = rownames(dat_fc2),
          border_gp = gpar(col = "black",lwd = 0.2) ,
          column_names_gp = gpar(fontsize=12),
         na_col="#2c7bb6"
)

ht1 = Heatmap(acc_ht[row_order(ht),],
          #col=circlize::colorRamp2(c(-1,-0.5,0,0.5,1),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          #col=circlize::colorRamp2(c(0,1,2,3,4),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          cluster_rows = F, cluster_columns = F,
          show_row_names = T, show_column_names = F,
          use_raster=FALSE,
          border = TRUE,
          #row_split=factor(time_levels, levels=time_levels),
          #column_split = rownames(dat_fc2),
          border_gp = gpar(col = "black",lwd = 0.4) ,
              heatmap_legend = list(title='TF_acc'),
          row_names_gp = gpar(fontsize=10),
         na_col="#2c7bb6"
)

ht2 = Heatmap(var_ht[row_order(ht),],
              #col = paletteContinuous(set = "solarExtra", n = 100),
          col=circlize::colorRamp2(c(-2,-1,0,1,2),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          #col=circlize::colorRamp2(c(0,1,2,3,4),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          cluster_rows = F, cluster_columns = F,
          show_row_names = T, show_column_names = F,
          use_raster=FALSE,
          border = TRUE,
          #row_split=factor(time_levels, levels=time_levels),
          #column_split = rownames(dat_fc2),
          border_gp = gpar(col = "black",lwd = 0.4) ,
              heatmap_legend = list(title='Motif Var'),
          row_names_gp = gpar(fontsize=10),
         na_col="#2c7bb6"
)

ht3 = Heatmap(all_causality2[row_order(ht),'cor', drop=F],
              col = paletteContinuous(set = "solarExtra", n = 100),
          #col=circlize::colorRamp2(c(-2,-1,0,1,2),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          #col=circlize::colorRamp2(c(0,1,2,3,4),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          cluster_rows = F, cluster_columns = F,
          show_row_names = T, show_column_names = F,
          use_raster=FALSE,
              width=unit(5,'mm'),
          #border = TRUE,
          #row_split=factor(time_levels, levels=time_levels),
          #column_split = rownames(dat_fc2),
          #border_gp = gpar(col = "black",lwd = 0.4) ,
          row_names_gp = gpar(fontsize=10),
              heatmap_legend = list(title='Correlation'),
         na_col="#2c7bb6"
)

options(repr.plot.width=12, repr.plot.height=12)
ht3+ht1+ht2

options(repr.plot.width=6, repr.plot.height=8)
set.seed(1234)
ht=Heatmap(tf_ht_order_data,
              #col = paletteContinuous(set = "solarExtra", n = 100),
          col=circlize::colorRamp2(c(-2,-1,0,1,2),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          #col=circlize::colorRamp2(c(0,1,2,3,4),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          cluster_rows =T, cluster_columns = F,
          show_row_names = T, show_column_names = F,
          use_raster=FALSE,
          border = TRUE,
          #row_split=factor(time_levels, levels=time_levels),
          #column_split = rownames(dat_fc2),
          border_gp = gpar(col = "black",lwd = 0.4) ,
              heatmap_legend = list(title='Motif Var'),
          row_names_gp = gpar(fontsize=10),
         na_col="#2c7bb6"
)
ht
ht_cl = cutree(as.hclust(row_dend(ht)),15)
ht_cl = sort(ht_cl)

htx = Heatmap(tf_ht_order_data[names(ht_cl),],
              #col = paletteContinuous(set = "solarExtra", n = 100),
          col=circlize::colorRamp2(c(-2,-1,0,1,2),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          #col=circlize::colorRamp2(c(0,1,2,3,4),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          cluster_rows = F, cluster_columns = F,
          show_row_names = F, show_column_names = F,
          use_raster=FALSE,
          border = TRUE,
          #row_split=factor(ht_cl, levels = c(8,2,6,3,1,4,5,7,9,10)),
          row_split=factor(ht_cl, levels = c(2,5,3,4,6,7,15,1,8,11,9,10,12,13,14)),
          #row_split=factor(ht_cl, levels = c(1,6,7,2,5,8,9,10,4,3)),
          #column_split = rownames(dat_fc2),
          border_gp = gpar(col = "black",lwd = 0.4) ,
              heatmap_legend = list(title='Motif Var'),
          row_names_gp = gpar(fontsize=10),
         na_col="#2c7bb6"
)
htx

saveRDS(list(var_ht=tf_ht_order_data,
             ht_cl = ht_cl,
             all_causality2=all_causality2), './placeholder_analysis/round_cluster02/merge/AT2_TF_var_ht.rds')

AT2_TF_var_ht = readRDS('./placeholder_analysis/round_cluster02/merge/AT2_TF_var_ht.rds')
var_ht = AT2_TF_var_ht$var_ht
ht_cl = AT2_TF_var_ht$ht_cl
all_causality2 = AT2_TF_var_ht$all_causality2

htx = Heatmap(tf_ht_order_data[names(ht_cl),],
              #col = paletteContinuous(set = "solarExtra", n = 100),
          col=circlize::colorRamp2(c(-2,-1,0,1,2),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          #col=circlize::colorRamp2(c(0,1,2,3,4),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          cluster_rows = F, cluster_columns = F,
          show_row_names = F, show_column_names = F,
          use_raster=FALSE,
          border = TRUE,
          #row_split=factor(ht_cl, levels = c(8,2,6,3,1,4,5,7,9,10)),
          row_split=factor(ht_cl, levels = c(2,5,3,4,6,7,15,8,1,11,9,10,12,13,14)),
          #row_split=factor(ht_cl, levels = c(1,6,7,2,5,8,9,10,4,3)),
          #column_split = rownames(dat_fc2),
          border_gp = gpar(col = "black",lwd = 0.4) ,
              heatmap_legend = list(title='Motif Var'),
          row_names_gp = gpar(fontsize=10),
         na_col="#2c7bb6"
)
htx

tmp_new_order = names(ht_cl)[unlist(row_order(htx))]

AT2_TF_var_ht = readRDS('./placeholder_analysis/round_cluster02/merge/AT2_TF_var_ht.rds')
var_ht = AT2_TF_var_ht$var_ht
ht_cl = AT2_TF_var_ht$ht_cl

#write.table(motif_meta, '../share_motifmeta.txt', quote=F, sep='\t', row.names = F)
motif_meta = read.table('../share_motifmeta.txt', sep='\t', header = T)

rownames(motif_meta) = motif_meta$motif

names(ht_cl)

library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)

var_order = unlist(row_order(htx))

sig_motif = enrich.mat %>%
    filter(p.adjust<0.05) %>%
    arrange(-fold.enrichment)

motif_region = StringToGRanges(unique(top_peaks$gene))

my_CollapseToLongestTranscript <- function (ranges) 
{
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(test = range.df$strand == "*", 
    yes = "+", no = range.df$strand)
  collapsed <- range.df[, .(unique(seqnames), min(start), 
    max(end), strand[[1]], gene_biotype[[1]], gene_name[[1]]), 
    "gene_id"]
  colnames(x = collapsed) <- c("gene_id", "seqnames", "start", 
    "end", "strand", "gene_biotype", "gene_name")
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(df = collapsed, 
    keep.extra.columns = TRUE)
  return(gene.ranges)
}


gene.coords <- my_CollapseToLongestTranscript(ranges = Annotation(srt))

tss = promoters(gene.coords, upstream = 2000, downstream = 0)

#tss = resize(gene.coords, width=1, fix='start')

hits = findOverlaps(motif_region, tss)

CRE_TssGene = data.frame(
    'CRE' = GRangesToString(motif_region[queryHits(hits)]),
    'Gene' = tss[subjectHits(hits)]$gene_name
)
#rownames(CRE_TssGene) = CRE_TssGene$CRE

CRE_TssGene = CRE_TssGene[!grepl('Rik', CRE_TssGene$Gene),]
CRE_TssGene = CRE_TssGene[!grepl('Gm', CRE_TssGene$Gene),]
CRE_TssGene = CRE_TssGene[!grepl('-', CRE_TssGene$Gene),]
CRE_TssGene = CRE_TssGene[!duplicated(CRE_TssGene$CRE),]

rownames(CRE_TssGene) = CRE_TssGene$CRE

sig_motifs = unique(sig_motif$motif)
motif_data = as.matrix(GetMotifData(srt))[CRE_TssGene$CRE,sig_motifs]

GRN_init = melt(motif_data) %>% 
    filter(value==TRUE)
GRN_init = GRN_init[,1:2]
colnames(GRN_init) = c('CRE', 'Motif')
GRN_init$CRE = as.vector(GRN_init$CRE)
GRN_init$Motif = as.vector(GRN_init$Motif)

GRN_init$closest_gene = CRE_TssGene[GRN_init$CRE, 'Gene']
GRN_init[,c('family', 'TF')] = motif_meta[GRN_init$Motif, c('family', 'TF')]
GRN_init = na.omit(GRN_init)
GRN_init = GRN_init %>% filter(TF%in%all_causality2$tf)
head(GRN_init)
dim(GRN_init)
cre_acces_pseudo = myAveragePeak(srt, 'Time', assay='CRE', slot='data')
cre_acces_pseudo = log2(cre_acces_pseudo+1)
cre_acces_pseudo = as.matrix(cre_acces_pseudo)
cre_acces_pseudo = myRowScale(cre_acces_pseudo, max = 2, min = -2,limit=TRUE)

meanVar_pseudo = myAveragePeak(srt, 'Time', assay='chromvar', slot='data')
meanVar_pseudo = as.matrix(meanVar_pseudo)
meanVar_pseudo = myRowScale(meanVar_pseudo, max = 2, min = -2,limit=TRUE)

pesutotime_levels = colnames(cre_acces_pseudo)[order(as.numeric(gsub('P','',colnames(cre_acces_pseudo))))]
GRN_init_cor = apply(GRN_init, 1, function(x){
    a = cre_acces_pseudo[x['CRE'], pesutotime_levels]
    b = meanVar_pseudo[x['Motif'], pesutotime_levels]
    xx=cor.test(a,b)
    return(c(xx$estimate, xx$p.value))
})
GRN_init_cor = as.data.frame(t(GRN_init_cor))
head(GRN_init_cor)
# keepmotifagainCREcenter50bpwithin range
cre_center = resize(granges(srt), width=100, fix='center')

motif_pwm = toPWM(pfm_filter)

library(motifmatchr)
motif_position = matchMotifs(
    pwms=motif_pwm,
    subject = cre_center,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    out = 'position'
)

motif_position_filter = GRangesList()
for(tmp_tf in names(motif_position)){
    i = motif_position[[tmp_tf]]
    i2 = i[i$score>quantile(i$score, 0.75)]
    motif_position_filter[[tmp_tf]] = i2
}

motif_pos_gr = unlist(motif_position_filter)

#motif_pos_gr = unlist(motif_position)

peak_motif_hits = findOverlaps(cre_center, motif_pos_gr)

cre_motif_df = data.frame(
    'CRE'=GRangesToString(granges(srt))[queryHits(peak_motif_hits)],
    'motif'=names(motif_pos_gr)[subjectHits(peak_motif_hits)]
)

cre_motif_df$combine_name = paste0(cre_motif_df$CRE, cre_motif_df$motif)


GRN_init_with_cor = cbind(GRN_init, GRN_init_cor)

GRN_init_with_cor = GRN_init_with_cor %>%
    mutate(adj = p.adjust(V2)) %>%
    filter(V2<0.01, cor>0) 

GRN_init_with_cor$combine_name = paste0(GRN_init_with_cor$CRE, GRN_init_with_cor$Motif)

GRN_init_with_cor = GRN_init_with_cor[GRN_init_with_cor$combine_name%in%cre_motif_df$combine_name,]

library(igraph)

GRN_init_with_cor2 = GRN_init_with_cor[,c('TF','CRE')]
GRN_init_with_cor2 = GRN_init_with_cor2 %>% unique()
g <- graph_from_data_frame(GRN_init_with_cor2[,c('TF','CRE')] , directed = TRUE)

write.table(GRN_init_with_cor2, glue('{output_dir}/AT2_TF_network_cytoscape.txt'), quote = F, sep='\t', row.names = F)

unique_TF = unique(GRN_init_with_cor2$TF)
unique_TF

options(repr.plot.width=8, repr.plot.height=12)
# unique_TF = c(unique(need_label_gene$TF), 'Nkx2-1','EGR3', 'SP9', 'SP8')
label_pos = na.omit(match(unique_TF,tmp_new_order))
ra = rowAnnotation(TF = anno_mark(at=label_pos, labels=unique_TF,padding=unit(8,'mm')))
#top_an = HeatmapAnnotation(Time=)
ht1 = Heatmap(var_ht[tmp_new_order,],
              #col = paletteContinuous(set = "solarExtra", n = 100),
          col=circlize::colorRamp2(c(-2,-1,0,1,2),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          #col=circlize::colorRamp2(c(0,1,2,3,4),c("#2c7bb6",'#abd9e9',"white","#fdae61","#d7191c")),
          cluster_rows = F, cluster_columns = F,
          show_row_names = F, show_column_names = F,
        right_annotation=ra,
          use_raster=T,
          border = TRUE,
                    #row_split=factor(ht_cl, levels = c(1,7,3,8,2,5,9,10,6,4)),
          #column_split = rownames(dat_fc2),
          border_gp = gpar(col = "black",lwd = 0.4) ,
              heatmap_legend = list(title='Motif Var'),
          row_names_gp = gpar(fontsize=10),
         na_col="#2c7bb6"
)
ht1

pdf(glue('{output_dir}/motifchange_heatmap_mark_transcription_factors.pdf'), width = 8,height = 12)
draw(ht1)
dev.off()

all_nodes = names(V(g))

top_n = 1:4
col_name_levels2 = col_name_levels[top_n]
nodes_color_value = sapply(all_nodes, function(x){
    if(x%in%GRN_init_with_cor$CRE){
        tmp_value = rowMeans(cre_acces_pseudo[x, col_name_levels2, drop=F])
        tmp_value=0
    }else{
        
        x2 = motif_meta[motif_meta$TF==x, 'motif']
        tmp_value = mean(rowMeans(meanVar_pseudo[x2,col_name_levels2, drop=F]))
    }
    return(tmp_value)
})
nodes_color_value_stage1 = unname(unlist(nodes_color_value))

top_n = 5:12
col_name_levels2 = col_name_levels[top_n]
nodes_color_value = sapply(all_nodes, function(x){
    if(x%in%GRN_init_with_cor$CRE){
        tmp_value = rowMeans(cre_acces_pseudo[x, col_name_levels2, drop=F])
        tmp_value=0
    }else{
        
        x2 = motif_meta[motif_meta$TF==x, 'motif']
        tmp_value = mean(rowMeans(meanVar_pseudo[x2,col_name_levels2, drop=F]))
    }
    return(tmp_value)
})
nodes_color_value_stage2 = unname(unlist(nodes_color_value))

top_n = 13:15
col_name_levels2 = col_name_levels[top_n]
nodes_color_value = sapply(all_nodes, function(x){
    if(x%in%GRN_init_with_cor$CRE){
        tmp_value = rowMeans(cre_acces_pseudo[x, col_name_levels2, drop=F])
        tmp_value=0
    }else{
        
        x2 = motif_meta[motif_meta$TF==x, 'motif']
        tmp_value = mean(rowMeans(meanVar_pseudo[x2,col_name_levels2, drop=F]))
    }
    return(tmp_value)
})
nodes_color_value_stage3 = unname(unlist(nodes_color_value))

all_node_features = data.frame('node_name'=names(V(g)),
            'node_type' = sapply(names(V(g)), function(x){ifelse(x %in% GRN_init_with_cor$CRE, 'CRE', 'TF')}),
            'CRE_value_stage1' = nodes_color_value_stage1,
            'CRE_value_stage2' = nodes_color_value_stage2,
            'CRE_value_stage3' = nodes_color_value_stage3
            
            )

node_features = data.frame(
    'degree' = degree(g),
    'degree_out' = degree(g, mode='out'),
    'degree_in' = degree(g, mode='in'),
    'betweenness' = betweenness(g, directed = FALSE),
    #'closeness' = closeness(g, mode='all', normalized = TRUE),
    'eigen_central' = eigen_centrality(g, directed = FALSE)$vector,
    'pagerank' = page_rank(g, directed = TRUE)$vector
)
node_features$node_type = 'TF'
node_features[rownames(node_features)%in%GRN_init_with_cor$CRE, 'node_type'] = 'CRE'
node_features$node_name=rownames(node_features)

node_features2 = melt(node_features, id=c('node_type','node_name'))

all_node_features[rownames(node_features), c('degree', 'degree_out', 'degree_in')] = node_features[,1:3]

hmm_ann = read.csv('./placeholder_analysis/round_cluster02/cCRE//all_merge_peaks_HMMann.csv')

hmm_ann$cre = gsub(':','-',hmm_ann$cre)

hmm_ann_filter = hmm_ann %>% filter(cre%in%node_features2$node_name)

hmm_ann_filter = hmm_ann_filter %>% filter(Tissue=='lung')

hmm_ann_filter_mat = reshape2::acast(hmm_ann_filter, cre~hmm_state)

hmm_states_cre = colnames(hmm_ann_filter_mat)[apply(hmm_ann_filter_mat, 1, which.max)]
names(hmm_states_cre) = rownames(hmm_ann_filter_mat)

all_node_features$hmm_state = hmm_states_cre[all_node_features$node_name]

motif_meta2 = motif_meta[, c('TF', 'family_abbr')] %>% unique()
rownames(motif_meta2) = motif_meta2$TF
head(motif_meta2 )

all_node_features[, 'family'] = motif_meta2[all_node_features$node_name, 'family_abbr']
all_node_features

write.table(all_node_features, glue('{output_dir}/AT2_TF_network_cytoscape_meta.txt'), quote = F, sep='\t', row.names = F)

node_features = data.frame(
    'degree' = degree(g),
    'degree_out' = degree(g, mode='out'),
    'degree_in' = degree(g, mode='in'),
    'betweenness' = betweenness(g, directed = FALSE),
    #'closeness' = closeness(g, mode='all', normalized = TRUE),
    'eigen_central' = eigen_centrality(g, directed = FALSE)$vector,
    'pagerank' = page_rank(g, directed = TRUE)$vector
)
node_features$node_type = 'TF'
node_features[rownames(node_features)%in%GRN_init_with_cor$CRE, 'node_type'] = 'CRE'
node_features$node_name=rownames(node_features)

node_features2 = melt(node_features, id=c('node_type','node_name'))

# hmm_ann = read.csv('./placeholder_analysis/round_cluster02/cCRE//all_merge_peaks_HMMann.csv')

# hmm_ann$cre = gsub(':','-',hmm_ann$cre)

# hmm_ann_filter = hmm_ann %>% filter(cre%in%node_features2$node_name)

# hmm_ann_filter = hmm_ann_filter %>% filter(Tissue=='lung')

# hmm_ann_filter_mat = reshape2::acast(hmm_ann_filter, cre~hmm_state)

# hmm_states_cre = colnames(hmm_ann_filter_mat)[apply(hmm_ann_filter_mat, 1, which.max)]
# names(hmm_states_cre) = rownames(hmm_ann_filter_mat)

node_features2$hmm_state = hmm_states_cre[node_features2$node_name]
node_features2$hmm_state[is.na(node_features2$hmm_state)] = 'other'

head(node_features2)

chromHMM_colors_simple <- c(
  "Pr-A" = "#006400",
  "Pr-W" = "#90EE90",
  "Pr-B" = "#D3D3D3",
  "Pr-F" = "#008000",
  "En-Sd" = "#FFFF00",
  "En-Sp" = "#FFFF00",
  "En-W" = "#FFFFE0",
  "En-Pd" = "#696969",
  "En-Pp" = "#696969",
  "Tr-S" = "#00008B",
  "Tr-P" = "#4169E1",
  "Tr-I" = "#ADD8E6",
  "Hc-P" = "#FA8072",
  "Hc-H" = "#FFC0CB",
  "NS"   = "#FFFFFF"
)
library(forcats)
options(repr.plot.width=6, repr.plot.height=8)
a = node_features2 %>% 
    filter(variable=='degree_in') %>%
    #filter(hmm_state!='NS') %>%
    filter(node_type!='TF') %>%
    filter(value>=2) %>%
    mutate(label_name=ifelse(value>=5, node_name, '')) %>%
    mutate(node_name=fct_reorder(node_name,value)) %>%
    ggplot(aes(y=node_name, x=value, fill=hmm_state))+
    geom_bar(stat='identity', color='black', size=0.1)+
    #geom_segment(aes(x=0, xend=value, yend=node_name), size=0.1)+
    #geom_point(size=2, shape=21)+
    geom_text(aes(x=4,label=node_name), color='black',force=4, force_pull = 0.1, size=label_size(6),
                   #nudge_x=2, nudge_y=-1
             )+
    scale_fill_manual(values = chromHMM_colors_simple)+
    theme_bw()+
    mytheme+
    theme(axis.text.y = element_blank(),
          axis.ticks.y=element_blank())
a

cre_order = rev(levels(a$data$node_name))

ggsave(glue('{output_dir}/AT2_networkCREin-degree(>1).pdf'), a,
       width=70, height=90, units='mm', dpi=600, bg='transparent')

write.table(node_features2, '../pycode/round_cluster02/cCRE/AT2_net_degree.txt')

node_features2$group='AT2'

node_features_AT1 = read.table('../pycode/round_cluster02/cCRE/AT1_net_degree.txt')
node_features_AT1$group='AT1'

node_features_all = rbind(node_features2, node_features_AT1)

node_features_all$node_name_new = paste0(node_features_all$node_name, ':', node_features_all$group)
node_features_all = node_features_all %>% 
    filter(variable=='degree')

node_features_all = node_features_all %>% arrange(node_type, value)
node_features_all$node_name_new = factor(node_features_all$node_name_new,node_features_all$node_name_new)

tmp_meta = motif_meta[motif_meta$TF%in%node_features_all$node_name,c('TF', 'family_abbr')] %>% unique()
rownames(tmp_meta) = tmp_meta$TF

node_features_all$TF_family = tmp_meta[node_features_all$node_name, 'family_abbr']

node_features_all_01=reshape2::dcast(node_features_all,
                                     node_name+hmm_state+TF_family~group, 
                                     value.var = 'value',
                                     fill=0)
head(node_features_all_01)

#node_features_all_01$AT1 = node_features_all_01$AT1 / 

node_features_all_01$group = 'TF'
node_features_all_01$group[node_features_all_01$hmm_state!='other'] = 'CRE'

a=node_features_all %>% 
    #filter(hmm_state!='NS') %>%
    #filter(node_type!='TF') %>%
    #filter(value>=2) %>%
    mutate(label_name=ifelse(value>=5, node_name, '')) %>%
    mutate(node_name=fct_reorder(node_name,value)) %>%
    ggplot(aes(x=node_name_new, y=value, fill=group))+
    geom_bar(stat='identity', color=NA)+
    #geom_segment(aes(x=0, xend=value, yend=node_name), size=0.1)+
    #geom_point(size=2, shape=21)+
    #geom_text(aes(x=4,label=node_name), color='black',force=4, force_pull = 0.1, size=label_size(6),
                   #nudge_x=2, nudge_y=-1
    #         )+
    #scale_fill_manual(values = chromHMM_colors_simple)+
    facet_wrap(~node_type, ncol=2, scales = 'free')+
    theme_bw()+
    mytheme+
    theme(axis.text.x = element_blank(),
          axis.ticks.x=element_blank())
a

ggsave(glue('{output_dir}/AT1_AT2_network_degree_comparison.pdf'), a,
       width=180, height=30, units='mm', dpi=600, bg='transparent')

a=node_features_all %>% 
    filter(node_type!='TF') %>%
    #filter(group=='AT1') %>%
    filter(value>=4) %>%
    #top_n(10,value) %>%
    mutate(label_name=ifelse(value>=5, node_name, '')) %>%
    mutate(node_name=fct_reorder(node_name,value)) %>%
    ggplot(aes(y=node_name, x=value, fill=hmm_state))+
    geom_bar(stat='identity', color='black', size=0.1)+
    #geom_segment(aes(x=0, xend=value, yend=node_name), size=0.1)+
    #geom_point(size=2, shape=21)+
    geom_text(aes(x=4,label=node_name), color='black',force=4, force_pull = 0.1, size=label_size(6),
                   #nudge_x=2, nudge_y=-1
             )+
    scale_fill_manual(values = chromHMM_colors_simple)+
    theme_bw()+
    mytheme+
    theme(axis.text.y = element_blank(),
          axis.ticks.y=element_blank())
a
ggsave(glue('{output_dir}/AT1_CRE_network_degree_comparison.pdf'), a,
       width=35, height=40, units='mm', dpi=600, bg='transparent')

node_features_all

node_features_all%>% 
    filter(node_type=='TF') %>%
    filter(node_name=='TEAD4')

tmp_node_features_all = node_features_all%>% 
    filter(node_type=='TF', group=='AT2') %>%
    
    #filter(value>=20) %>%
    top_n(7,value)

tmp_node_features_all$stage = ifelse(tmp_node_features_all$node_name%in%c('SP4','HINFP'),  'Stage2', 'Stage1')
tmp_node_features_all

a=tmp_node_features_all %>%
    mutate(label_name=ifelse(value>=5, node_name, '')) %>%
    mutate(node_name=fct_reorder(node_name,value)) %>%
    ggplot(aes(y=node_name, x=value, fill=stage))+
    geom_bar(stat='identity', color='black', size=0.1)+
    #geom_segment(aes(x=0, xend=value, yend=node_name), size=0.1)+
    #geom_point(size=2, shape=21)+
    geom_text(aes(x=4,label=node_name), color='black',force=4, force_pull = 0.1, size=label_size(6),
                   #nudge_x=2, nudge_y=-1
             )+
    #scale_fill_manual(values = chromHMM_colors_simple)+
    scale_fill_npg()+
    theme_bw()+
    mytheme+
    theme(axis.text.y = element_blank(),
          axis.ticks.y=element_blank())
a
ggsave(glue('{output_dir}/AT2_TF_network_degree_comparison.pdf'), a,
       width=35, height=40, units='mm', dpi=600, bg='transparent')

