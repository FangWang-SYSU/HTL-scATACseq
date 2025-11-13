
### 1.单个样本的质量控制过程
#' @description  
#' SingleSample.QC.scATAC
#' 
# 调用 R markdown
#' @param datapath 提供数据集的名称和地址
#' @param outPath 数据输出的地址
#' @param sampleName description
#' @param min.cells default min.cells=10
#' @param min.features default min.features = 100
#' @param Uniquefragment default Uniquefragment =1000
#' @param TSS.Enrichment.min default TSS.Enrichment.min = 3
#' @param NS.sign default NS.sign = 2,
#' @param fragmentsInPeaks peaks区域的fragment数量
#' @param pct_reads_in_peaks peak区域的fragment的比例
#' @param QC.features 展示的指标
#' @param dims default  dims = 2:50
#' @param resolution 聚类时的分辨率，default resolution = 0.8
#' @param projectName  default projectName = NULL
#' @param default F; 是否手动的注视细胞类型

SingleSampleQCscATACHtml <- function(dataPath = NULL,outPath = NULL,  Species = 'hg38',
                                   sampleName = NULL,
                                   min.cells = 10,
                                   min.features = 100,
                                   Uniquefragment = 1000,
                                   TSS.Enrichment.min = 3,
                                   NS.sign = 2,
                                   fragmentsInPeaks = NULL,
                                   pct_reads_in_peaks = 30,
                                   QC.features = NULL,
                                   dims = 2:50,
                                   resolution = 0.8,
                                   peak.gene = c('Ptprc','Cd3d','Cdh5','Epcam'),
                                   projectName = NULL,
                                   CellInfoPath = NULL,
                                   CodePath = NULL,
                                   output_file = NULL,
                                   output_dir = NULL){
 
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RMDFile/1.SingleSampleQualityControl.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
}




#' @description
#' A short description...
#' 
# 调用 R markdown
#' @param obj seurat 对象
#' @param sampleName 项目名称
#' @param outPath 数据存储地址
#' @param filterDoublet 是否使用doublet方法过滤双细胞，默认过滤
#' @param TSS.enrichment.min 过滤指标1
#' @param passed_filters.min 过滤指标2
#' @param nFeature_ATAC.min 过滤指标1
#' @param nFeature_ATAC.max 过滤指标1
#' @param nCount_ATAC.min 过滤指标1
#' @param nCount_ATAC.max 过滤指标1
#' @param peak_region_fragments 过滤指标1
#' @param pct_reads_in_peaks 过滤指标
#' @param output_file 
#' @param output_dir  
#' 
scATACQCfilterHtml <- function(obj,
                               sampleName=NULL,
                               Species='hg38',
                               dims = 2:30,
                               outPath = NULL,
                               filterDoublet = T,
                               fragPath = NULL, 
                               TSS.enrichment.min = NULL, 
                               passed_filters.min = NULL,
                               nFeature_ATAC.min = NULL, 
                               nFeature_ATAC.max = NULL,
                               nCount_ATAC.min = NULL, 
                               nCount_ATAC.max = NULL,
                               peak_region_fragments = NULL,
                               pct_reads_in_peaks = NULL,
                               CellInfoPath = NULL,
                               CodePath = NULL,
                               output_file=NULL,
                               output_dir = NULL){
    
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"/RmdFile/2.SingleSample_QCFilter.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
}


 
#' @description
#' OptinalHarmony
#' @param seuratObj description
#' @param RunUmap 是否使用UMAP进行降维，当输入数据没有进行UMAP降维时，必须设置为Ture，默认为FALSE
#' @param groupBy description
#' @param dims description
#' @param assay description
#' @param output_file description
#' @param output_dir description
OptinalHarmonyHtml <- function(seuratObj,
                               RunUmap = F,
                               groupBy = NULL,
                               dims = 2:30,
                               assay = 'ATAC',
                               CodePath = NULL,
                               outPath = NULL,
                               output_file = NULL,
                               output_dir = NULL){
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"/RmdFile/Optional_Harmony.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
}




#' @description  ClusteringAndTyping
#' clustering过程1
# 调用 R markdown
#' @param obj 提供数据集 
#' @param outPath 数据输出的路径
#' @param focus.features 关注的特征，样本的分组信息，用于输出UMAP图展示不同分组信息的分布
#' @param assay 用于聚类的assay
#' @param Clustering 是否进行聚类分析，默认为不进行，需提前对数据进行聚类分析，如果设置为T，则需提供dims等信息
#' @param dims 聚类的维度，默认为2:50
#' @param marker.peaks
#' @param max.cells.per.ident 计算差异peaks时，每组随机选择多少个细胞，为了节约时间，默认为500; 正常分析时可以设置 max.cells.per.ident = Inf
#' @param topN 想要展示的top 个数
#' @param cluster.heatmap.topN 绘制cluster水平的热图时，每一个cluster展示的特征
#' @param nLabel  绘制cluster水平的热图时，每一个cluster展示名字的个数
#' @param PrintMarker 是否展示特定的marker gene，输入基因，将展示在RNA cluster-level热图中展示，如果不设置，将输出每个cluster topN的基因名
#' @param output_file html文件的输出名
#' @param output_dir html文件的输出地址


scATACClusteringHtml <- function(object = NULL,
                                outPath = NULL,
                                focus.features = NULL,group_by = NULL,
                                assay = NULL,
                                Clustering = F,
                                dims = 2:50,
                                marker.peaks = NULL,
                                max.cells.per.ident = 500,
                                topN = 5,
                                cluster.heatmap.topN = NULL,
                                nLabel = 3,
                                PrintMarker = NULL,
                                CodePath = NULL,
                                output_file=NULL,
                                output_dir = NULL){
    
    # 调用rmd文档
    rmarkdown::render( input = file.path(CodePath,"RmdFile/3.ClusterAndTyping.Rmd"), 
                       output_dir = output_dir,
                       output_file = output_file)
}

#' @description
#' A short description...
#' @param obj description
#' @param ncol description
#' @param markers 展示的基因组区域，向量，可以为chr1-1243324-34324123，也可以为指定基因名
#' @param plot.expr 是否在侧边增加表达，默认不增加，如果设置如TURE则增加表达，设置为T时，输入的必须是基因名

CovergePlotHtml <- function(obj, ncol = 4, groupBy = NULL,
                        markers = NULL,
                        plot.expr = F,
                        assay = NULL,
                        CodePath = NULL,    
                        output_dir,
                        output_file
){
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RmdFile/Optional_CoveragePlot.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
}


 

 
#' @description  CellType.Typing.SeurObj
# 调用 R markdown
#' @param obj 提供数据集 
#' @param outPath 数据输出的路径
#' @param focus.group 关注的特征，样本的分组信息，用于输出UMAP图展示不同分组信息的分布
#' @param pre.cell 预测的细胞类型结果 data.frame 行名为cluster，列名分别为cluster和CellType
#' @param atac_assay_name  
#' @param cell.heatmap.topN 绘制cell水平的热图时，每一个cluster展示的特征
#' @param cluster.heatmap.topN  绘制cluster水平的热图时，每一个cluster展示名字的个数
#' @param nLabel  如果不设置，将输出每个cluster topN的基因名
#' @param output_file html文件的输出名
#' @param output_dir html文件的输出地址

CellType.Typing.SeurObj <- function(obj = NULL,
                                    outPath = NULL,
                                    focus.group = NULL,
                                    plot.Time = F,
                                    markerGene = NULL,
                                    atac_assay_name = 'ATAC',
                                    cell.heatmap.topN = 10,
                                    cluster.heatmap.topN = NULL,
                                    nLabel = 5, 
                                    CodePath = NULL,
                                    output_file=NULL,
                                    output_dir = NULL){
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RmdFile/4.Typing.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
}


#' @description 
#' @param object
#' @param markerPeak
#' @param Enrichment.Ref.peaks.celltype
#' @param outPath
#' @param CodePath
#' @param output_file
CellType.Typing.SeurObj <- function(object = NULL,
                                    markerPeak = NULL,
                                    Enrichment.Ref.peaks.celltype = F,
                                    outPath = NULL,
                                    CodePath = NULL,
                                    output_file=NULL,
                                    output_dir = NULL){
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RmdFile/4.CellTypePrediction.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
}





#' @description 
#' @param object
#' @param outPath
#' @param Species 物种信息，默认为人类，hg38
#' @param PeaksGroup.list
#' @param 
cCREs.Calculated.SeurObj <- function(object = NULL,
                                    outPath = NULL,
                                    Species = 'hg38',
                                    PeaksGroup.list =NULL,
                                    Ref.cCREs = NULL,
                                    mainType = 'CellType',
                                    other.focus.group = NULL,
                                    atac_assay_name = 'peaks',
                                    marker.peaks = NULL,
                                    CodePath = NULL,
                                    output_file=NULL,
                                    output_dir = NULL){
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RmdFile/5.cCREs_calculated.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
}



#' @description 
#' @param object
#' @param outPath
#' @param PeaksGroup.list
cCREs.Enrichment.SeurObj <- function(cCREs.dat = NULL,Species ='hg38',
                                    outPath = NULL,
                                    CodePath = NULL,
                                    output_file=NULL,
                                    output_dir = NULL){
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RmdFile/5.cCREs_Enrichment.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
}



#' @description
#' 共可及性分析
#' @param object description
#' @param assay description
#' @param ClusterBy
#' @param expression.assay RNA数据使用的assay名称
#' @param PP.link description
#' @param PG.link description
#' @param region 基因组区域（"chr1-207880000-207920000"）/某一特定基因(CD34);
#' @param CodePath
#' @param outPath
#' @param output_file
#' @param output_dir
cCREs.CoAccessible.Html <- function(object,
                                    assay = 'peaks',
                                    ClusterBy = 'CellType',
                                    expression.assay = "RNA",
                                    PP.link = peakToPeak.link,
                                    PG.link = PeakToGene.link,
                                    region = c('CD34','IL7R'),
                                    CodePath = NULL,
                                    outPath = NULL,
                                    output_file=NULL,
                                    output_dir = NULL){
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RmdFile/6.Co_accessibility.Rmd"),
                      output_dir = output_dir,
                      output_file = output_file)
}







#' @description 
#' @param object
#' @param assay ATAC peak数据使用的assay名称
#' @param ClusterBy
#' @param diff.Peaks 默认为null，计算的marker.peaks 存储在输出文件夹的DiffAccessPeak文件夹中，以diff.peak.csv命名
#' @param SpeCompare

MotifAnalysisHtml <- function(object,
                              assay = 'peaks',
                              ClusterBy = 'CellType',
                              diff.Peaks = NULL,
                              outPath = NULL,
                              CodePath = NULL,
                              output_file=NULL,
                              output_dir = NULL){
    
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RmdFile/7.Motif_Analysis.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
    
}

#' @description 
#' @param object
#' @param assay ATAC peak数据使用的assay名称
#' @param ClusterBy
#' @param diff.Peaks 默认为null，计算的marker.peaks 存储在输出文件夹的DiffAccessPeak文件夹中，以diff.peak.csv命名
#' @param SpeCompare

MotifGroupCompare.Html <- function(object,
                              assay = 'peaks',
                              ClusterBy = 'CellType',
                              diff.Peaks = NULL,
                              outPath = NULL,
                              CodePath = NULL,
                              output_file=NULL,
                              output_dir = NULL){
    
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RmdFile/7.Motif_GroupCompare.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
    
}

## 调用Html文件，输出结果
#' @description 
#' @param projHeme
#' @param assay ATAC peak数据使用的assay名称
#' @param ClusterBy
#' @param diff.Peaks 默认为null，计算的marker.peaks 存储在输出文件夹的DiffAccessPeak文件夹中，以diff.peak.csv命名
#' @param SpeCompare

TrajectoryAnalysisHtml <- function(projHeme,
                              trajectory = NULL,
                              trajectory_name = 'trajectory',
                              groupBy = NULL,
                              feature = c("Dntt","Socs2"),
                              colorBy = "GeneScoreMatrix", ## 可以选择任何ArchR中具备的Matrix，如GeneScoreMatrix, MotifMatrix
                              outPath = NULL,
                              CodePath = NULL,
                              output_file=NULL,
                              output_dir = NULL){
    
    # 调用rmd文档
    rmarkdown::render(input = file.path(CodePath,"RmdFile/8_trajectoriesArchR.Rmd"), 
                      output_dir = output_dir,
                      output_file = output_file)
    
} 