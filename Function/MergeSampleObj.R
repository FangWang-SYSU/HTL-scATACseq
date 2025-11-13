
#' @description
#' A short description...
#' @param merge.method 合并的方法; merge.method = c('commonPeaks','direct')，我们推荐使用commonPeaks的方法，默认merge.method = 'commonPeaks',
#' @param objPath description
#' @param sample_Name description
#' @param projectName description
#' @param Cellranger.path description
#' @param Species 物种，默认为人，如果处理小鼠数据，需甚至Species = 'mm10'


Merge.Process <- function(objPath,
                          sample_Name,
                          merge.method = 'commonPeaks',
                          projectName = '',
                          Cellranger.path = NULL,
                          Species = 'hg38'
){
    # ## 判断传入参数是否正确
    # if(addMeta ==T){
    #     if(is.null(sampleInfo)){
    #         stop('Error: Parameter sampleInfo is null. if you want to add metadata, you shoule input a dataframe or a matrix which include sample infomation. if not, you shoule set addMeta =F ')
    #     }
    #     if(is.matrix(sampleInfo)){
    #         sampleInfo = as.data.frame(sampleInfo)
    #     }
    #     
    #     if(!is.data.frame(sampleInfo)){
    #         stop('Error: Parameter sampleInfo must be a data.frame or a matirx')
    #     }
    # }
    
    ## 合并
    if(merge.method =='direct'){
        obj <- merge.Obj.direct(objPath = objPath,## object对象所在的位置
                                sample_Name = sample_Name,projectName = projectName)
    }else if(merge.method =='commonPeaks'){
        obj <- Merge.by.common.peaks(sample_Name = sample_Name,## 需要合并样本的名字，向量，包括所有需要合并样本的名字
                                     Cellranger.path = Cellranger.path, ## cellranger所在的位置
                                     objPath = objPath,## object对象所在的位置
                                     projectName = projectName ## 该项目的名称
        ) 
    }
    
    obj <- .AddGeneAnnota(obj,Species = Species)
    
    # if(addMeta == T){
    #     rownames(sampleInfo) = sampleInfo$sampleID
    #     sampleInfo = sampleInfo[obj$sampleID,]
    #     obj = AddMetaData(object = obj,metadata = sampleInfo)
    #     
    # }
    return(obj)
}


################




### Merge
#'@description Merge.by.common.peaks
#'A short description...
#' @param sample_Name 需要合并到样本
#' @param Cellranger.path  文件所在地址
#' @param FilterSampleObj description
#' @param projectName description
#' 

Merge.by.common.peaks <- function(sample_Name,
                                  Cellranger.path = NULL,
                                  objPath = NULL,
                                  projectName = NULL){
    
    library(Signac)
    library(Seurat)
    library(GenomicRanges)
    library(future)
    plan("multicore", workers = 4)
    options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
    
    # -----------------Peaks
    j = 1
    for(i in sample_Name){
        #i = sample_Name[1]
        peaks.path  = file.path(Cellranger.path,i,'outs/filtered_peak_bc_matrix/peaks.bed')
        peaks <- read.table(file = peaks.path,col.names = c("chr", "start", "end"))
        gr <- makeGRangesFromDataFrame(peaks)
        
        if(j==1){
            gr.all = gr
        }else{
            gr.all = c(gr.all,gr)
        }
        j = j+1
    }
    combined.peaks <- GenomicRanges::reduce(gr.all)
    peakwidths <- width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
    
    
    # -----------------
    objList = list()
    for(i in sample_Name){
        sam_obj =  readRDS(file.path(objPath,paste0(i,'.RDS')))
        fragPath = file.path(Cellranger.path,i,'outs/fragments.tsv.gz')
        # Metadata
        singlecells <- sam_obj@meta.data
        # Fragments
        frags <- CreateFragmentObject(path = fragPath,cells = rownames(singlecells))
        # Counts
        obj.counts <- FeatureMatrix( fragments = frags,features = combined.peaks,cells = rownames(singlecells))
        # Assay
        Frag.assay <- CreateChromatinAssay(obj.counts, fragments = frags)
        # Object
        merge.obj <- CreateSeuratObject(Frag.assay, assay = "ATAC", meta.data=singlecells)
        merge.obj$sampleID = i
        objList[[i]] = merge.obj
    }
    
    combinedObj <- merge(x = objList[[1]],y = objList[2:length(objList)] ,add.cell.ids = paste0(projectName,'_',sample_Name))
    return(combinedObj)
}

#' @description merge.Obj.direct
#' 合并obj对象
#' @param data.list list 对象，每一个list为一个obj对象
#' @param projectName 项目名称，该名称会添加在合并后的barcode前面
#' @param sample_Name 样本名称， 

merge.Obj.direct = function(objPath,sample_Name,projectName = NULL){
    objList = list()
    for( i in sample_Name){
        objList[[i]] = readRDS(paste0(objPath,i,'.RDS'))
    }
    combineObj <- Merge.Seurat.obj(data.list = objList,projectName = projectName,sampleName = sample_Name)
    return(combineObj)
}


#' @description
#' 合并obj对象
#' @param data.list list 对象，每一个list为一个obj对象
#' @param projectName 项目名称，该名称会添加在合并后的barcode前面
#' @param sampleName 样本名称，该名称会添加在合并后的barcode前面

Merge.Seurat.obj <- function(data.list,projectName = "Project",sampleName){
    if(length(data.list)>1){
        obj = suppressMessages(merge(data.list[[1]],data.list[2:length(data.list)],add.cell.ids = paste0(projectName,'_',sampleName)))    
    }else{
        obj = data.list[[1]] 
    };
    return(obj)
}


#' @description
#' 合并obj对象
#' @param data.list list 对象，每一个list为一个obj对象
#' @param projectName 项目名称，该名称会添加在合并后的barcode前面
#' @param sampleName 样本名称，该名称会添加在合并后的barcode前面

.AddGeneAnnota <- function(object,Species){
   
    
    if(Species =='hg38'){
        library(EnsDb.Hsapiens.v86)
        ensdb = EnsDb.Hsapiens.v86
    }else if(Species =='mm10'){
        library(EnsDb.Mmusculus.v79)
        ensdb = EnsDb.Mmusculus.v79
    }
    annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
    seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
    genome(annotation) <- 'GRCh38'
    object@assays[["ATAC"]]@annotation <- annotation
    return(object)
}



#' @description getGeneActivity
#' 计算基因活性
#' @param obj 
#' 
getGeneActivity <- function(object,assay = 'ATAC'){
    gene.activities <- GeneActivity(object,assay =assay )
    # add the gene activity matrix to the Seurat object as a new assay and normalize it
    object[['RNA']] <- CreateAssayObject(counts = gene.activities)
    object <- NormalizeData(
        object = object,
        assay = 'RNA',
        normalization.method = 'LogNormalize',
        scale.factor = median(object$nCount_RNA))
    return(object)
}


