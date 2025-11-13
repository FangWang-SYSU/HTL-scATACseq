library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(ggsci)
library(ComplexHeatmap)
library(JASPAR2024)
library(TFBSTools)
library(patchwork)
library(reshape2)
library(viridis)
library(ggpubr)
library(cowplot)
library(tidyr)


#' @description 
#' @param obj 
#' 
AddMotifSelf <- function(obj,Species = 'hg38',assay = NULL){
    library(JASPAR2024)
    library(TFBSTools)
    library(patchwork)
    
    if(Species =='hg38'){
        library(BSgenome.Hsapiens.UCSC.hg38)
        genome.s = BSgenome.Hsapiens.UCSC.hg38
        species = 9606
    }else if(Species =='mm10'){
        library(BSgenome.Mmusculus.UCSC.mm10)
        genome.s = BSgenome.Mmusculus.UCSC.mm10
        species = 10090
    }
    if(is.null(assay)){
        assay = DefaultAssay(obj)
    }
    mainChrom = standardChromosomes(genome.s)
    keep.peaks <- which(as.character(seqnames(granges(obj))) %in% mainChrom)
    obj[[assay]] <- subset(obj[[assay]], features = rownames(obj[[assay]])[keep.peaks])
    
    jaspar <- JASPAR2024::JASPAR2024()
    pfm <- getMatrixSet(
        x = jaspar@db,opts = list(collection = "CORE",species = species , all_versions = FALSE))
    obj <- AddMotifs(object = obj,genome = genome.s,pfm = pfm,verbose = F)
    return(obj)
    
} 


AddChromVARSelf <- function(obj,Species = 'hg38'){
    if(Species =='hg38'){
        library(BSgenome.Hsapiens.UCSC.hg38)
        genome.s = BSgenome.Hsapiens.UCSC.hg38
    }else if(Species =='mm10'){
        library(BSgenome.Mmusculus.UCSC.mm10)
        genome.s = BSgenome.Mmusculus.UCSC.mm10
    }
    
    obj_chromVAR <- RunChromVAR(object = obj,genome = genome.s)
    return(obj_chromVAR)
    
}




#' @description 
#' @param object seurat object
#' @param assay default: 'peaks', must be a Chrom object in seurat object that which include motif matrix
#' @param foreground.peaks
#' @param background.peaks
#' @return enriched.motifs
MotifEnrichSelf <- function(object,
                            assay = 'peaks',
                            method = 'hyper',
                             foreground.peaks = NULL,
                             background.peaks = NULL
                            ){
    if(is.null(background.peaks)){
        # test enrichment
        enriched.motifs <- FindMotifs(object = object, features = foreground.peaks)
    }else{
        meta.feature <- GetAssayData(object, assay = assay, layer = "meta.features")
        peaks.matched <- MatchRegionStats(meta.feature = meta.feature[background.peaks, ],
                                          query.feature = meta.feature[foreground.peaks, ],
                                          n = 50000)
        # test enrichment
        enriched.motifs <- FindMotifs(object = object, features = foreground.peaks, background = peaks.matched,method = method)
    }
    return(enriched.motifs)
}

SetIfNull <- function(x, y) {
    if (is.null(x = x)) {
        return(y)
    } else {
        return(x)
    }
}

FindMotifs <- function(
        object,
        features,
        background = 40000,
        assay = NULL,
        verbose = TRUE,
        method = 'hyper',
        p.adjust.method = "BH",
        ...
) {
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    background <- SetIfNull(x = background, y = rownames(x = object))
    if (is(object = background, class2 = "numeric")) {
        if (verbose) {
            message("Selecting background regions to match input ",
                    "sequence characteristics")
        }
        meta.feature <- GetAssayData(
            object = object,
            assay = assay,
            slot = "meta.features"
        )
        mf.choose <- meta.feature[
            setdiff(x = rownames(x = meta.feature), y = features), , drop = FALSE
        ]
        missing.features <- setdiff(x = features, y = rownames(x = meta.feature))
        if (length(x = missing.features) > 0) {
            warning(
                "The following features were not found in the assay: ",
                missing.features,
                "\nRemoving missing features", immediate. = TRUE)
            features <- intersect(x = features, y = rownames(x = meta.feature))
        }
        mf.query <- meta.feature[features, , drop = FALSE]
        background <- MatchRegionStats(
            meta.feature = mf.choose,
            query.feature = mf.query,
            regions = features,
            n = background,
            verbose = verbose,
            ...
        )
    }
    
    if (verbose) {
        msg <- ifelse(
            test = length(x = features) > 1,
            yes = " regions",
            no = " region"
        )
        message("Testing motif enrichment in ", length(x = features), msg)
    }
    if (length(x = features) < 10) {
        warning("Testing motif enrichment using a small number of regions is ",
                "not recommended")
    }
    
    motif.all <- GetMotifData(
        object = object, assay = assay, slot = "data"
    )
    
    motif.names <- GetMotifData(
        object = object, assay = assay, slot = "motif.names"
    )
    
    query.motifs <- motif.all[features, , drop = FALSE]
    background.motifs <- motif.all[background, , drop = FALSE]
    query.counts <- colSums(x = query.motifs)
    background.counts <- colSums(x = background.motifs)
    percent.observed <- query.counts / length(x = features) * 100
    percent.background <- background.counts / length(x = background) * 100
    fold.enrichment <- percent.observed / percent.background
    p.list <- vector(mode = "numeric")
    
    # for (i in seq_along(along.with = query.counts)) {
    #     p.list[[i]] <- phyper(
    #         q = query.counts[[i]] - 1,
    #         m = background.counts[[i]],
    #         n = nrow(x = background.motifs) - background.counts[[i]],
    #         k = length(x = features),
    #         lower.tail = FALSE
    #     )
    # }
    
    if (method == 'hyper') {
        for (i in seq_along(along.with = query.counts)) {
            p <- -phyper(
                q = query.counts[[i]] - 1,
                m = background.counts[[i]],
                n = nrow(x = background.motifs) - background.counts[[i]],
                k = length(x = features),
                lower.tail = FALSE,
                log.p = TRUE)
            p = p/log(10)
            p.list[[i]] = p
        }
        
    }else if (method == 'fisher') {
        for (i in seq_along(along.with = query.counts)) {
            matrix = matrix(c(query.counts[[i]],length(x = features) -query.counts[[i]],background.counts[[i]],nrow(x = background.motifs) - background.counts[[i]]), nrow = 2, byrow = TRUE)
            rownames(matrix) <- c("Foreground", "Background")
            colnames(matrix) <- c("Target Gene", "Non-target Gene")
            # 进行双尾Fisher精确检验
            fisher.result <- fisher.test(matrix, alternative = "two.sided")
            p.list[[i]] = -log10(fisher.result[["p.value"]])
        }
    }
    
    
    results <- data.frame(
        motif = names(x = query.counts),
        observed = query.counts,
        background = background.counts,
        percent.observed = percent.observed,
        percent.background = percent.background,
        fold.enrichment = fold.enrichment,
        pvalue =  10^(-p.list),
        mlog10p = p.list,
        motif.name = as.vector(
            x = unlist(x = motif.names[names(x = query.counts)])
        ),
        #p.adjust = p.adjust(p = p.list, method = p.adjust.method),
        stringsAsFactors = FALSE
    )
    results$mlog10Padj <- pmax(results$mlog10p - log10(ncol(results)), 0)
    results$p.adjust <- 10^(-results$mlog10Padj)
    
    if (nrow(x = results) == 0) {
        return(results)
    } else {
        return(results[order(results[, 7], -results[, 6]), ])
    }
}



select.motif <- function(enriched.motifs,
                         n = 20,
                         p.cut = 20,
                         MotifRef = NULL){
    
    motif.mat = enriched.motifs[,c('motif','fold.enrichment','p.adjust','Cluster', 'motif.name')]
    length(unique(motif.mat$motif.name))
    
    motif.mat$`-log10P` <- -log10(motif.mat$p.adjust)
    motif.mat[which(motif.mat$`-log10P`  > 100),'-log10P'] <- 100
    motif.mat$`-log10P` =  ifelse(motif.mat$fold.enrichment<1,-motif.mat$`-log10P`,motif.mat$`-log10P`)
    #### //在所有组中都显著富集的motif
    motif.mat2 = subset(motif.mat,`-log10P` >= 90)
    nClus = length(unique(motif.mat2$Cluster))
    temp <- as.data.frame(table(motif.mat2$motif.name))
    temp = temp[temp$Freq<4, ]
    motif.mat = motif.mat[motif.mat$motif.name%in%temp$Var1,]
    
    
    ## 选择
    top.motifMat <- motif.mat %>% group_by(motif.name) %>% slice_max(n=1,order_by = `-log10P`,with_ties = F)
    #top.motifMat = sig.motifMat
    motif = top.motifMat%>%group_by(Cluster)%>%slice_max(n=n ,order_by = fold.enrichment, with_ties = F)
    
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

.rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
    z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
    if(limit){
        z[z > max] <- max
        z[z < min] <- min
    }
    return(z)
}
