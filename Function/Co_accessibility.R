library(Seurat)
library(Signac)
library(SeuratWrappers)
library(cicero)
library(monocle3)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(ggthemes)
library(ArchR)


#' @description 
#' @param object
#' @param assay
#' @param genome default genome = 'all', 如果特定分析某个/几个染色体上的peaks共可及性情况，则可使genome = c('chr1','chr2')
#' @param 
#' @return a GRanges about peaks to peaks links, which can add to seurat object by function: Links(object) = links.gr
CalculateCoAccConSelf <- function(object,
                                  assay = 'peaks',
                                  genome = 'all',
                                  wind=5e+05){
  
    # convert to CellDataSet format and make the cicero object
    cds <- as.cell_data_set(x = object)
    cds <- monocle3::detect_genes(cds)
    cds <- cds[Matrix::rowSums(exprs(cds)) != 0,] 
    cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
    cicero_cds <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$UMAP)
    
    ## Find Cicero connections
    genomeAll <- seqlengths(object@assays[[assay]]@annotation)
    
    if(length(genome) == 1){
        genome = genomeAll
    }else{
        genome = genomeAll[genome]
    }
    
    genome.df <- data.frame("chr" = names(genome), "length" = genome)
    conns <- run_cicero(cicero_cds, genomic_coords = genome.df, sample_num = 100,window =wind)
    
    ## Find cis-co-accessible networks (CCANs)
    CCANs <- generate_ccans(conns)
    links.gr <- ConnectionsToLinks(conns = conns, ccans = CCANs)
    return(list('links.gr'=links.gr,
        'CCANs'=CCANs,
        'conns'=conns))
}

