# change Fragments path'

ChangeFragmentsPath <- function(object,assay = 'ATAC',old,new){
    n = length(object@assays[[assay]]@fragments)
    for(i in 1:n){
        path = object@assays[[assay]]@fragments[[i]]@path
        path = sub(old,new,path)
        object@assays[[assay]]@fragments[[i]]@path = path
    }
    return(object)
}

