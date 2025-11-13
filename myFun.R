suppressMessages(suppressWarnings({options(warn=-1)}))
suppressMessages(suppressWarnings({library(Signac)}))
suppressMessages(suppressWarnings({library(Seurat)}))
suppressMessages(suppressWarnings({library(glue)}))
suppressMessages(suppressWarnings({library(EnsDb.Mmusculus.v79)}))
suppressMessages(suppressWarnings({library(ggsci)}))
suppressMessages(suppressWarnings({base_dir='/Volumes/wangxin/鼠肺项目/code/'}))
suppressMessages(suppressWarnings({CodePath = glue('{base_dir}/ATAC_pipeline_zyr/Code/')}))
suppressMessages(suppressWarnings({source(file.path(CodePath,'Function/PipelineFuncHtml.R'))}))
suppressMessages(suppressWarnings({source(file.path(CodePath,'Function/LoadsingleSample.func.R'))}))
suppressMessages(suppressWarnings({source(file.path(CodePath,'Function/MotifAnalysis.func.R'))}))
suppressMessages(suppressWarnings({source(file.path(CodePath,'Function/Signac.Heatmap.R'))}))
suppressMessages(suppressWarnings({source(file.path(CodePath,'Function/ClusterAndTypeing.func.R'))}))
suppressMessages(suppressWarnings({source(file.path(CodePath,'Function/cCREsStatistics.func.R'))}))
suppressMessages(suppressWarnings({source(file.path(CodePath,'Function/Signac.Heatmap.R'))}))
suppressMessages(suppressWarnings({source(file.path(CodePath,'Function/MergeSampleObj.R'))}))
suppressMessages(suppressWarnings({source(file.path(CodePath,'Function/PredictCellType.R'))}))
suppressMessages(suppressWarnings({library(patchwork)}))
suppressMessages(suppressWarnings({library(ggrepel)}))
suppressMessages(suppressWarnings({library(dplyr)}))
suppressMessages(suppressWarnings({library(ggplot2)}))

maintype_color=c('Immune'='#D7191C',
                 'Epithelium'='#1F78B4',
                 'Endothelium'='#FEE08B',
                 'Stromal'='#A6761D', 
                 'CNS'='#7570B3')
col_map_func <- function(dat){
  res = list()
  for(n in colnames(dat)){
    tmp_name = sort(unique(dat[, n]))
    tmp_col = pal_igv()(length(tmp_name))
    names(tmp_col) = tmp_name
    res[[n]] = tmp_col
  }
  return(res)
}
give.n = function(x) {
  la.df = data.frame(y = mean(x)+ max(x) / 10,label = round(mean(x), 2));
  return(la.df)
}
my_process_srt <- function(object, res=0.8, dim=30,min.cutoff='q5'){
  #DefaultAssay(object) <- 'ATAC'
  object <- RunTFIDF(object)
  object <- FindTopFeatures(object, min.cutoff = min.cutoff,verbose = F)
  object <- RunSVD(object,verbose = F)
  object <- RunUMAP(object = object, reduction = 'lsi', dims = 2:dim)
  object <- FindNeighbors(object = object, reduction = 'lsi', dims = 2:dim )
  object <- FindClusters(object = object, verbose = FALSE, algorithm = 3,resolution=res)
  return(object)
}
sample_levels = c('0709_P0',
                  '1009_P0',
                  '0726_P1',
                  '0815_P1',
                  '0828_P2',
                  '0606_P3',
                  '0903_P4',
                  '0819_P5',
                  '0731_P6',
                  '0821_P7',
                  '0425_P8',
                  '1017_P8',
                  '0823_P9',
                  '1018_P9',
                  '0905_P10',
                  '0906_P11',
                  '0806_P12',
                  '0722_P13',
                  '1022_P13',
                  '0617_P14')
time_color = c(
  '0709_P0' ='#313695',
  '1009_P0' ='#313695',
  '0726_P1' ='#4575B4',
  '0815_P1' ='#4575B4',
  '0828_P2' ='#619CC3',
  '0606_P3' ='#74ADD1',
  '0903_P4' ='#89C4DC',
  '0819_P5' ='#ABD9E9',
  '0731_P6' ='#D0EBF5',
  '0821_P7' ='#E0F3F8',
  '0425_P8' ='#FFFFBF',
  '1017_P8' ='#FFFFBF',
  '0823_P9' ='#FEE090',
  '1018_P9' ='#FEE090',
  '0905_P10'='#FDBB84',
  '0906_P11'='#FDAE61',
  '0806_P12'='#F46D43',
  '0722_P13'='#D73027',
  '1022_P13'='#D73027',
  '0617_P14'='#A50026')

time_color2 = c(
  'P0' ='#313695',
  'P1' ='#4575B4',
  'P2' ='#619CC3',
  'P3' ='#74ADD1',
  'P4' ='#89C4DC',
  'P5' ='#ABD9E9',
  'P6' ='#D0EBF5',
  'P7' ='#E0F3F8',
  'P8' ='#FFFFBF',
  'P9' ='#FEE090',
  'P10'='#FDBB84',
  'P11'='#FDAE61',
  'P12'='#F46D43',
  'P13'='#D73027',
  'P14'='#A50026')
time_levels = paste0('P', 0:14)

myFindRegion <- function (object, region, sep = c("-", "-"), assay = NULL, extend.upstream = 0, 
                         extend.downstream = 0) 
{
  if (!is(object = region, class2 = "GRanges")) {
    region <- tryCatch(expr = suppressWarnings(expr = StringToGRanges(regions = region, 
                                                                      sep = sep)), error = function(x) {
                                                                        region <- LookupGeneCoords(object = object, assay = assay, 
                                                                                                   gene = region)
                                                                        return(region)
                                                                      })
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(x = region, upstream = extend.upstream, 
                                           downstream = extend.downstream))
  return(region)
}

myCoveragePlotSingle <- function(obj_cells,obj_peaks, region, 
                                     extend.upstream = 0, 
                                     extend.downstream = 0,
                                     ymax=NULL,
                                     window=100,heights = NULL){
  # plot
  rg <- myFindRegion(object = obj_cells, region = region, 
                     sep = c("-", "-"),
                     extend.upstream = extend.upstream,
                     extend.downstream = extend.downstream)
  a1=CoveragePlot(obj_cells, 
                  region = rg,
                  ymax=ymax,
                  window=window, 
                  extend.upstream = 0, 
                  extend.downstream = 0,
                  annotation = FALSE,
                  peaks = FALSE)

  a2=AnnotationPlot(obj_cells, rg, extend.upstream = 0, extend.downstream = 0)
  a3=PeakPlot(obj_peaks, region=rg)
  return(CombineTracks(plotlist = list(a1,a2,a3),heights = NULL))
}


# myCoveragePlotByBwSingle <- function(obj, region, bigwigs,
#                                extend.upstream = 0, 
#                                extend.downstream = 0,
#                                window=100,heights = NULL,ymax=NULL){
#   # plot
#   rg <- myFindRegion(object = obj, region = region, 
#                          sep = c("-", "-"),
#                          extend.upstream = extend.upstream,
#                      extend.downstream = extend.downstream)
#   a1=BigwigTrack(region=rg,
#                  bigwig=bigwigs,
#                  smooth = window,
#                  ymax=ymax,
#                  extend.upstream = extend.upstream,
#                  extend.downstream = extend.downstream,
#                  type='coverage'
#   )
#   a1$data[,'group'] = a1$data$bw
#   a1$data[,'coverage'] = a1$data$score
#   a2=AnnotationPlot(obj, rg,extend.upstream = 0, extend.downstream = 0)
#   a3=PeakPlot(obj, region=rg)
#   a1 = a1 + theme(legend.position = 'none', axis.text.y = element_blank())
#   a2 = a2 + theme(legend.position = 'none')
#   a3 = a3 + theme(legend.position = 'none')
#   return(CombineTracks(plotlist = list(a1,a2,a3),heights = NULL))
# }

myCoveragePlotMultiple <- function(obj_cells,obj_peaks, region, 
                               ymax=NULL,window=2000, heights = NULL,
                               extend.upstream = 0, 
                               extend.downstream = 0,...){
  if (length(x = region) == 1) {
    region <- list(region)
  }
  plot.list = lapply(X = seq_along(region), FUN = function(x) {
    myCoveragePlotSingle(obj_cells,obj_peaks, region[[x]], 
                             extend.upstream = extend.upstream, 
                             extend.downstream = extend.downstream,
                             ymax=ymax,window=window,heights = heights)
  })
  return(wrap_plots(plot.list, ...))
}

myCoverageModify<-function(plot_res, 
                           region,
                           ncol=NULL,
                           group_color=NULL, 
                           genes_color='black', 
                           peaks_color='red',
                           line_color_with_group=TRUE,...){
  if(is.null(ncol)){
    row_first = c(1)
  }else{
    row_first = seq(1,length(region),ncol)
  }
  if(is.null(group_color)){
    fill_theme = scale_fill_igv()
    color_theme = scale_color_igv()
  }else{
    fill_theme = scale_fill_manual(group_color)
    color_theme = scale_color_manual(group_color)
  }
  for(gene_plot in 1:length(plot_res)){
    tmp_layer = plot_res[[gene_plot]][[2]]$layers
    gene_title=region[gene_plot]
    for(ix in 1:length(tmp_layer)){
      if('GeomText'%in%class(tmp_layer[[ix]]$geom)){
        y = tmp_layer[[ix]]$data
        gene_title = y[which.max(y$width), 'gene_name']
      }
    }
    # gene_names = plot_res[[gene_plot]][[2]]$layers[[4]]$data$gene_name
    # gene_names = intersect(region, gene_names)
    # if(length(gene_names)>=1){
    #   gene_title=gene_names
    # }else{
    #   gene_title=''
    # }
    #print(gene_title)
    if(gene_plot%in%row_first){
      # coverage 1
      y_upper = max(plot_res[[gene_plot]][[1]]$data$coverage)
      if(line_color_with_group){
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(aes(yintercept=0,color=group), size=0.8)+
          color_theme+
          labs(title=gene_title)+
          #theme_void()+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            strip.placement = 'outside',
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            strip.text = element_text(size=8, face="bold"),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }else{
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(yintercept=0, size=0.5)+
          labs(title=gene_title)+
          theme(
            panel.spacing.y = unit(x = 0, units = "line"),
            strip.placement = 'outside',
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_text(size=8, face="bold"),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }
      # gene
      a2 = plot_res[[gene_plot]][[2]]
      label_ix = 4
      for(ix in 1:length(a2$layers)){
        if('GeomText'%in%class(a2$layers[[ix]]$geom)){
          label_ix=ix
        }else{
          if(ix==1){
            a2$layers[[ix]]$aes_params$size=2
          }else{
            a2$layers[[ix]]$aes_params$size=0.1
          }
        }
        
      }
      # a2$layers[[1]]$aes_params$size=2
      # a2$layers[[2]]$aes_params$size=0.1
      # a2$layers[[3]]$aes_params$size=0.1
      a2$layers[[label_ix]] = NULL
      #gene_y_lims = layer_scales(a2)$y$range$range
      tmp_da = a2$layer[[1]]$data
      dodge_y = as.numeric(tmp_da[tmp_da$gene_name==gene_title, 'dodge'][1])
      gene_y_lims = c(dodge_y-0.1,dodge_y+0.1)
      a2=a2+
        scale_color_manual(values = c(genes_color))+
        scale_size_manual(values=0.1)+
        scale_y_continuous(limits = gene_y_lims)+
        theme(panel.border = element_rect(fill=NA, size=0.5),
              strip.placement = 'inside',
              axis.line = element_blank(),
              axis.title.y = element_text(angle = 0, vjust = 0.5,
                                          size=8, face="bold",
                                          margin=margin(l=0,r=-2,t=0,b=0, unit = 'line'))
        )
      
      plot_res[[gene_plot]][[2]] = a2
      
        
      # peaks 3
      a3 = plot_res[[gene_plot]][[3]]
      if(length(a3$layers)>0){
          a3$layers[[1]]$aes_params$size=3
      }
      a3 =a3+
        scale_color_manual(values = c(peaks_color))+
        #scale_y_continuous(limits = c(-0,0))+
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(fill=NA, size=0.5),
          strip.placement = 'inside',
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5,
                                      size=8, face="bold",
                                      margin=margin(l=0,r=-2,t=0,b=0, unit = 'line'))
        )
      plot_res[[gene_plot]][[3]] = a3
    }else{
      # coverage 1
      y_upper = max(plot_res[[gene_plot]][[1]]$data$coverage)
      if(line_color_with_group){
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(aes(yintercept=0,color=group), size=0.8)+
          color_theme+
          labs(title=gene_title)+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            #strip.placement = 'outside',
            strip.text.y.left = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }else{
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(yintercept=0, size=0.5)+
          labs(title=gene_title)+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            #strip.placement = 'outside',
            strip.text.y.left = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }
      # gene
      a2 = plot_res[[gene_plot]][[2]]
      label_ix = 4
      for(ix in 1:length(a2$layers)){
        if('GeomText'%in%class(a2$layers[[ix]]$geom)){
          label_ix=ix
        }else{
          if(ix==1){
            a2$layers[[ix]]$aes_params$size=2
          }else{
            a2$layers[[ix]]$aes_params$size=0.1
          }
        }
        
      }
      # a2$layers[[1]]$aes_params$size=2
      # a2$layers[[2]]$aes_params$size=0.1
      # a2$layers[[3]]$aes_params$size=0.1
      a2$layers[[label_ix]] = NULL
      #gene_y_lims = layer_scales(a2)$y$range$range
      tmp_da = a2$layer[[1]]$data
      dodge_y = as.numeric(tmp_da[tmp_da$gene_name==gene_title, 'dodge'][1])
      gene_y_lims = c(dodge_y-0.1,dodge_y+0.1)
      a2=a2+
        scale_color_manual(values = c(genes_color))+
        scale_size_manual(values=0.1)+
        scale_y_continuous(limits = gene_y_lims)+
        theme(panel.border = element_rect(fill=NA, size=0.5),
              strip.placement = 'inside',
              axis.line = element_blank(),
              #axis.title.y = element_text(angle = 0, vjust = 0.5,
              #                            margin=margin(l=0,r=-4,t=0,b=0, unit = 'line')),
              axis.title.y=element_blank()
        )
      
      plot_res[[gene_plot]][[2]] = a2
      # peaks 3
      a3 = plot_res[[gene_plot]][[3]]
      if(length(a3$layers)>0){
          a3$layers[[1]]$aes_params$size=3
      }
      a3 =a3+
        scale_color_manual(values = c(peaks_color))+
        #scale_y_continuous(limits = c(-0,0))+
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_rect(fill=NA, size=0.5),
          strip.placement = 'inside',
          axis.line = element_blank(),
          axis.title.y=element_blank()
        )
      plot_res[[gene_plot]][[3]] = a3
    }
  }
  wp = list()
  for(i in 1:length(plot_res)){
    wp[[i]] = plot_res[[i]]
  }
  return(wrap_plots(wp, ncol = ncol))
}
myCoveragePlot <- function(obj, region, group.by, ymax='q75',
                           window=2000, heights = c(20,1,1),
                           ncol=NULL,
                           group_color=NULL, genes_color='black', peaks_color='red',line_color_with_group=TRUE,...
){
  plot_res =CoveragePlot(obj, region = region, 
                         group.by = group.by,
                         ymax=ymax, window=window, heights = heights,ncol=ncol,...)
  if(is.null(ncol)){
    row_first = c(1)
  }else{
    row_first = seq(1,length(region),ncol)
  }
  if(is.null(group_color)){
    fill_theme = scale_fill_igv()
    color_theme = scale_color_igv()
  }else{
    fill_theme = scale_fill_manual(group_color)
    color_theme = scale_color_manual(group_color)
  }
  for(gene_plot in 1:length(plot_res)){
    tmp_layer = plot_res[[gene_plot]][[2]]$layers
    gene_title=region[gene_plot]
    for(ix in 1:length(tmp_layer)){
        if('GeomText'%in%class(tmp_layer[[ix]]$geom)){
          y = tmp_layer[[ix]]$data
          gene_title = y[which.max(y$width), 'gene_name']
        }
    }
    # gene_names = plot_res[[gene_plot]][[2]]$layers[[4]]$data$gene_name
    # gene_names = intersect(region, gene_names)
    # if(length(gene_names)>=1){
    #   gene_title=gene_names
    # }else{
    #   gene_title=''
    # }
    #print(gene_title)
    if(gene_plot%in%row_first){
      # coverage 1
      y_upper = max(plot_res[[gene_plot]][[1]]$data$coverage)
      if(line_color_with_group){
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(aes(yintercept=0,color=group), size=0.8)+
          color_theme+
          labs(title=gene_title)+
          #theme_void()+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            strip.placement = 'outside',
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            strip.text = element_text(size=8, face="bold"),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }else{
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(yintercept=0, size=0.5)+
          labs(title=gene_title)+
          theme(
            panel.spacing.y = unit(x = 0, units = "line"),
            strip.placement = 'outside',
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_text(size=8, face="bold"),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }
      
      # gene
      a2 = plot_res[[gene_plot]][[2]]
      label_ix = 4
      for(ix in 1:length(a2$layers)){
          if('GeomText'%in%class(a2$layers[[ix]]$geom)){
              label_ix=ix
          }else{
              if(ix==1){
                  a2$layers[[ix]]$aes_params$size=2
              }else{
                  a2$layers[[ix]]$aes_params$size=0.1
              }
          }
      
      }
      # a2$layers[[1]]$aes_params$size=2
      # a2$layers[[2]]$aes_params$size=0.1
      # a2$layers[[3]]$aes_params$size=0.1
      a2$layers[[label_ix]] = NULL
      #gene_y_lims = layer_scales(a2)$y$range$range
      tmp_da = a2$layer[[1]]$data
      dodge_y = as.numeric(tmp_da[tmp_da$gene_name==gene_title, 'dodge'][1])
      gene_y_lims = c(dodge_y-0.1,dodge_y+0.1)
      a2=a2+
        scale_color_manual(values = c(genes_color))+
        scale_size_manual(values=0.1)+
        scale_y_continuous(limits = gene_y_lims)+
        theme(panel.border = element_rect(fill=NA, size=0.5),
              strip.placement = 'inside',
              axis.line = element_blank(),
              axis.title.y = element_text(angle = 0, vjust = 0.5,
                                          size=8, face="bold",
                                          margin=margin(l=0,r=-2,t=0,b=0, unit = 'line'))
        )
      
      plot_res[[gene_plot]][[2]] = a2
      
      # peaks 3
      a3 = plot_res[[gene_plot]][[3]]
      a3$layers[[1]]$aes_params$size=3
      a3 =a3+
        scale_color_manual(values = c(peaks_color))+
        #scale_y_continuous(limits = c(-0,0))+
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(fill=NA, size=0.5),
          strip.placement = 'inside',
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5,
                                      size=8, face="bold",
                                      margin=margin(l=0,r=-2,t=0,b=0, unit = 'line'))
        )
      plot_res[[gene_plot]][[3]] = a3
    }else{
      # coverage 1
      y_upper = max(plot_res[[gene_plot]][[1]]$data$coverage)
      if(line_color_with_group){
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(aes(yintercept=0,color=group), size=0.8)+
          color_theme+
          labs(title=gene_title)+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            #strip.placement = 'outside',
            strip.text.y.left = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }else{
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(yintercept=0, size=0.5)+
          labs(title=gene_title)+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            #strip.placement = 'outside',
            strip.text.y.left = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }
      
      # gene
      a2 = plot_res[[gene_plot]][[2]]
      label_ix = 4
      for(ix in 1:length(a2$layers)){
        if('GeomText'%in%class(a2$layers[[ix]]$geom)){
          label_ix=ix
        }else{
          if(ix==1){
            a2$layers[[ix]]$aes_params$size=2
          }else{
            a2$layers[[ix]]$aes_params$size=0.1
          }
        }
        
      }
      # a2$layers[[1]]$aes_params$size=2
      # a2$layers[[2]]$aes_params$size=0.1
      # a2$layers[[3]]$aes_params$size=0.1
      a2$layers[[label_ix]] = NULL
      #gene_y_lims = layer_scales(a2)$y$range$range
      tmp_da = a2$layer[[1]]$data
      dodge_y = as.numeric(tmp_da[tmp_da$gene_name==gene_title, 'dodge'][1])
      gene_y_lims = c(dodge_y-0.1,dodge_y+0.1)
      a2=a2+
        scale_color_manual(values = c(genes_color))+
        scale_size_manual(values=0.1)+
        scale_y_continuous(limits = gene_y_lims)+
        theme(panel.border = element_rect(fill=NA, size=0.5),
              strip.placement = 'inside',
              axis.line = element_blank(),
              #axis.title.y = element_text(angle = 0, vjust = 0.5,
              #                            margin=margin(l=0,r=-4,t=0,b=0, unit = 'line')),
              axis.title.y=element_blank()
        )
      
      plot_res[[gene_plot]][[2]] = a2
      
      # peaks 3
      a3 = plot_res[[gene_plot]][[3]]
      a3$layers[[1]]$aes_params$size=3
      a3 =a3+
        scale_color_manual(values = c(peaks_color))+
        #scale_y_continuous(limits = c(-0,0))+
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_rect(fill=NA, size=0.5),
          strip.placement = 'inside',
          axis.line = element_blank(),
          axis.title.y=element_blank()
        )
      plot_res[[gene_plot]][[3]] = a3
    }
  }
  wp = list()
  for(i in 1:length(plot_res)){
    wp[[i]] = plot_res[[i]]
  }
  return(wrap_plots(wp, ncol = ncol))
}


myHeatmapPeaks <- function(marker.peaks, obj=NULL, group='seurat_clusters',
                      top=10,assay='ATAC',
                      min_log2FC=0.25,max_pval=0.05,
                      color=circlize::colorRamp2(c(-2,0,2),c("#045a8d","white","#a50f15")),
                      row_fontsize=3,column_fontsize=3,
                      cell_meta=NULL,mat=NULL,filter_pseudogene=FALSE,
                           lwd=0.05
                      ){
  if(!is.null(obj)){
    cell_meta = obj@meta.data
    mat=obj@assays[[assay]]$data
  }
  if(filter_pseudogene){
      marker.peaks = filter_pseudoGene(marker.peaks)
    
  }

  top_peaks = marker.peaks %>% 
    dplyr::filter(avg_log2FC >=min_log2FC&p_val_adj <=max_pval) %>%
    group_by(cluster) %>% 
    top_n(top,avg_log2FC) %>%
    arrange(cluster, -avg_log2FC)
  
  cell_meta = cell_meta[order(cell_meta[,group]),]
  top_mat = mat[top_peaks$gene, rownames(cell_meta)]
  top_mat = t(apply(top_mat, 1, scale))
  ht=Heatmap(as.matrix(top_mat),
            col=color,
            cluster_rows = F, cluster_columns = F,
            show_row_names = F, show_column_names = F,
            #top_annotation = HeatmapAnnotation(df=cell_meta$leiden),
            column_split = cell_meta[,group],
            row_split = top_peaks$cluster,
            use_raster=FALSE,
            column_gap = unit(0.0, 'mm'),
            row_gap = unit(0.0, 'mm'),
            border = TRUE,
            border_gp = gpar(col = "black",lwd = lwd),
            row_title_gp = gpar(fontsize=row_fontsize),
            column_title_gp = gpar(fontsize=column_fontsize)
  )
  return(ht)
}

myHeatmapPeaksGenes <- function(marker.peaks, obj=NULL, group='seurat_clusters',
                                top=10,assay='ATAC',
                                min_log2FC=0.25,max_pval=0.05,
                                color=circlize::colorRamp2(c(-2,0,2),c("#045a8d","white","#a50f15")),
                                peak.distance=NULL,
                                region_type='all',
                                gene_class = NULL,
                                gene_col_na = 'gray',
                                gene_name_fontsize=2,
                                row_fontsize=12, column_fontsize=12,
                                show_gene_name_by_class=FALSE,
                                cell_meta=NULL,mat=NULL,
                                lwd=0.05,
                                filter_pseudoGene=TRUE
){
    if(!is.null(obj)){
    cell_meta = obj@meta.data
    mat=obj@assays[[assay]]$data
  }

  #cell_meta = obj@meta.data
  #mat=obj@assays[[assay]]$data
  # filter
  if(filter_pseudoGene){
      marker.peaks = filter_pseudoGene(marker.peaks)
  }
  
  if(region_type=='gene_body'){
    marker.peaks = subset(marker.peaks,isin_body == 1)
  }else if(region_type=='gene_promoter'){
    marker.peaks = subset(marker.peaks,isin_promoter == 1)
  }else if(region_type=='genes'){
    marker.peaks = subset(marker.peaks,isin_promoter == 1 | isin_body == 1)
  }
  if(!is.null(peak.distance)){
    marker.peaks = subset(marker.peaks,distance <= peak.distance)
  }
  top_peaks = marker.peaks %>% 
    dplyr::filter(avg_log2FC >=min_log2FC&p_val_adj <=max_pval) %>%
    group_by(cluster) %>% 
    top_n(top,avg_log2FC) %>%
    arrange(cluster, -avg_log2FC)
  
  cell_meta = cell_meta[order(cell_meta[,group]),]
  top_peaks$cluster = factor(top_peaks$cluster, 
                           levels = cell_meta[,group][!duplicated(cell_meta[,group])])
  top_peaks = top_peaks %>% arrange(cluster)

  top_mat = mat[top_peaks$gene, rownames(cell_meta)]
  top_mat = t(apply(top_mat, 1, scale))
  
  if(show_gene_name_by_class){
    rownames(top_mat) = paste0(top_peaks$gene_name,'-->', gene_class[top_peaks$gene_name])
    show_gene_col_legend=FALSE
  }else{
    rownames(top_mat) = top_peaks$gene_name
    show_gene_col_legend = TRUE
  }
  lg = list(Legend(title = "Peaks",
                   col_fun = circlize::colorRamp2(c(-2,0,2),c("#045a8d","white","#a50f15")),
                   grid_height = unit(2, "mm"),
                   grid_width = unit(2, "mm"),
                   labels_gp = gpar(fontsize = 6),
                   title_gp = gpar(fontsize = 7, fontface = "bold"))
  )
  if(is.null(gene_class)){
    gene_col = rep('black', nrow(top_mat))
  }else{
    tmp_gene_class = gene_class[top_peaks$gene_name]
    tmp_vec = unique(na.omit(tmp_gene_class))
    gene_name_col=pal_igv()(length(tmp_vec))
    names(gene_name_col) = unique(tmp_vec)
    
    gene_col = gene_name_col[tmp_gene_class]
    gene_col[is.na(gene_col)] = gene_col_na
    if(show_gene_col_legend){
      lg[[2]]=Legend(labels = names(gene_col), 
                     title = "Class", type = "points", 
                     legend_gp = gpar(col = gene_col),
                     background = "white",
                     grid_height = unit(2, "mm"),
                     grid_width = unit(2, "mm"),
                     labels_gp = gpar(fontsize = 6),
                     title_gp = gpar(fontsize = 7, fontface = "bold"))
    }
  }
  
  ht=Heatmap(top_mat,
             col=color,
             cluster_rows = F, cluster_columns = F,
             show_row_names = T, show_column_names = F,
             #top_annotation = HeatmapAnnotation(df=cell_meta$leiden),
             column_split = cell_meta[,group],
             row_split = top_peaks$cluster,
             use_raster=FALSE,
             column_gap = unit(0.0, 'mm'),
             row_gap = unit(0.0, 'mm'),
             row_names_gp = gpar(col=gene_col,fontsize=gene_name_fontsize),
             border = TRUE,
             show_heatmap_legend = F,
             border_gp = gpar(col = "black",lwd = lwd),
             row_title_gp = gpar(fontsize=row_fontsize),
             column_title_gp = gpar(fontsize=column_fontsize)
  )
  return(list('ht'=ht,'lgd'=lg))
}

get_promoter <- function(object=EnsDb.Hsapiens.v86,upstream = 2000,downstream = 0){
  gene.ranges <- genes(object)
  gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
  gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')
  gene.promoters <- promoters(gene.ranges, upstream = upstream,downstream = downstream)
  return(gene.promoters)
}
get_genebody <- function(object=EnsDb.Hsapiens.v86){
  gene.ranges <- genes(object)
  gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
  gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')
  return(gene.ranges)
}

add_gene_info <- function(obj, marker.peaks, geonome=EnsDb.Mmusculus.v79, promoter_range=2000){
  # 1.添加最近距离的基因
  uniquePeaks = unique(marker.peaks$gene)
  closest_genes <- ClosestFeature(obj, regions = uniquePeaks)
  marker.peaks_gene=left_join(marker.peaks, closest_genes, by = join_by('gene' == 'query_region'))
  
  # 2.添加所在启动子区域的基因
  promoter = get_promoter(geonome,upstream = promoter_range,downstream = 0)
  seqlevels(promoter) <- paste0('chr', seqlevels(promoter))
  closest_genes_promoter <- ClosestFeature(obj, regions = uniquePeaks, annotation = promoter)
  closest_genes_promoter = closest_genes_promoter[closest_genes_promoter$distance==0,]
  marker.peaks_gene$isin_promoter=0
  marker.peaks_gene[marker.peaks_gene$gene_name%in%closest_genes_promoter$gene_name, 'isin_promoter'] = 1
  
  # 3.添加基因body区域基因
  genebody = get_genebody(geonome)
  seqlevels(genebody) <- paste0('chr', seqlevels(genebody))
  closest_genes_body <- ClosestFeature(obj, regions = uniquePeaks, annotation = genebody)
  closest_genes_body = closest_genes_body[closest_genes_body$distance==0,]
  marker.peaks_gene$isin_body=0
  marker.peaks_gene[marker.peaks_gene$gene_name%in%closest_genes_body$gene_name, 'isin_body'] = 1
  return(marker.peaks_gene)
}

filter_pseudoGene <- function(marker.peaks){
  marker.peaks = marker.peaks[!grepl('^Gm', marker.peaks$gene_name),]
  marker.peaks = marker.peaks[!grepl('^n-', marker.peaks$gene_name),]
  marker.peaks = marker.peaks[!grepl('^[0-9]', marker.peaks$gene_name),]
  marker.peaks = marker.peaks[!grepl('^AU[0-9]', marker.peaks$gene_name),]
  marker.peaks = marker.peaks[!grepl('^BC[0-9]', marker.peaks$gene_name),]
  marker.peaks = marker.peaks[!grepl('\\|', marker.peaks$gene_name),]
  return(marker.peaks)
}

myDimPlot<-function(umap, groupby='leiden', label=TRUE, group_color=NULL, point_size=0.3){
    if(is.null(group_color)){
        group_color = ggsci::pal_igv()(50)
    }
    scale_color = scale_color_manual(values=group_color)
    scale_fill = scale_fill_manual(values=group_color)
    
    if(label){
      label_pos = umap %>% 
          dplyr::group_by(!!sym(groupby)) %>%
          dplyr::summarise(UMAP_1=median(UMAP_1), UMAP_2=median(UMAP_2)) %>%
          as.data.frame()
      a=umap %>%
          ggplot(aes(x=UMAP_1, y=UMAP_2))+
          geom_point(aes_string(fill=groupby),shape=21,size=point_size, stroke=NA)+
          geom_label_repel(data=label_pos, mapping = aes_string(label=groupby,fill=groupby))+
          scale_fill+
          theme_classic()+
          theme(
            axis.line = element_blank(),
            axis.ticks = element_blank(), # 刻度不显示
            axis.text = element_blank(), # 刻度text不显示
            axis.title = element_blank(),
            legend.background = element_rect(fill = "white", size = 1, colour = "white"),
          )+
          guides(fill = guide_legend(override.aes = list(size=4, stroke=NA,shape=21)))+ # legend大小
          labs('title'='')+
          geom_segment(aes(x=min(umap$UMAP_1), y=min(umap$UMAP_2), xend=min(umap$UMAP_1)+2, yend=min(umap$UMAP_2)),
                       colour="black", size=0.5,arrow = arrow(length=unit(0.2,"cm")))+ 
          geom_segment(aes(x = min(umap$UMAP_1), y = min(umap$UMAP_2), xend = min(umap$UMAP_1), yend=min(umap$UMAP_2)+2),
                       colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"))) +
          annotate("text", x = min(umap$UMAP_1)+4, y = min(umap$UMAP_2), label = "UMAP1",
                   color="black",size = 4) + 
          annotate("text", x = min(umap$UMAP_1), y = min(umap$UMAP_2)+4, label = "UMAP2",
                   color="black",size = 4, angle=90) 
    }else{
        a=umap %>%
          ggplot(aes(x=UMAP_1, y=UMAP_2))+
          geom_point(aes_string(fill=groupby),shape=21,size=point_size, stroke=NA)+
          scale_fill+
          theme_classic()+
          theme(
            axis.line = element_blank(),
            axis.ticks = element_blank(), # 刻度不显示
            axis.text = element_blank(), # 刻度text不显示
            axis.title = element_blank(),
            legend.background = element_rect(fill = "white", size = 1, colour = "white"),
          )+
          guides(fill = guide_legend(override.aes = list(size=4, stroke=NA,shape=21)))+ # legend大小
          labs('title'='')+
          geom_segment(aes(x=min(umap$UMAP_1), y=min(umap$UMAP_2), xend=min(umap$UMAP_1)+2, yend=min(umap$UMAP_2)),
                       colour="black", size=0.5,arrow = arrow(length=unit(0.2,"cm")))+ 
          geom_segment(aes(x = min(umap$UMAP_1), y = min(umap$UMAP_2), xend = min(umap$UMAP_1), yend=min(umap$UMAP_2)+2),
                       colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"))) +
          annotate("text", x = min(umap$UMAP_1)+4, y = min(umap$UMAP_2), label = "UMAP1",
                   color="black",size = 4) + 
          annotate("text", x = min(umap$UMAP_1), y = min(umap$UMAP_2)+4, label = "UMAP2",
                   color="black",size = 4, angle=90) 
    }
    return(a)
}



myCoverageModify_Heat<-function(plot_res,
                                region_title,
                                axis_x_title=NULL,
                                group_color=NULL, 
                                title_angle=0, 
                                title_hjust=0.5, 
                                title_size=12,
                               group_size=12,
                               axis_x_title_size=12){
    if(is.null(group_color)){
      fill_theme = scale_fill_igv()
      color_theme = scale_color_igv()
    }else{
      fill_theme = scale_fill_manual(values=group_color)
      color_theme = scale_color_manual(values=group_color)
    }
    if(is.null(axis_x_title)){
        axis_x_title = rep(NULL, length(plot_res))
    }
    for(gene_plot in 1:length(plot_res)){
        gene_title = region_title[gene_plot]
        y_upper = max(plot_res[[gene_plot]]$data$coverage)
        if(gene_plot==1){
            plot_res[[gene_plot]]=plot_res[[gene_plot]]+
              fill_theme+
              scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                                 expand=c(0,0))+
              geom_hline(yintercept = y_upper, size=0.5)+
              geom_hline(aes(yintercept=0,color=group), size=0.8)+
              color_theme+
              labs(title=gene_title, x=axis_x_title[gene_plot])+
            theme_void()+
              theme(
                panel.spacing.y = unit(x = 0, units = "line"),
                  panel.spacing.x = unit(x = 0, units = "line"),
                strip.placement = 'outside',
                axis.title.y = element_blank(),
                  axis.title.x=element_text(size=axis_x_title_size, face="bold"),
                axis.ticks = element_blank(),
                axis.text.x=element_blank(),
                axis.line = element_line(size = 0.2),
                strip.text = element_text(size=group_size, face="bold"),
                plot.title = element_text(angle=title_angle, hjust=title_hjust,size=title_size, face="bold.italic"),
                plot.margin=unit(unit(x=c(0,0,0,0), units='line')),
                  legend.position='none'
                  
              )
        }else{
            plot_res[[gene_plot]]=plot_res[[gene_plot]]+
              fill_theme+
              scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                                 expand=c(0,0))+
              geom_hline(yintercept = y_upper, size=0.5)+
              geom_hline(aes(yintercept=0,color=group), size=0.8)+
              color_theme+
              labs(title=gene_title, x=axis_x_title[gene_plot])+
            theme_void()+
              theme(
                panel.spacing.y = unit(x = 0, units = "line"),
                  panel.spacing.x = unit(x = 0, units = "line"),
                strip.placement = 'outside',
                  strip.text.y.left = element_blank(),
                axis.title.y = element_blank(),
                  axis.title.x=element_text(size=axis_x_title_size, face="bold"),
                axis.ticks = element_blank(),
                axis.text.x=element_blank(),
                axis.line = element_line(size = 0.2),
                strip.text = element_text(size=group_size, face="bold"),
                plot.title = element_text(angle=title_angle, hjust=title_hjust,size=title_size, face="bold.italic"),
                  #plot.background=element_rect(fill='red'),
                  legend.position='none',
                  plot.margin=unit(unit(x=c(0,0,0,0), units='line'))
              )
        }
    }
  wp = list()
  for(i in 1:length(plot_res)){
    wp[[i]] = plot_res[[i]]
  }
  return(wrap_plots(wp, ncol = length(plot_res)))
}