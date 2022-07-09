#' Scatter Plots of Expression vs. Time.
#'
#' @param x DCP_Rhythmicity() output
#' @param genes.plot a vector of character strings. Names of genes to be plotted. Names should be same in nomenclature as gname in the data list. If NULL, the top 10 most rhythmic genes from group I will be used
#' @param Info1 character string. Used in the plot title for group I.
#' @param Info2 character string. Used in the plot title for group II (if exist).
#' @param filename character string of the filename for exporting. If NULL, plots are not saved.
#' @param file.width height of the export plot
#' @param file.height width of the export plot
#' @param text.size the size of annotation texts.
#'
#' @return A list of ggplot2 plots.
#' @export
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#'
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' rhythm.plots = DCP_ScatterPlot(rhythm.res)
#' #to display plot in Rstudio use DCP_PlotDisplay()
#' #display the first five plots
#' DCP_PlotDisplay(rhythm.plots, id = 1:5)
DCP_ScatterPlot = function(x, genes.plot = NULL,
                          Info1 = "gI", Info2 = "gII",
                          text.size = 12,
                          filename = NULL, file.width = 8, file.height = 8){
  if(all(c("x1", "x2", "rhythm.joint")%in%names(x))){
    x1 = x$x1
    x2 = x$x2
  }else{
    x1 = x
    x2 = NULL
  }
  if(is.null(x2)){
    stopifnot("x1$data must be dataframe" = is.data.frame(x1$data))
    stopifnot("x1 does not contain result from DCP_Rhythmicity()" = !is.null(x1$rhythm))
    stopifnot("Number of samples in data does not match that in time. " = ncol(x1$data)==length(x1$time))
    stopifnot("Please input the gene labels x1$gname. " = !is.null(x1$gname))
    stopifnot("Number of gnames does not match number of genes in data. " = nrow(x1$data)==length(x1$gname))
  }else{
    stopifnot("x1$data must be dataframe and x2$data must be dataframe. " = (is.data.frame(x1$data)&is.data.frame(x2$data)))
    stopifnot("x1 or x2 do not contain result from DCP_Rhythmicity()" = !(is.null(x1$rhythm)|is.null(x2$rhythm)))
    stopifnot("Number of samples in data does not match that in time. " = (ncol(x1$data)==length(x1$time))&(ncol(x2$data)==length(x2$time)))
    stopifnot("Please input the gene labels x1$gname and x2$gname. " = !(is.null(x1$gname)|is.null(x2$gname)))
    stopifnot("Number of gnames does not match number of genes in data. " = (nrow(x1$data)==length(x1$gname))&(nrow(x2$data)==length(x2$gname)))
    gname.overlap = intersect(x1$gname, x2$gname)
    stopifnot("There are no overlapping genes between x1$gname and x2$gname. " = length(gname.overlap)>0)
  }

  if(is.null(genes.plot)){
    genes.plot = x1$gname[order(x1$rhythm$pvalue)][1:ifelse(length(x1$gname)>10, 10, length(x1$gname))]
  }
  P = x1$P
  p1.list = lapply(1:length(genes.plot), function(a){
    gene1 = genes.plot[a]
    t1 = x1$time
    expr1 = as.numeric(x1$data[match(gene1, x1$gname),])
    apar1 = x1$rhythm[match(gene1, x1$rhythm$gname),]

    myylim = c(min(expr1, apar1$M-apar1$A), max(expr1, apar1$M+apar1$A))
    # tod1=t1; period=P
    p1 = circadianDrawing_one(t1, expr1, apar1, gene1, P, specInfo1=Info1, myylim, text.size)
    return(p1)
  })

  names(p1.list) = genes.plot
  if(!is.null(x2)){

    p2.list0 =
      lapply(1:length(genes.plot), function(a){
        gene2 = genes.plot[a]
        t2 = x2$time
        expr1 = as.numeric(x1$data[match(gene2, x1$gname),])
        expr2 = as.numeric(x2$data[match(gene2, x2$gname),])
        apar1 = x1$rhythm[match(gene2, x1$rhythm$gname),]
        apar2 = x2$rhythm[match(gene2, x2$rhythm$gname),]
        myylim = c(min(expr1, expr2, apar1$M-apar1$A, apar2$M-apar2$A), max(expr1, expr2, apar1$M+apar1$A, apar2$M+apar2$A))
        p1 = suppressMessages(p1.list[[a]]+ggplot2::ylim(myylim[1], myylim[2]))
        p2 = circadianDrawing_one(t2, expr2, apar2, gene2, P, specInfo1=Info2, myylim, text.size)
        return(list(p1 = p1,
                    p2 = p2))
      })
    p1.list = lapply(p2.list0, `[[`, 1)
    p2.list = lapply(p2.list0, `[[`, 2)
    names(p1.list) = names(p2.list) = genes.plot
  }

  if(!is.null(filename)){
    grDevices::pdf(filename, file.width, file.height)
    if(is.null(x2)){
      lapply(1:length(p1.list), function(a){
        print(p1.list[[a]])
      })
    }else{
      lapply(1:length(p1.list), function(a){
        gridExtra::grid.arrange(p1.list[[a]], p2.list[[a]], ncol = 2)
      })
    }

    grDevices::dev.off()
  }

  if(is.null(x2)){
    return(p1.list)
  }else{
    return(list(x1 = p1.list,
                x2 = p2.list))
  }

}

#' Displaying DCP_ScatterPlot outputs
#'
#' @param x DCP_ScatterPlot output
#' @param id integer. Indexes for plots to display.
#'
#' @export
DCP_PlotDisplay = function(x = DCP_ScatterPlot(x, genes.plot = NULL,
                                             Info1 = "gI", Info2 = "gII",
                                             filename = NULL, height = 8, width = 8),
                           id = NULL){
  if(length(x[[1]])==1){
    return(list(gridExtra::grid.arrange(x[[1]][[1]], x[[2]][[1]], ncol = 2)))
  }else{
    if(is.null(id)){
      return(lapply(1:length(x[[1]]), function(a.g){
        gridExtra::grid.arrange(x[[1]][[a.g]], x[[2]][[a.g]], ncol = 2)
      }))
    }else{
      return(lapply(id, function(a.g){
        gridExtra::grid.arrange(x[[1]][[a.g]], x[[2]][[a.g]], ncol = 2)
      }))
    }
  }
}

circadianDrawing_one = function(tod1, expr1, apar1, gene1, period,
                                specInfo1=NULL, myylim, text.size = 12){

  geneName <- gene1
  fun.cosinor = function(x){
    apar1$M+apar1$A*cos(2*pi/period*x+apar1$phase)
  }
  amain1 <- paste0(specInfo1, " " ,geneName,
                   "\n",  " p=", round(apar1$pvalue, 4), ", ", "R2= ", round(apar1$R2, 2), ", A= ", round(apar1$A, 2), ", phase= ", round(apar1$peak, 1),
                   "\n", "sigma= ", round(apar1$sigma, 2))

  df = data.frame(Time = tod1, Expression = expr1)
  p1 = ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Expression))+
    ggplot2::geom_point(size = 2)+
    ggplot2::geom_function(fun = fun.cosinor, color = "red", size = 2)+
    ggplot2::ylim(myylim[1], myylim[2])+
    ggplot2::ggtitle(amain1)+
    # ggplot2::guides(color=ggplot2::guide_legend(title=category1.label))+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="bottom",
                   plot.title = ggplot2::element_text(size=text.size+2, hjust = 0.5))

  return(p1)
}

#' Heatmaps of Rhythmicity Signals
#'
#' @param x DCP_Rhythmicity() output
#' @param genes.plot a vector of character strings. Names of genes to be plotted. Names should be same in nomenclature as gname in the data list. If NULL, the top 100 most rhythmic genes from group I will be used
#' @param col_breaks color breaks
#' @param col_low_mid_high color choice for the heatmaps
#' @param Info1 character string. Used in the plot title for group I.
#' @param Info2 character string. Used in the plot title for group II (if exist).
#' @param ... other argument to pass to \code{ComplexHeatmap::Heatmap()}
#' @param filename character string of the filename for exporting. If NULL, plots are not saved.
#' @param file.width height of the export plot
#' @param file.height width of the export plot
#'
#' @export
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#'
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' DCP_PlotHeatmap(rhythm.res)
DCP_PlotHeatmap = function(x, genes.plot = NULL,
                          col_breaks = c(seq(-2,0,length=100),
                                         seq(.1,1,length=200),
                                         seq(1.1,2, length = 200)),
                          col_low_mid_high = c("blue","yellow","gold"),
                          Info1 = "group I", Info2 = "group II", ...,
                          filename = NULL, file.width = 8, file.height = 4){
  standardize = function(gene.i){
    x = as.numeric(gene.i)
    names(x) = names(gene.i)
    s.x = (x-mean(x))/stats::sd(x)
    s.x[which(s.x>2)] = 2
    s.x[which(s.x<(-2))] = -2
    return(s.x)
  }

  if(all(c("x1", "x2", "rhythm.joint")%in%names(x))){
    x1 = x$x1
    x2 = x$x2
  }else{
    x1 = x
    x2 = NULL
  }

  if(is.null(genes.plot)){
    genes.plot = x1$rhythm$gname[order(x1$rhythm$pvalue)][1:ifelse(nrow(x1$rhythm)>100, 100, nrow(x1$rhythm))]
    warning()
  }
  stopifnot("x1 is not given" = !is.null(x1))

  if(is.null(x2)){
    if(!all(genes.plot%in%x1$gname)){
      warning("Not all input genes are present in x1. The heatmap will not show the missing genes. ")
    }

    genes.plot = genes.plot[genes.plot%in%x1$gname]
    rhythm.x1 = x1$rhythm[match(genes.plot, x1$rhythm$gname), ]
    mat.x1 = x1$data[match(genes.plot, x1$gname), ]
    rownames(mat.x1) = genes.plot
    mat.x1 = mat.x1[order(rhythm.x1$peak), order(x1$time)]
    mat.x1 = t(apply(mat.x1, 1, standardize))

    p1 = ComplexHeatmap::Heatmap(mat.x1, name = Info1, col = circlize::colorRamp2(col_breaks, grDevices::colorRampPalette(col_low_mid_high)(n = 500)),
                                 cluster_rows = FALSE, cluster_columns = FALSE, ...)
    if(!is.null(filename)){
      grDevices::pdf(filename, width = file.width, height = file.height)
      print(p1)
      grDevices::dev.off()
    }

    return(p1)
  }else{
    if(!(all(genes.plot%in%x1$gname)&all(genes.plot%in%x2$gname))){
      warning("Not all input genes are present in x1 or x2. The heatmap will not show genes in neither x1 nor x2. The side by side heatmaps might contain missingness.")
    }
    genes.plot = genes.plot[genes.plot%in%x1$gname|genes.plot%in%x2$gname]
    rhythm.x1 = x1$rhythm[match(genes.plot, x1$rhythm$gname), ]
    mat.x1 = x1$data[match(genes.plot, x1$gname), ]
    rownames(mat.x1) = genes.plot
    mat.x1 = mat.x1[order(rhythm.x1$peak), order(x1$time)]
    mat.x1 = t(apply(mat.x1, 1, standardize))

    rhythm.x2 = x2$rhythm[match(genes.plot, x2$rhythm$gname), ]
    mat.x2 = x2$data[match(genes.plot, x2$gname), ]
    rownames(mat.x2) = genes.plot
    mat.x2 = mat.x2[order(rhythm.x1$peak), order(x2$time)]
    mat.x2 = t(apply(mat.x2, 1, standardize))

    p1 = ComplexHeatmap::Heatmap(mat.x1, name = Info1, col = circlize::colorRamp2(col_breaks, grDevices::colorRampPalette(col_low_mid_high)(n = 500)),
                                 cluster_rows = FALSE, cluster_columns = FALSE, ...)
    p2 = ComplexHeatmap::Heatmap(mat.x2, name = Info2, col = circlize::colorRamp2(col_breaks, grDevices::colorRampPalette(col_low_mid_high)(n = 500)),
                                 cluster_rows = FALSE, cluster_columns = FALSE, ...)
    if(!is.null(filename)){
      grDevices::pdf(filename, width = file.width, height = file.height)
      print(p1+p2)
      grDevices::dev.off()
    }
    return(p1+p2)
  }


  # grab_grob <- function(){
  #   gridGraphics::grid.echo()
  #   grid::grid.grab()
  # }
  # heatmap.plots = lapply(list(mat.x1, mat.x2), function(a){
  #   gplots::heatmap.2(mat.x1, Rowv = NA, Colv = NA, main = paste0("\n\n", Info1),
  #                     ylab = info.gene,  col = colorRampPalette(c("blue","yellow","gold"))(n = 499), symkey=FALSE, trace = "none",
  #                     cexRow = 0.5, key = TRUE,density.info = "none", breaks = col_breaks, dendrogram = "none", margins = c(6, 7),
  #                     lmat = rbind(c(4, 4, 3, 3),c(2, 1, 1, 1)),lhei = c(1, 4), lwid = c(0.5,0.5, 1,2))
  #   grab_grob()
  # })

}




#' Histogram of peak time
#' Make a circos histogram plot for peak time.
#'
#' @param x DCP_Rhythmicity() output
#' @param TOJR toTOJR output. If NULL, rhythm.joint object in x will be used.
#' @param RhyBothOnly For two-group output, plot only RhyBoth genes or all rhythmic genes in separate groups?
#' @param sig.cut A list. Used only for single-group plot, only genes satisfying sig.cut will be plotted. If NULL then genes all genes in regardless of significance in rhythmicity will be plotted. \itemize{
#' \item param parameter used for the cutoff. Should be a column in rhythmicity estimates result (e.g. x$rhythm for one-group analysis, or x$x1$rhythm for two-group analysis).
#' \item fun character string. Either "<" or ">"
#' \item val numeric. The value used for the cutoff}
#' @param time.start numeric. What time do you want the phase start? Default is -6, which is midnight if time is in ZT scale.
#' @param Info1 character string. Used in the plot title for group I
#' @param Info2 character string. Used in the plot title for group II (if exist).
#' @param GroupSplit logical. For two-group output, should the histogram of each group be plotted separately?
#' @param filename character string. The filename for plot export. If NULL, the plot will not be saved.
#' @param file.width width of the export plot
#' @param file.height height of the export plot
#' @param color.hist Input color. The length of the vector should be the same of the number of the groups.
#' @param cir.y.breaks numeric. A vector for breaks for the angles. Should start with phase.start and end with phase.start+period
#' @param single.binwidth numeric. The binwidth for plotting peak histogram.
#' @param axis.text.size Size for the axis text.
#' @param legend.position One of "left”, "top", "right", "bottom", or "none"
#'
#' @return
#' @export
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(2, 3), A2=c(2, 3),
#' phase1=c(0, pi/4), phase2=c(pi/2, pi*3/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' # Make a two-group plot:
#' DCP_PlotPeakHist(rhythm.res)
#' # Make a one-group plot
#' DCP_PlotPeakHist(rhythm.res$x1)
DCP_PlotPeakHist = function(x, TOJR = NULL, RhyBothOnly = FALSE, sig.cut = list(param = "pvalue", fun = "<", val = 0.05),
                           time.start = -6,
                           Info1 = "groupI", Info2 = "groupII", GroupSplit = FALSE,
                           filename = NULL, file.width = 8, file.height = 8,
                           color.hist = NULL,
                           cir.y.breaks = seq(-6,18, 4),
                           single.binwidth = 1,
                           axis.text.size = 12,
                           legend.position="right"){

  studyType = To.studyType(x)

  if(studyType == "One"){

    period = x$P
    a.min = time.start
    a.max = time.start+period
    if(is.null(color.hist)){
      color.hist = "#3374b0"
    }

    if(is.null(sig.cut)){

      warning("sig.cut input is NULL. All genes will be plotted. ")
      peak.df = x$rhythm

    }else{

      peak.df = x$rhythm
      CheckSigCut(peak.df, sig.cut, "x$rhythm")
      xx = x$rhythm[, sig.cut$param]
      peak.df = x$rhythm[.Primitive(sig.cut$fun)(xx, sig.cut$val), ]

    }

    peak.df$peak = adjust.circle(peak.df$peak, time.start, period)
    pp.file = paste0(filename, "_", Info1, "_PeakHist")

    pp = ggplot2::ggplot(data = peak.df,ggplot2::aes(y = peak))+
      ggplot2::geom_histogram(binwidth = single.binwidth, fill = color.hist) +
      ggplot2::scale_y_continuous(breaks = cir.y.breaks, limits = c(a.min,a.max),
                                  labels = cir.y.breaks) +
      ggplot2::xlab("") + ggplot2::ylab(paste0("Phase in ", Info1)) +
      ggplot2::ggtitle(paste0("Histogram of ohase in", Info1))+
      ggplot2::theme_bw()+
      ggplot2::theme(aspect.ratio = 1, axis.line = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = axis.text.size),
                     axis.text.y = ggplot2::element_text(size = axis.text.size),
                     panel.border = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size=axis.text.size*0.8),
                     legend.text = ggplot2::element_text(size=axis.text.size*0.8),
                     plot.title = ggplot2::element_text(size=axis.text.size+2, hjust = 0.5)) +
      ggplot2::coord_polar(theta="y", start=0)
  }else{

    stopifnot("x$x1$P is not equal to x$x2$P. " = x$x1$P==x$x2$P)
    period = x$x1$P
    a.min = time.start
    a.max = time.start+period

    if(is.null(color.hist)){
      color.hist = c("#F8766D", "#00BFC4")
    }

    if(is.null(TOJR)){
      stopifnot("There is no TOJR input, and x$rhythm.joint$TOJR is also NULL." = !is.null(x$rhythm.joint$TOJR))
      warning("There is no TOJR input, x$rhythm.joint$TOJR will be used for joint rhythmicity category. ")
      TOJR = x$rhythm.joint$TOJR
      rhyI = x$rhythm.joint$gname[TOJR == "rhyI"];
      rhyII = x$rhythm.joint$gname[TOJR == "rhyII"];
      rhyboth = x$rhythm.joint$gname[TOJR == "both"];
      rhyI.both = c(rhyI, rhyboth);
      rhyII.both = c(rhyII, rhyboth);

    }else{
      warning("TOJR input will be used joint rhythmicity category. ")
      rhyI = x$gname_overlap[TOJR == "rhyI"]
      rhyII = x$gname_overlap[TOJR == "rhyII"]
      rhyboth = x$gname_overlap[TOJR == "both"]
      rhyI.both = c(rhyI, rhyboth)
      rhyII.both = c(rhyII, rhyboth)
    }

    if(RhyBothOnly){
      peak.df.long = data.frame(peak = c(peak.select(x, rhyboth, "x1"), peak.select(x, rhyboth, "x2")),
                                group = factor(rep(c(Info1, Info2), each = length(rhyboth))))
    }else{
      peak.df.long = data.frame(peak = c(peak.select(x, rhyI.both, "x1"), peak.select(x, rhyII.both, "x2")),
                                group = factor(c(rep(Info1, length(rhyI.both)),
                                                 rep(Info2, length(rhyII.both))
                                )))
    }

    peak.df.long$peak = adjust.circle(peak.df.long$peak, time.start, period)

    pp.file = paste0(filename, "_", Info1, "_", Info2, "_PeakHist")
    pp = ggplot2::ggplot(data = peak.df.long, ggplot2::aes(y = peak, fill = group))+
      #binwidth = single.binwidth, fill = "#3374b0"
      ggplot2::geom_histogram(binwidth = single.binwidth, position = 'identity', alpha = 0.6)+
      ggplot2::xlab(paste0("")) + ggplot2::ylab("") +
      ggplot2::ggtitle(paste0("Histogram of peak time"))+
      ggplot2::scale_y_continuous(breaks = cir.y.breaks, limits = c(a.min,a.max), labels = cir.y.breaks) +
      # ggplot2::scale_y_continuous(breaks = cir.y.breaks, labels = cir.y.breaks) +
      ggplot2::scale_fill_manual(
        values = color.hist,
        aesthetics = "fill")+
      ggplot2::theme_bw()+
      ggplot2::theme(aspect.ratio = 1, axis.line = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = axis.text.size),
                     axis.text.y = ggplot2::element_text(size = axis.text.size),
                     panel.border = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size=axis.text.size*0.8),
                     legend.text = ggplot2::element_text(size=axis.text.size*0.8),
                     legend.position = legend.position,
                     plot.title = ggplot2::element_text(size=axis.text.size+2, hjust = 0.5)) +
      ggplot2::coord_polar(theta="y", start=0)+
      {if(GroupSplit)ggplot2::facet_wrap(~group)}+
      {if(GroupSplit)ggplot2::theme(legend.position="none")}
  }
  if(!is.null(filename)){
    grDevices::pdf(paste0(pp.file, ".pdf"), width = file.width, height = file.height)
    print(pp)
    grDevices::dev.off()
  }else{
    print(pp)
  }
  return(pp)
}

#' Radar plot of peak time
#' Make a circos histogram plot for peak time.
#'
#' @param x DCP_Rhythmicity() output
#' @param TOJR toTOJR output. If NULL, rhythm.joint object in x will be used.
#' @param RhyBothOnly For two-group output, plot only RhyBoth genes or all rhythmic genes in separate groups?
#' @param sig.cut A list. Used only for single-group plot, only genes satisfying sig.cut will be plotted. If NULL then genes all genes in regardless of significance in rhythmicity will be plotted. \itemize{
#' \item param parameter used for the cutoff. Should be a column in x.
#' \item fun character string. Either "<", or ">"
#' \item val numeric. The value used for the cutoff}
#' @param time.start numeric. What time do you want the phase start? Default is -6, which is midnight if time is in ZT scale.
#' @param Info1 character string. Used in the plot title for group I
#' @param Info2 character string. Used in the plot title for group II (if exist).
#' @param filename character string. The filename for plot export. If NULL, the plot will not be saved.
#' @param file.width width of the export plot
#' @param file.height height of the export plot
#' @param color Input color. The length of the vector should be the same of the number of the groups.
#' @param single.binwidth numeric. The binwidth for plotting peak histogram.
#' @param axis.text.size Size for the axis text.
#' @param legend.position One of "left”, "top", "right", "bottom", or "none"
#'
#' @return
#' @export
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(2, 3), A2=c(2, 3),
#' phase1=c(0, pi/4), phase2=c(pi/2, pi*3/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' # Make a two-group plot:
#' DCP_PlotPeakRadar(rhythm.res)
#' # Make a one-group plot
#' DCP_PlotPeakRadar(rhythm.res$x1)
DCP_PlotPeakRadar = function(x, TOJR = NULL, RhyBothOnly = FALSE, sig.cut = list(param = "pvalue", fun = "<", val = 0.05),
                             time.start = -6,
                             Info1 = "groupI", Info2 = "groupII",
                             filename = NULL, file.width = 8, file.height = 8,
                             color = NULL,
                             single.binwidth = 2,
                             axis.text.size = 12,
                             legend.position="right"){

  to.df.radar = function(peak.vec = peak.df$peak, group = "gI", bin.width = 1, phase.start, period){
    a.break = seq(phase.start, phase.start+period, by = bin.width)
    # ranges = cut(peak.vec, breaks = a.break)
    # a.df = as.data.frame(table(ranges))

    ranges = lapply(peak.vec, function(a){cut(a, breaks = a.break)})
    a.df = lapply(ranges, function(a){as.data.frame(table(a))})
    df.grid = a.df[[1]][, 1]
    a.df = cbind.data.frame(lapply(a.df, `[[`, 2))
    # a.df.grid1 = a.break[-length(a.break)]
    # a.df.grid2 = a.break[-1]
    # a.df.grid = (a.df.grid1+a.df.grid2)/2
    # a.df$ranges = a.df.grid
    a.max = max(a.df)
    a.df2 = as.data.frame(t(a.df)); colnames(a.df2) = df.grid; rownames(a.df2) = NULL
    df.radar = a.df2/a.max
    df.radar = cbind.data.frame(data.frame(group = group),
                                df.radar)
    return(list(df = df.radar,
                max = a.max))
  }

  studyType = To.studyType(x)

  if(studyType == "One"){

    period = x$P
    a.min = time.start
    a.max = time.start+period
    if(is.null(color)){
      color = "#3374b0"
    }

    if(is.null(sig.cut)){

      warning("sig.cut input is NULL. All genes will be plotted. ")
      peak.df = x$rhythm

    }else{

      peak.df = x$rhythm
      CheckSigCut(peak.df, sig.cut, "x$rhythm")
      xx = x$rhythm[, sig.cut$param]
      peak.df = x$rhythm[.Primitive(sig.cut$fun)(xx, sig.cut$val), ]

    }

    peak.df$peak = adjust.circle(peak.df$peak, time.start, period)
    pp.file = paste0(filename, "_", Info1, "_PeakRadar")

    df.radar = to.df.radar(peak.vec = list(peak.df$peak), group = "gI", bin.width = single.binwidth, a.min, period)

    pp =     suppressMessages(
      ggradar::ggradar(df.radar$df,
                       values.radar = c(0, round(max(df.radar$max)/2), max(df.radar$max)),
                       # font.radar = "roboto",
                       grid.label.size = axis.text.size*0.8,  # Affects the grid annotations (0%, 50%, etc.)
                       axis.label.size = axis.text.size*0.5, # Afftects the names of the variables
                       group.point.size = axis.text.size*0.1   # Simply the size of the point
      )+
        ggplot2::scale_color_manual(#name = "color",
          # breaks = ,
          values = color,
          # labels = sig.color.breaks
        )+
        ggplot2::theme(
          legend.position = "none"
          # legend.justification = c(1, 0),
          # legend.text = ggplot2::element_text(size = 14),
          # legend.key = ggplot2::element_rect(fill = NA, color = NA),
          # legend.background = ggplot2::element_blank()
        ) +
        ggplot2::labs(title = paste0("Radar plot of peak time")) +
        ggplot2::theme(
          # plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
          # panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
          # # plot.background = element_rect(fill = "grey", color = "#fbf9f4"),
          # panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
          # plot.title.position = "plot", # slightly different from default
          plot.title = ggplot2::element_text(
            # family = "lobstertwo",
            size = axis.text.size*2.5,
            face = "bold",
            color = "#2a475e",
            hjust=0.5
          )
        )
    )


  }else{

    stopifnot("x$x1$P is not equal to x$x2$P. " = x$x1$P==x$x2$P)
    period = x$x1$P
    a.min = time.start
    a.max = time.start+period

    if(is.null(color)){
      color = c("#F8766D", "#00BFC4")
    }

    if(is.null(TOJR)){
      stopifnot("There is no TOJR input, and x$rhythm.joint$TOJR is also NULL." = !is.null(x$rhythm.joint$TOJR))
      warning("There is no TOJR input, x$rhythm.joint$TOJR will be used for joint rhythmicity category. ")
      TOJR = x$rhythm.joint$TOJR
      rhyI = x$rhythm.joint$gname[TOJR == "rhyI"];
      rhyII = x$rhythm.joint$gname[TOJR == "rhyII"];
      rhyboth = x$rhythm.joint$gname[TOJR == "both"];
      rhyI.both = c(rhyI, rhyboth);
      rhyII.both = c(rhyII, rhyboth);
    }else{
      warning("TOJR input will be used joint rhythmicity category. ")
      rhyI = x$gname_overlap[TOJR == "rhyI"]
      rhyII = x$gname_overlap[TOJR == "rhyII"]
      rhyboth = x$gname_overlap[TOJR == "both"]
      rhyI.both = c(rhyI, rhyboth)
      rhyII.both = c(rhyII, rhyboth)
    }

    if(RhyBothOnly){
      peak.df.long = data.frame(peak = c(peak.select(x, rhyboth, "x1"), peak.select(x, rhyboth, "x2")),
                                group = factor(rep(c(Info1, Info2), each = length(rhyboth))),
                                group2 = factor(rep(c(1, 2), each = length(rhyboth))))
    }else{
      peak.df.long = data.frame(peak = c(peak.select(x, rhyI.both, "x1"), peak.select(x, rhyII.both, "x2")),
                                group = factor(c(rep(Info1, length(rhyI.both)),
                                                 rep(Info2, length(rhyII.both))
                                )),
                                group2 = factor(c(rep(1, length(rhyI.both)),
                                                 rep(2, length(rhyII.both))
                                )))
    }

    peak.df.long$peak = adjust.circle(peak.df.long$peak, time.start, period)


    pp.file = paste0(filename, "_", Info1, "_", Info2, "_PeakHist")
    df.radar = to.df.radar(peak.vec = list(peak.df.long[peak.df.long$group2==1, "peak"],
                                           peak.df.long[peak.df.long$group2==2, "peak"]),
                           group = c(Info1, Info2), bin.width = single.binwidth, a.min, period)

    pp =     suppressMessages(
      ggradar::ggradar(df.radar$df,
                       values.radar = c(0, round(max(df.radar$max)/2), max(df.radar$max)),
                       # font.radar = "roboto",
                       grid.label.size = axis.text.size*0.8,  # Affects the grid annotations (0%, 50%, etc.)
                       axis.label.size = axis.text.size*0.5, # Afftects the names of the variables
                       group.point.size = axis.text.size*0.1   # Simply the size of the point
      )+
        ggplot2::scale_color_manual(#name = "color",
          # breaks = ,
          values = color,
          # labels = sig.color.breaks
        )+
        ggplot2::theme(
          legend.position = legend.position,
          legend.justification = c(1, 0),
          legend.text = ggplot2::element_text(size = axis.text.size),
          legend.key = ggplot2::element_rect(fill = NA, color = NA),
          legend.background = ggplot2::element_blank()
        ) +
        ggplot2::labs(title = paste0("Radar plot of peak time")) +
        ggplot2::theme(
          # plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
          # panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
          # # plot.background = element_rect(fill = "grey", color = "#fbf9f4"),
          # panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
          # plot.title.position = "plot", # slightly different from default
          plot.title = ggplot2::element_text(
            # family = "lobstertwo",
            size = axis.text.size*2.5,
            face = "bold",
            color = "#2a475e",
            hjust=0.5
          )
        )
    )
  }
  if(!is.null(filename)){
    grDevices::pdf(paste0(pp.file, ".pdf"), width = file.width, height = file.height)
    print(pp)
    grDevices::dev.off()
  }else{
    print(pp)
  }
  return(pp)

}


#' Circos plot of peak time shift
#' Make a circos plot for peak time shift between two groups
#'
#' @param x DCP_Rhythmicity() output.
#' @param TOJR toTOJR output. If NULL, rhythm.joint object in x will be used.
#' @param dPhase DCP_DiffPar() output with phase shift tested.
#' @param color.cut list. Genes will be plotted according a cutoff. Used only when dPhase is not NULL. \itemize{
#' \item param parameter used for the cutoff. Should be a column in dPhase.
#' \item fun character string. Either "<", ">", or "="
#' \item val numeric. The value used for the cutoff.
#' \item color.sig color for points that pass the color.cut criterion.
#' \item color.none color for points that do not pass the color.cut criterion.}
#' @param color.df data.frame. A more advanced color coding. The data frame should have two columns named color and label. label will be displayed as legend. The color vector should have the same order as genes in dPhase.
#' @param time.start numeric. What time do you want the phase start? Default is -6, which is midnight if time is in ZT scale.
#' @param Info1 character string. Used in the plot title for group I
#' @param Info2 character string. Used in the plot title for group II (if exist).
#' @param filename character string. The filename for plot export. If NULL, the plot will not be saved.
#' @param file.width width of the export plot
#' @param file.height height of the export plot
#' @param concordance.ref The radius where the concordance reference line be plotted away from \eqn{\Delta}peak = 0.
#' @param cir.x.breaks numeric. A vector for breaks for the radius (phase difference). Should only contains values from -period/2 to period/2 and it is recommended that the break is equal spaced.
#' @param cir.y.breaks numeric. A vector for breaks for the angles. Should start with time.start and end with time.start+period
#' @param axis.text.size numeric. Size for the axis text.
#' @param legend.position One of "left”, "top", "right", "bottom", or "none"
#' @param color.diff.refband color of the reference band around \eqn{\Delta}peak = 0.
#' @param color.diff.xlim color of the start and end of the phase difference range.
#' @param color.diff.baseline color of the reference line for \eqn{\Delta}peak = 0.
#'
#' @return
#' @export
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(2, 3), A2=c(2, 3),
#' phase1=c(0, pi/4), phase2=c(pi/2, pi*3/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' rhythm.diffPar = DCP_DiffPar(rhythm.res, "phase")
#' #make a plot with genes with differential phase p-value<0.05 in a different color
#' DCP_PlotPeakDiff(rhythm.res, NULL, rhythm.diffPar)
DCP_PlotPeakDiff = function(x, TOJR = NULL, dPhase = NULL,
                            color.cut = list(param = "pvalue", fun = "<", val = 0.05, color.sig = "#b33515", color.none = "dark grey"),
                            color.df = NULL,
                            Info1 = "groupI", Info2 = "groupII",
                            filename = NULL, file.width = 8, file.height = 8,
                            time.start = -6, concordance.ref = 4,
                            cir.x.breaks = seq(-12,12, 4), cir.y.breaks = seq(-6,18, 4),
                            axis.text.size = 12, legend.position="right",
                            color.diff.refband = "darkgreen",
                            color.diff.xlim = "grey40",
                            color.diff.baseline = "blue"){

  studyType = To.studyType(x)
  stopifnot("DCP_PlotPeakDiff can only be plotted for DCP_Rhythmicity() of two-group analysis. " = studyType == "Two")

  stopifnot("x$x1$P is not equal to x$x2$P. " = x$x1$P==x$x2$P)
  period = x$x1$P

  a.min = time.start
  a.max = time.start+period


  if(is.null(dPhase)){
    if((!is.null(color.cut))|(!is.null(color.df))){
      warning("color.cut or color.df input is ignored. color.cut or color.df is only used when dPhase is given. ")
      color.cut = NULL; color.df = NULL
    }
    if(is.null(TOJR)){

      stopifnot("There is no input for TOJR or dPhase, and x$rhythm.joint$TOJR is also NULL." = !is.null(x$rhythm.joint$TOJR))
      TOJR = x$rhythm.joint$TOJR
      warning("There is no input for TOJR or dPhase, x$rhythm.joint$TOJR will be used for extracting RhyBoth genes. ")
      rhyboth = x$rhythm.joint$gname[TOJR == "both"];
      print(1)

    }else{

      warning("dPhase is not given. TOJR will be used to identify RhyBoth genes. ")
      rhyboth = x$gname_overlap[TOJR == "both"]

    }
  }else{
    if((!is.null(color.cut))&(!is.null(color.df))){
      warning("Both color.cut and color.df are non-NULL, color.df will be used. ")
      color.cut = NULL
    }
    warning("dPhase is given. Genes contained in dPhase will be used as RhyBoth in regardless of any TOJR input.")
    rhyboth = dPhase$gname

  }

  peak.df = data.frame(peak1 = peak.select(x, rhyboth, "x1"),
                       peak2 = peak.select(x, rhyboth, "x2"))
  peak.df$peak1 = adjust.circle(peak.df$peak1, time.start, period)
  peak.df$peak2 = adjust.circle(peak.df$peak2, time.start, period)
  peak.df$delta.peak = peak.df$peak2 - peak.df$peak1
  peak.df$delta.peak[peak.df$delta.peak>(period/2)] = peak.df$delta.peak[peak.df$delta.peak>(period/2)] -period
  peak.df$delta.peak[peak.df$delta.peak< -(period/2)] = peak.df$delta.peak[peak.df$delta.peak< -(period/2)] +period
  pp.file = paste0(filename, "_", Info1, "_", Info2, "_PeakDiff")

  if(!is.null(dPhase)){
    if(!is.null(color.cut)){
      CheckSigCut(dPhase, color.cut, "dPhase")
      xx = ifelse(rep(color.cut$param == "delta.peak", nrow(dPhase)), abs(dPhase[, color.cut$param]), dPhase[, color.cut$param])
      sig.color = .Primitive(color.cut$fun)(xx, color.cut$val)
      sig.color = ifelse(sig.color, "sig", "none") #red: #b33515
      legend.label = ifelse(color.cut$param == "delta.peak", paste0("|phase difference|", color.cut$fun, color.cut$val), c(paste(unlist(color.cut[1:3]), collapse="")))
      #paste(unlist(color.cut), collapse="")
    }else if(!is.null(color.df)){
      sig.color = color.df$label
      sig.color.breaks = names(table(color.df$label))
      sig.color.values.ind = apply(table(color.df$color, color.df$label), 2, function(a){which(a!=0)})
      if(class(sig.color.values.ind)=="list"){
        stop("Please make sure that label and color in color.df are one-to-one matched. ")
      }
      sig.color.values = rownames(table(color.df$color, color.df$label))[sig.color.values.ind]
    }else{
      sig.color = NULL
    }
  }else{
    sig.color = NULL
  }


  peak.df$delta.peak2 = ifelse(peak.df$delta.peak>=cir.x.breaks[1], peak.df$delta.peak, peak.df$delta.peak+period)
  cir.x.breaks2 = ifelse(cir.x.breaks>=cir.x.breaks[1], cir.x.breaks, cir.x.breaks+period)
  if(cir.x.breaks[1]==utils::tail(cir.x.breaks, 1)){cir.x.breaks2[length(cir.x.breaks2)] = cir.x.breaks2[length(cir.x.breaks2)] + period}
  highlight.radius = period/6
  if(sum(cir.x.breaks2%%period==0)){#if 0 is in the given axis
    highlight.center =  cir.x.breaks2[which(cir.x.breaks2%%period==0)]
  }else{
    cir.x.breaks2.resid = cir.x.breaks2-period
    cir.x.breaks2.resid.min.ind = which.min(abs(cir.x.breaks2.resid))
    highlight.center = cir.x.breaks2[cir.x.breaks2.resid.min.ind]-cir.x.breaks2.resid[cir.x.breaks2.resid.min.ind]
  }

  InRange = function(x, y){
    x.min = min(x); x.max = max(x)
    if(y<=x.max&y>=x.min){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }

  if(InRange(cir.x.breaks2, highlight.center+highlight.radius)&InRange(cir.x.breaks2, highlight.center-highlight.radius)){
    rects = data.frame(start = highlight.center-highlight.radius, end = highlight.center+highlight.radius, group = 1)
  }else if(InRange(cir.x.breaks2, highlight.center+highlight.radius)){
    rects1 = data.frame(start = min(cir.x.breaks), end = highlight.center+highlight.radius, group = 1)
    rects2 = data.frame(start = highlight.center-highlight.radius+period, end = max(cir.x.breaks), group = 2)
    rects = rbind.data.frame(rects1, rects2)
  }else{
    InRange(cir.x.breaks2, highlight.center-highlight.radius)
    rects1 = data.frame(start = highlight.center-highlight.radius, end = max(cir.x.breaks), group = 1)
    rects2 = data.frame(start = min(cir.x.breaks), end = highlight.center+highlight.radius-period, group = 2)
    rects = rbind.data.frame(rects1, rects2)
  }


  reference.unit = abs(cir.y.breaks[1])+abs(utils::tail(cir.y.breaks, 1))
  pp = ggplot2::ggplot(data = peak.df,ggplot2::aes(x = delta.peak2, y = peak1))+
    ggplot2::geom_rect(data=rects, inherit.aes=FALSE, ggplot2::aes(xmin=start, xmax=end, ymin=min(cir.y.breaks),
                                                                   ymax=max(cir.y.breaks), group=group), color="transparent", fill=color.diff.refband, alpha=0.1)+
    ggplot2::geom_vline(xintercept = c(cir.x.breaks2[1], utils::tail(cir.x.breaks2, 1)), linetype="dashed", color=color.diff.xlim) +
    # ggplot2::geom_vline(xintercept = c(-1*concordance.ref,concordance.ref), linetype="dashed", color="darkgreen",size=1, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = highlight.center,  color = color.diff.baseline, alpha = 0.6) +
    # ggplot2::geom_hline(yintercept = cir.y.grid, color="grey40") +
    ggplot2::geom_point(size=0.4, ggplot2::aes(color = sig.color), alpha = 0.8) +
    {if(!is.null(color.df)) ggplot2::scale_color_manual(name = "color",
                                                        breaks = sig.color.breaks,
                                                        values = sig.color.values,
                                                        labels = sig.color.breaks)}+
    {if(!is.null(color.cut)) ggplot2::scale_color_manual(name = "color",
                                                         breaks = c("sig"),
                                                         values = c("sig" = color.cut$color.sig, "none"= color.cut$color.none),
                                                         labels=c(paste(unlist(color.cut[1:3]), collapse="")))}+
    ggplot2::geom_text(data=data.frame(x=cir.x.breaks2, y=sum(cir.y.breaks[1:2])/2, label=cir.x.breaks),
                       ggplot2::aes(x=x, y=y, label = label), nudge_x = -0.2, vjust = -1, size=axis.text.size*1/3) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    # ggplot2::ylab(paste0("Angles: peak time in ", Info1, "\n",
    #                      "Radius: ", "peak difference (", Info2, "-", Info1, ")"))+
    ggplot2::ggtitle(paste0("Phase difference (", Info2, " - ", Info1, ")"))+
    ggplot2::scale_x_continuous(breaks = cir.x.breaks2, limits = c(cir.x.breaks2[1]-1.5, utils::tail(cir.x.breaks2, 1)+2),
                                expand = c(0, -2)) +
    ggplot2::scale_y_continuous(breaks = cir.y.breaks, limits = c(a.min,a.max),
                                labels = cir.y.breaks) +
    ggplot2::theme_bw()+
    ggplot2::theme(aspect.ratio = 1, axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = axis.text.size),
                   axis.text.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size=axis.text.size*0.8),
                   legend.text = ggplot2::element_text(size=axis.text.size*0.8),
                   plot.title = ggplot2::element_text(size=axis.text.size+2, hjust = 0.5)) +
    ggplot2::annotate("segment", x = utils::tail(cir.x.breaks2, 1)+2, xend = utils::tail(cir.x.breaks2, 1)+2,
                      y = mean(cir.y.breaks)-1.5*reference.unit/24, yend = mean(cir.y.breaks)-0.5*reference.unit/24, arrow = ggplot2::arrow(length = ggplot2::unit(axis.text.size, "pt")),
                      color = color.diff.xlim)+ #the angle annotation
    ggplot2::geom_text(data = data.frame(x = utils::tail(cir.x.breaks2, 1)+2,
                                         y = mean(cir.y.breaks)-1.5*reference.unit/24),
                       ggplot2::aes(x = x, y = y),
                       vjust = 3.5, #hjust = 0.5,
                       nudge_y = -0.5,
                       label = as.expression(bquote(phi~"in"~.(Info1))))+
    # paste0("peak in ", Info1)
    ggplot2::annotate("segment", x = cir.x.breaks2[1], xend = utils::tail(cir.x.breaks2, 1)+2,
                      y = sum(cir.y.breaks[1:2])/2, yend = sum(cir.y.breaks[1:2])/2,
                      arrow = ggplot2::arrow(length = ggplot2::unit(axis.text.size, "pt")), color = color.diff.xlim)+ # the radius annotation
    ggplot2::geom_point(data = data.frame(x = cir.x.breaks2, y = sum(cir.y.breaks[1:2])/2),
                        ggplot2::aes(x = x, y = y),
                        color = color.diff.xlim, size = 1, shape = 20) + #to create tick like marks
    ggplot2::geom_text(data = data.frame(x = utils::tail(cir.x.breaks2, 1)+2,
                                         y = sum(cir.y.breaks[1:2])/2),
                       ggplot2::aes(x = x, y = y),
                       vjust = -1,
                       label = expression(Delta~phi))+
    ggplot2::coord_polar(theta="y", start=0, clip = "off")



  if(!is.null(filename)){
    grDevices::pdf(paste0(pp.file, ".pdf"), width = file.width, height = file.height)
    print(pp)
    grDevices::dev.off()
  }else{
    print(pp)
  }
  return(pp)
}



# Functions for phase plots -----------------------------------------------
adjust.circle = function(x, a.min = -6,  period = 24){
  a.max = a.min + period
  x[x<a.min] = x[x<a.min]+period
  x[x>a.max] = x[x>a.max]-period
  return(x)
}

To.studyType = function(x){
  if(all(c("x1", "x2", "gname_overlap", "rhythm.joint")%in%names(x))){
    type = "Two"
    stopifnot("x should be output of DCP_Rhythmicity" = (!is.null(x$x1$rhythm))&(!is.null(x$x2$rhythm)))
  }else if(all(c("rhythm", "P")%in%names(x))){
    type = "One"
  }else{
    stop("x should be output of DCP_Rhythmicity")
  }
  return(type)
}

CheckSigCut = function(peak.df, sig.cut, a.message){

  isColor <- function(x){
    res <- try(grDevices::col2rgb(x),silent=TRUE)
    return(!"try-error"%in%class(res))
  } #from stackflow: https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation

  a.message2 = paste0("sig.cut$param should be a column of the ", a.message)
  if(!sig.cut$param%in%colnames(peak.df)){
    stop(a.message2)
  }
  stopifnot("color.cut$fun should be one of '<', '>', '='. " = sig.cut$fun %in%c("<", ">", "="))
  stopifnot("color.cut$fun of '<' or '>' should only ve used with a numeric color.cut$val. " = sig.cut$fun %in%c("<", ">")&(is.numeric(sig.cut$val)))
  stopifnot("color.cut$color.sig should be valid color specifications. see grDevices::col2rgb(). " = isColor(sig.cut$color.sig))
  stopifnot("color.cut$color.none should be valid color specifications. see grDevices::col2rgb(). " = isColor(sig.cut$color.none))

}

peak.select = function(x, select.gnames, group){
  return(x[[group]]$rhythm$peak[x[[group]]$rhythm$gname%in%select.gnames])
}
