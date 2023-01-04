#' Differential Rhythmicity Parameter Test
#'
#' This function performs differential rhythmicity parameter tests.
#' @param x outpuf of DCP_Rhythmicity(x1, x2)
#' @param Par One of "A", "phase", "M", "A&phase", or "A&phase&M"
#' @param TOJR toTOJR output. If NULL, rhythm.joint object in x will be used.
#' @param alpha cutoff of p-values for significant differential parameter shift.
#' @param p.adjust.method input for p.adjust() in R package \code{stat}
#' @param parallel.ncores integer. Number of cores used if using parallel
#' computing with \code{mclapply()}. Not functional for windows system.
#'
#' @return A dataframe of differential rhythm parameter test results.
#' @export
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' rhythm.diffPar = DCP_DiffPar(rhythm.res, Par = "A&phase")
DCP_DiffPar = function(x, Par = c("A"), TOJR=NULL, alpha = 0.05,
                       p.adjust.method = "BH", parallel.ncores = 1){
  stopifnot('Par should be one of "A", "phase", "M", "A&phase", or "A&phase&M".'
            = Par %in% c("A", "phase", "M", "A&phase", "A&phase&M"))
  stopifnot(
    "x should be output of DCP_Rhythmicity() with both x1 and x2 not NULL" =
      all(c("x1", "x2", "gname_overlap", "rhythm.joint")%in%names(x)))

  if(is.null(TOJR)){
    overlap.g = x$rhythm.joint$gname[x$rhythm.joint$TOJR == "both"]
    stopifnot('There is no RhyBoth genes, please use a lower threshold for TOJR. ' = length(overlap.g)>0)
  }else{
    stopifnot('The input number of types of joint rhythmicity does not
              match that of overlapping genes in two groups ' =
                length(x$rhythm.joint$gname)==length(TOJR)
    )
    overlap.g = x$gname_overlap[TOJR == "both"]
  }

  x1 = x$x1
  x2 = x$x2
  x1.overlap = x1$data[match(overlap.g, x1$gname), ]
  x2.overlap = x2$data[match(overlap.g, x2$gname), ]
  t1 = x1$time
  t2 = x2$time
  x1.rhythm = x1$rhythm[match(overlap.g, x1$rhythm$gname), ]
  x2.rhythm = x2$rhythm[match(overlap.g, x2$rhythm$gname), ]
  stopifnot("x$x1$P is not equal to x$x2$P. " = x$x1$P==x$x2$P)
  period = x$x1$P

  x.list = lapply(seq_along(overlap.g), function(a){
    list(x1.time = t1,
         x2.time = t2,
         y1 = as.numeric(x1.overlap[a, ]),
         y2 = as.numeric(x2.overlap[a, ]))
  })

  if(Par == "A"|Par == "phase"|Par == "M"){
    Par2 = switch(Par, "A" = "amplitude", "phase" = "phase", "M" = "basal")

    test_diffPar = parallel::mclapply(seq_along(x.list), function(a){
      out.diffPar = LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                                           x.list[[a]]$x2.time, x.list[[a]]$y2,
                                           period , FN = TRUE, type=Par2)
      one.row = data.frame(gname = overlap.g[a],
                           Par1 = ifelse(Par2=="phase", x1.rhythm[, "peak"][a],
                                         x1.rhythm[, Par][a]),
                           Par2 = ifelse(Par2=="phase", x2.rhythm[, "peak"][a],
                                         x2.rhythm[, Par][a]),
                           delta.Par = ifelse(Par2=="phase",
                                              x2.rhythm[, "peak"][a]-x1.rhythm[, "peak"][a],
                                              x2.rhythm[, Par][a]-x1.rhythm[, Par][a]),
                           pvalue = out.diffPar$pvalue
      )
      #notice that the peak estimate from diffCircadian and DCP_Rhythmicity is peak_dfCircadain + peak_CP = 30
      return(one.row)
    }, mc.cores = parallel.ncores)
    diffPar.tab = do.call(rbind.data.frame, test_diffPar)
    diffPar.tab$qvalue= stats::p.adjust(diffPar.tab$pvalue, p.adjust.method)
    colnames(diffPar.tab)[2:4] = gsub("Par", Par, colnames(diffPar.tab)[2:4])

  }else{
    data =  cbind.data.frame(x1.overlap, x2.overlap)
    time = c(t1, t2)
    group = c(rep(1, ncol(x1.overlap)), rep(2, ncol(x2.overlap)))
    group.1 = group == 1
    group.2 = group == 2

    design.full = data.frame(
      M1 = 1,
      g2 = group.2,
      inphase = cos(2*pi/period*time), outphase = sin(2*pi/period*time),
      g2.inphase = group.2*cos(2*pi/period*time),
      g2.outphase = group.2*sin(2*pi/period*time))
    fit.full = limma::lmFit(data,
                            stats::model.matrix(~g2+inphase+outphase+g2.inphase+g2.outphase,
                                                data = design.full))
    fit.full = limma::eBayes(fit.full)

    if(Par == "A&phase"){
      diff.top = limma::topTable(fit.full, coef = 5:6, n = nrow(data), sort.by = "none")
      test_diffPar = parallel::mclapply(seq_along(x.list), function(a){
        out.diffA = LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                                           x.list[[a]]$x2.time, x.list[[a]]$y2,
                                           period , FN = TRUE, type="amplitude")
        out.diffphase= LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                                              x.list[[a]]$x2.time, x.list[[a]]$y2,
                                              period , FN = TRUE, type="phase")
        one.row = data.frame(
          A1 = x1.rhythm$A[a],
          A2 = x2.rhythm$A[a],
          delta.A = x2.rhythm$A[a]-x1.rhythm$A[a],
          p.delta.A = out.diffA$pvalue,
          phase1 = x1.rhythm$peak[a],
          phase2 = x2.rhythm$peak[a],
          delta.phase = x2.rhythm$peak[a]-x1.rhythm$peak[a],
          p.delta.phase = out.diffphase$pvalue
        )
        return(one.row)
      }, mc.cores = parallel.ncores)
      diffPar.tab = do.call(rbind.data.frame, test_diffPar)
      all.tab = data.frame(gname = overlap.g,
                           diff.top[, c("P.Value", "adj.P.Val")], diffPar.tab)
      colnames(all.tab)[2:3] = c("p.overall", "q.overall")
      all.tab$q.delta.A = stats::p.adjust(all.tab$p.delta.A, p.adjust.method)
      all.tab$q.delta.phase = stats::p.adjust(all.tab$p.delta.phase, p.adjust.method)
      all.tab$post.hoc.A.By.q = all.tab$q.overall<alpha&all.tab$q.delta.A<alpha.adjust(alpha, 2, method = "sidak")
      all.tab$post.hoc.phase.By.q = all.tab$q.overall<alpha&all.tab$q.delta.phase<alpha.adjust(alpha, 2, method = "sidak")
      all.tab$post.hoc.A.By.p = all.tab$p.overall<alpha&all.tab$p.delta.A<alpha.adjust(alpha, 2, method = "sidak")
      all.tab$post.hoc.phase.By.p = all.tab$p.overall<alpha&all.tab$p.delta.phase<alpha.adjust(alpha, 2, method = "sidak")
      diffPar.tab = all.tab[, c("gname", "p.overall", "q.overall","post.hoc.A.By.p","post.hoc.phase.By.p", "post.hoc.A.By.q", "post.hoc.phase.By.q",
                                "A1", "A2", "delta.A", "p.delta.A", "q.delta.A", "phase1", "phase2", "delta.phase", "p.delta.phase", "q.delta.phase")]
    }else if(Par == "A&phase&M"){
      diff.top = limma::topTable(fit.full, coef = c(2, 5:6), n = nrow(data), sort.by = "none")
      test_diffPar = parallel::mclapply(seq_along(x.list), function(a){
        out.diffA = LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period , FN = TRUE, type="amplitude")
        out.diffphase= LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period , FN = TRUE, type="phase")
        out.diffM= LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period , FN = TRUE, type="basal")
        one.row = data.frame(
          A1 = x1.rhythm$A[a],
          A2 = x2.rhythm$A[a],
          delta.A = x2.rhythm$A[a]-x1.rhythm$A[a],
          p.delta.A = out.diffA$pvalue,
          phase1 = x1.rhythm$peak[a],
          phase2 = x2.rhythm$peak[a],
          delta.phase = x2.rhythm$peak[a]-x1.rhythm$peak[a],
          p.delta.phase = out.diffphase$pvalue,
          M1 = x1.rhythm$M[a],
          M2 = x2.rhythm$M[a],
          delta.M = x2.rhythm$M[a]-x1.rhythm$M[a],
          p.delta.M = out.diffM$pvalue
        )
        return(one.row)
      }, mc.cores = parallel.ncores)
      diffPar.tab = do.call(rbind.data.frame, test_diffPar)
      all.tab = data.frame(gname = overlap.g, diff.top[, c("P.Value", "adj.P.Val")], diffPar.tab)
      colnames(all.tab)[2:3] = c("p.overall", "q.overall")
      all.tab$q.delta.A = stats::p.adjust(all.tab$p.delta.A, p.adjust.method)
      all.tab$q.delta.phase = stats::p.adjust(all.tab$p.delta.phase, p.adjust.method)
      all.tab$q.delta.M = stats::p.adjust(all.tab$p.delta.M, p.adjust.method)
      all.tab$post.hoc.A.By.q = all.tab$q.overall<alpha&all.tab$q.delta.A<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.phase.By.q = all.tab$q.overall<alpha&all.tab$q.delta.phase<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.M.By.q = all.tab$q.overall<alpha&all.tab$q.delta.M<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.A.By.p = all.tab$p.overall<alpha&all.tab$p.delta.A<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.phase.By.p = all.tab$p.overall<alpha&all.tab$p.delta.phase<alpha.adjust(alpha, 3, method = "sidak")
      all.tab$post.hoc.M.By.p = all.tab$p.overall<alpha&all.tab$p.delta.M<alpha.adjust(alpha, 3, method = "sidak")
      diffPar.tab = all.tab[, c("gname", "p.overall", "q.overall","post.hoc.A.By.p","post.hoc.phase.By.p", "post.hoc.M.By.p",
                                "post.hoc.A.By.q", "post.hoc.phase.By.q", "post.hoc.M.By.q",
                                "A1", "A2", "delta.A", "p.delta.A", "q.delta.A", "phase1", "phase2", "delta.phase", "p.delta.phase", "q.delta.phase",
                                "M1", "M2", "delta.M", "p.delta.M", "q.delta.M")]
    }

  }

  colnames(diffPar.tab) = gsub("phase", "peak", colnames(diffPar.tab))
  diffPar.tab$P = period
  out = diffPar.tab
  return(out)
}
