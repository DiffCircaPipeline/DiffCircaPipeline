#' Rhythmicity Analysis with Cosinor Model
#'
#' This function either takes single-group data and performs only rhythmicity analysis, or takes two-group data and also categorize the genes into types of joint rhythmicity (TOJR).
#'
#' @param x1 group I data. A list with the following components: \itemize{
#' \item data: data.frame genes in rows and samples in columns.
#' \item time: time of expression. Should be in the same order as samples in the data.
#' \item gname: labels for gene. Should be in the same order as rows in the data.
#' }
#' @param x2 group II data. Components are same as x1.
#' @param method character string specifying the algorithm used for joint rhythmicity categorization. Should be one of "Sidak_FS", "Sidak_BS", "VDA", "AWFisher". Default "Sidak_FS" and is recommended.
#' @param period numeric. The length of the oscillation cycle. Default is 24 for circadian rhythm.
#' @param alpha numeric. Threshold for rhythmicity p-value in joint rhythmicity categorization. If CI = TRUE, (1-alpha) confidence interval for parameters will be returned.
#' @param alpha.FDR numeric. Threshold for rhythmicity p-value in joint rhythmicity categorization adjusted for global FDR control.
#' @param CI logical. Should confidence interval for A, phase and M be returned?
#' @param p.adjust.method input for p.adjust() in R package `stat`.
#' @param parallel.ncores integer. Number of cores used if using parallel computing with \code{mclapply()}. Not functional for windows system.
#'
#' @return A list of original x input with rhythmicity analysis estimates. If given two data sets, types of joint rhythmicity will also be available as the list component rhythm.joint.
#' @export
#' @details The methods "Sidak_FS" and "Sidak_BS" implement selective sequential model selection with Sidak adjusted p-value. "FS" represents forward stop, and "BS" basic stop, respectively (Fithian, W., et. al., 2015). The method "Sidak_FS" has better type I error control compared to venn diagram analysis (VDA) and adaptively weighted fisher's method (AWFisher).
#' @references Fithian, W., Taylor, J., Tibshirani, R., & Tibshirani, R. (2015). Selective sequential model selection. arXiv preprint arXiv:1512.02565.
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#'
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
DCP_Rhythmicity = function(x1, x2=NULL, method = "Sidak_FS", period = 24, alpha = 0.05, alpha.FDR = 0.05, CI = FALSE, p.adjust.method = "BH", parallel.ncores  = 1){

  x1 = CP_OneGroup(x1, period, alpha, CI, p.adjust.method)
  # x1 = CP_OneGroup(x1, period=24, alpha=0.05, CI=FALSE, p.adjust.method="BH")
  if(is.null(x2)){
    return(x1)
  }else{
    gname.overlap = intersect(x1$gname, x2$gname)
    stopifnot("There are no overlapping genes between x1$gname and x2$gname. " = length(gname.overlap)>0)

    x2 = CP_OneGroup(x2, period, alpha, CI, p.adjust.method)
    # x2 = CP_OneGroup(x2, period=24, alpha=0.05, CI=FALSE, p.adjust.method="BH")

    pM = data.frame(pG1 = x1$rhythm[match(gname.overlap, x1$rhythm$gname), "pvalue"],
                    pG2 = x2$rhythm[match(gname.overlap, x2$rhythm$gname), "pvalue"])#pG1 is p-value for group 1;
    action1 = apply(pM, 1, which.min)
    action2 = ifelse(action1==1, 2, 1)
    action = data.frame(action1, action2)
    pv = data.frame(pS1 = sapply(1:length(action1), function(i){pM[i, action1[i]]}),
                    pS2 = sapply(1:length(action2), function(i){pM[i, action2[i]]}))#pS1 is p-value for step 1;
    # The model selection procedure---
    TOJR = unlist(parallel::mclapply(1:nrow(action), function(i){
      SeqModelSel(action[i, ], pv[i, ], alpha, method)
    }, mc.cores = parallel.ncores))
    #NAME:TOJR:type of joint rhythmicity

    x = list(x1 = x1, x2 = x2,
             gname_overlap = gname.overlap,
             rhythm.joint = cbind.data.frame(gname.overlap, action, pM, TOJR))
    colnames(x$rhythm.joint) = c("gname", "action1", "action2", "pG1", "pG2", "TOJR")
    x$rhythm.joint$TOJR.FDR = toTOJR(x, method, alpha.FDR, adjustP = TRUE, p.adjust.method, parallel.ncores)
    return(x)
  }
}

CP_OneGroup = function(x1, period = 24, alpha = 0.05, CI = FALSE, p.adjust.method = "BH"){
  stopifnot("x1$data must be dataframe" = is.data.frame(x1$data))
  stopifnot("Number of samples in data does not match that in time. " = ncol(x1$data)==length(x1$time))
  stopifnot("Please input the gene labels x1$gname. " = !is.null(x1$gname))
  stopifnot("Number of gnames does not match number of genes in data. " = nrow(x1$data)==length(x1$gname))
  data = x1$data
  time = x1$time
  gname = x1$gname

  #design matrix
  design.vars = data.frame(inphase = cos(2 * pi * time / period),
                           outphase = sin(2 * pi * time / period))
  design = stats::model.matrix(~inphase+outphase, data = design.vars)
  fit = limma::lmFit(data, design)
  # fit = limma::eBayes(fit, trend = TRUE, robust = TRUE)
  fit = limma::eBayes(fit) #the default setting has better power.
  top = limma::topTable(fit, coef = 2:3, n = nrow(data), sort.by = "none")
  m.top = limma::topTable(fit, coef = 1, n = nrow(data), sort.by = "none")
  all.top = limma::topTable(fit, coef = 1:3, n = nrow(data), sort.by = "none")

  A.hat = apply(top[, 1:2], 1, function(i){sqrt(i[1]^2 + i[2]^2)})
  phase.hat = apply(top[, 1:2], 1, function(i){get_phase(i[1], i[2])$phase})
  # x.predict = t(design%*%t(all.top[, 1:3]))
  # RSS = apply((data-x.predict)^2, 1, sum)/(ncol(data)-3)
  RSS = fit$sigma^2
  TSS = apply(data, 1, function(i){sum((i-mean(i))^2)})

  R2 = 1-RSS*(length(time)-3)/TSS

  if(CI){
    # fit@.Data[[7]] is variance
    CI.m.hat.radius = sapply(1:length(fit$sigma), function(i){
      calculate_CI.M(fit@.Data[[7]], A.t = matrix(c(1, 0, 0), nrow = 1),
                     r.full = 3, ncol(data), alpha, fit$sigma[i])
    })
    se.hat.A.phase = t(sapply(1:length(fit$sigma), function(i){
      calculate_CI_A.phase.Taylor(fit@.Data[[7]],
                                  A.t = rbind(c(0, 1, 0),
                                              c(0, 0, 1)),
                                  phase.hat[i], A.hat[i], fit$sigma[i])
    }))
    se.A.hat = unlist(se.hat.A.phase[,"se.A.hat"]); se.phase.hat = unlist(se.hat.A.phase[, "se.phase.hat"])
    se.phase.hat = ifelse(se.phase.hat>(pi/2), pi/2, se.phase.hat)
    se.qt = stats::qt(1-alpha/2, ncol(data)-3)

    rhythm = data.frame(gname = gname,
                        M = m.top$logFC,
                        A = A.hat,
                        phase = phase.hat,
                        peak = (period-period*phase.hat/(2*pi))%%period,
                        M.ll = m.top$logFC-CI.m.hat.radius, M.ul = m.top$logFC+CI.m.hat.radius,
                        A.ll = A.hat-se.qt*se.A.hat, A.ul = A.hat+se.qt*se.A.hat,
                        phase.ll = phase.hat - se.qt*se.phase.hat, phase.ul = phase.hat + se.qt*se.phase.hat,
                        pvalue = top$P.Value,
                        qvalue = top$adj.P.Val,
                        sigma = fit$sigma, R2 = R2)

  }else{
    rhythm = data.frame(gname = gname,
                        M = m.top$logFC,
                        A = A.hat,
                        phase = phase.hat,
                        peak = (period-period*phase.hat/(2*pi))%%period,
                        pvalue = top$P.Value,
                        qvalue = top$adj.P.Val,
                        sigma = fit$sigma, R2 = R2)

  }
  x1$rhythm = rhythm
  x1$P = period
  return(x1)


}

#' Title Types of Joint Rhythmicity (TOJR)
#'
#' Categorize genes to four types of joint rhythmicity (TOJR).
#' @param x one of the following two: \itemize{
#' \item output of DCP_rhythmicity(x1, x2), both x1 and x2 are not NULL. (A list with the rhythm.joint component).
#' \item A list of with two outputs from DCP_rhythmicity(x1, x2 = NULL), each using data from a different group. }
#' @param method character string specifying the algorithm used for joint rhythmicity categorization. Should be one of "Sidak_FS", "Sidak_BS", "VDA", "AWFisher".
#' @param alpha integer. Threshold for rhythmicity p-value in joint rhythmicity categorization.
#' @param adjustP logic. Should joint rhythmicity categorization be based on adjusted p-value?
#' @param p.adjust.method input for p.adjust() in R package \code{stat}
#' @param parallel.ncores integer. Number of cores used if using parallel computing with \code{mclapply()}. Not functional for windows system.
#'
#' @return Joint rhythmicity categories of genes
#' @export
#'
#' @examples
#' #Re-calculate TOJR for DCP_rhythmicity(x1, x2) output with q-value cutoff 0.1
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' TOJR.new = toTOJR(rhythm.res, alpha = 0.1, adjustP = TRUE)
#'
#' #Calculate TOJR for two DCP_rhythmicity(x1, x2 = NULL) outputs with p-value cutoff 0.05
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res1 = DCP_Rhythmicity(x1 = x[[1]])
#' rhythm.res2 = DCP_Rhythmicity(x1 = x[[2]])
#' TOJR = toTOJR(x = list(x1 = rhythm.res1, x2 = rhythm.res2),
#' alpha = 0.05, adjustP = FALSE)
toTOJR = function(x, method = "Sidak_FS", alpha = 0.05, adjustP = TRUE, p.adjust.method = "BH", parallel.ncores = 1){

  if(is.null(x$rhythm.joint)){
    stopifnot("Please see examples for correct x input" = (length(x)==2)&(!is.null(x[[1]]$rhythm))&(!is.null(x[[2]]$rhythm)))
    x1 = x[[1]]
    x2 = x[[2]]
    gname.overlap = intersect(x1$gname, x2$gname)

    pM = data.frame(pG1 = x1$rhythm[match(gname.overlap, x1$rhythm$gname), "pvalue"],
                    pG2 = x2$rhythm[match(gname.overlap, x2$rhythm$gname), "pvalue"])#pG1 is p-value for group 1;
    action1 = apply(pM, 1, which.min)
    action2 = ifelse(action1==1, 2, 1)
    action = data.frame(action1, action2)
    pv = data.frame(pS1 = sapply(1:length(action1), function(i){pM[i, action1[i]]}),
                    pS2 = sapply(1:length(action2), function(i){pM[i, action2[i]]}))#pS1 is p-value for step 1;

  }else{
    pM = as.data.frame(x$rhythm.joint[, c("pG1", "pG2")])
    action1 = x$rhythm.joint$action1
    action2 = x$rhythm.joint$action2
    action = data.frame(action1, action2)
    pv = data.frame(pS1 = sapply(1:length(action1), function(i){pM[i, action1[i]]}),
                    pS2 = sapply(1:length(action2), function(i){pM[i, action2[i]]}))#pS1 is p-value for step 1;
  }

  if(adjustP){
    #method1.2: stratified BH with group
    qM = data.frame(p1 = stats::p.adjust(pM$pG1, p.adjust.method), p2 = stats::p.adjust(pM$pG2, p.adjust.method)) #q.value from each group
    q.action1 = apply(qM, 1, which.min)
    q.action2 = ifelse(q.action1==1, 2, 1)
    qv = data.frame(qS1 = sapply(1:length(q.action1), function(i){qM[i, q.action1[i]]}),
                    qS2 = sapply(1:length(q.action2), function(i){qM[i, q.action2[i]]}))
    q.action = data.frame(q.action1, q.action2)
    TOJR_adj = unlist(parallel::mclapply(1:nrow(q.action), function(i){
      SeqModelSel(q.action[i, ], qv[i, ], alpha, method)
    }, mc.cores = parallel.ncores))
  }else{
    TOJR_adj = unlist(parallel::mclapply(1:nrow(action), function(i){
      # SeqModelSel(action[i, ], pv[i, ], alpha, method)
      SeqModelSel(action[i, ], pv[i, ], alpha, "Sidak_FS")
    }, mc.cores = parallel.ncores))
  }

  return(TOJR_adj)
}

# Functions for CI --------------------------------------------------------
get_phase = function(b1.x = beta1.hat, b2.x = beta2.hat){
  ph.x = atan(-b2.x/b1.x)
  #adjust ph.x
  if(b2.x>0){
    if(ph.x<0){
      ph.x = ph.x+2*pi
    }else if (ph.x>0){
      ph.x = ph.x+pi
    }
  }else if(b2.x<0){
    if(ph.x<0){
      ph.x = ph.x+pi
    }
  }else{
    ph.x = 88#88 means one the the beta estimate is 0. I did not account for such senerio because it is rare and complicated
  }
  return(list(phase = ph.x,
              tan = -b2.x/b1.x))
}

calculate_CI.M = function(XX.inv = mat.S.inv, A.t = matrix(c(1, 0, 0, 1, 0, 0), nrow = 1),
                          r.full = 6, n = length(tod), alpha = 0.05, sigma.hat = sigma.hat){
  CI.m.hat.radius = stats::qt(1-alpha/2, n-r.full)*sigma.hat*sqrt(A.t%*%XX.inv%*%t(A.t))
  return(CI.m.hat.radius)
}

calculate_CI_A.phase.Taylor = function(XX.inv = mat.S.inv,
                                       A.t = rbind(c(0, 0, 1, 0, 0, 0),
                                                   c(0, 0, 0, 1, 0, 0)),
                                       phase.hat = phase1.hat, A.hat = A1.hat,
                                       sigma.hat = sigma.hat){
  #notice that this has been changed from the one_consinor_OLS
  var.new = A.t%*%XX.inv%*%t(A.t)
  var.beta1 = var.new[1, 1]; var.beta2 = var.new[2, 2]; var.beta1.beta2 = var.new[1, 2]
  se.A.hat = sigma.hat*sqrt(var.beta1*cos(phase.hat)^2
                            -2*var.beta1.beta2*sin(phase.hat)*cos(phase.hat)
                            +var.beta2*sin(phase.hat)^2)
  se.phase.hat = sigma.hat*sqrt(var.beta1*sin(phase.hat)^2
                                +2*var.beta1.beta2*sin(phase.hat)*cos(phase.hat)
                                +var.beta2*cos(phase.hat)^2)/A.hat
  return(list(se.A.hat = se.A.hat,
              se.phase.hat = se.phase.hat))

}

calculate_CI_A.phase.Scheffe = function(XX.inv = mat.S.inv, A.t = rbind(c(0, 1, 0, 0, 1, 0),
                                                                        c(0, 0, 1, 0, 0, 1)),
                                        sigma2.hat = sigma2.hat,
                                        est = est,r.full = 6, n = length(tod), alpha, CItype = "conservative"){
  #CItype = "conservative" or "PlugIn"
  # XX.inv = mat.S.inv;
  # A.t = rbind(c(0, 1, 0),
  #             c(0, 0, 1))
  # r.full=3

  q = Matrix::rankMatrix(A.t)[[1]]
  est2 = A.t%*%est
  beta1.hat = est2[1]
  beta2.hat = est2[2]
  B.inv = solve(A.t%*%XX.inv%*%t(A.t))
  B11 = B.inv[1, 1]; B12 = B.inv[1, 2]; B22 = B.inv[2, 2]
  R = q*sigma2.hat*stats::qf(1-alpha, q, n-r.full)

  C1 = -(B11*beta1.hat+B12*beta2.hat)/(B22*beta2.hat+B12*beta1.hat)
  C2 = -(R-2*B12*beta1.hat*beta2.hat-B11*beta1.hat^2-B22*beta2.hat^2)/(B22*beta2.hat+B12*beta1.hat)
  D1 = B22*C1^2+B11+2*B12*C1
  D2 = 2*B22*C1*C2+2*B12*C2-(B12*beta1.hat+B22*beta2.hat)*C1-B11*beta1.hat-B12*beta2.hat
  D3 = B22*C2^2-(B12*beta1.hat+B22*beta2.hat)*C2

  #calculate CI of phi
  #check if 0 is in ellipse
  zero.in.ellipse = B11*beta1.hat^2+2*B12*beta1.hat*beta2.hat+B22*beta2.hat^2 < R
  if(CItype == "conservative"){
    if(!zero.in.ellipse){
      delta.poly = D2^2-4*D1*D3
      if(delta.poly<0){
        phi.lower.limit = list(tan = -99, phase = -99)
        phi.upper.limit = list(tan = -99, phase = -99)
      }else{
        phi.beta1.roots = c((-D2-sqrt(delta.poly))/(2*D1), (-D2+sqrt(delta.poly))/(2*D1))
        phi.beta2.roots = C1*phi.beta1.roots+C2

        # phi.limit1 = get_phase(phi.beta1.roots[1], phi.beta2.roots[1])
        # phi.limit2 = get_phase(phi.beta1.roots[2], phi.beta2.roots[2])

        #change to adjusted phi limits
        phi.limits = get_phaseForCI(b1.x = beta1.hat, b2.x = beta2.hat,
                                    b1.r1 = phi.beta1.roots[1], b2.r1 = phi.beta2.roots[1],
                                    b1.r2 = phi.beta1.roots[2], b2.r2 = phi.beta2.roots[2])
        phi.lower.limit = phi.limits$phi.lower.limit
        phi.upper.limit = phi.limits$phi.upper.limit
      }
    }else{
      phi.lower.limit = list(tan = 99, phase = 99)
      phi.upper.limit = list(tan = 99, phase = 99)
    }

    #calculate CI of A
    #calculate normal ellipse parameters
    ellipse.parameters = solve.ellipse.parameters(a = B11, b = B22, c = 2*B12,
                                                  d = -2*(B11*beta1.hat+B12*beta2.hat),
                                                  e = -2*(B12*beta1.hat+B22*beta2.hat),
                                                  f = B11*beta1.hat^2+B22*beta2.hat^2+2*B12*beta1.hat*beta2.hat-q*sigma2.hat*stats::qf(1-alpha, q, n-r.full))
    angle.point.newOrigin.major =
      get.angle_point.newOrigin.major(x0 = ellipse.parameters$x0,
                                      y0 = ellipse.parameters$y0,
                                      theta.rotate = ellipse.parameters$theta.rotate)
    r.point.to.center = sqrt(ellipse.parameters$x0^2+ellipse.parameters$y0^2)
    x.new = r.point.to.center*cos(angle.point.newOrigin.major)
    y.new = r.point.to.center*sin(angle.point.newOrigin.major)

    #check the position of the new origin to the ellipse
    if(x.new==0&y.new==0){
      A.limit1 = ellipse.parameters$minor
      A.limit2 = ellipse.parameters$major
    }else if(x.new==0){
      A.limit1 = abs(ellipse.parameters$minor-y.new)
      A.limit2 = ellipse.parameters$minor+y.new
    }else if(y.new==0){
      A.limit1 = abs(ellipse.parameters$major-x.new)
      A.limit2 = ellipse.parameters$major+x.new
    }else if(x.new!=0&y.new!=0){
      x.new = abs(x.new)
      y.new = abs(y.new)
      fun.t = function(t){
        (ellipse.parameters$major*x.new/(t+ellipse.parameters$major^2))^2+
          (ellipse.parameters$minor*y.new/(t+ellipse.parameters$minor^2))^2-1
      }

      # root.upper = 50
      # while(fun.t(-ellipse.parameters$minor^2+0.0001)*fun.t(root.upper)>0){
      #   root.upper = root.upper+50
      # }
      #    root1 = uniroot(fun.t, c(-ellipse.parameters$minor^2, root.upper),extendInt="downX")$root
      root1 = stats::uniroot(fun.t, c(-ellipse.parameters$minor^2, -ellipse.parameters$minor^2+50),extendInt="downX")$root
      x.root1 = ellipse.parameters$major^2*x.new/(root1+ellipse.parameters$major^2)
      y.root1 = ellipse.parameters$minor^2*y.new/(root1+ellipse.parameters$minor^2)
      dmin = sqrt((x.new-x.root1)^2+(y.new-y.root1)^2)

      # root.lower = -50
      # while(fun.t(-ellipse.parameters$major^2-0.0001)*fun.t(root.lower)>0){
      #   root.lower = root.lower-50
      # }
      #root2 = uniroot(fun.t, c(root.lower, -ellipse.parameters$major^2-0.0001), extendInt="yes")$root
      root2 = stats::uniroot(fun.t, c(-ellipse.parameters$major^2-50, -ellipse.parameters$major^2), extendInt="upX")$root
      # root2 = uniroot(fun.t, c(-ellipse.parameters$major^2-50, -ellipse.parameters$major^2), extendInt="yes")$root
      # root2 = uniroot(fun.t, c(-ellipse.parameters$major^2-50, -ellipse.parameters$major^2), extendInt="no")$root
      x.root2 = ellipse.parameters$major^2*x.new/(root2+ellipse.parameters$major^2)
      y.root2 = ellipse.parameters$minor^2*y.new/(root2+ellipse.parameters$minor^2)
      dmax = sqrt((x.new-x.root2)^2+(y.new-y.root2)^2)

      A.limit1 = min(dmin, dmax)
      A.limit2 = max(dmin, dmax)
    }
  }else if(CItype=="PlugIn"){
    if(!zero.in.ellipse){
      delta.poly = D2^2-4*D1*D3
      if(delta.poly<0){
        phi.lower.limit = list(tan = -99, phase = -99)
        phi.upper.limit = list(tan = -99, phase = -99)
      }else{
        #get the roots through the rootSolve package
        #install.packages("rootSolve")
        phi.beta1.roots = c((-D2-sqrt(delta.poly))/(2*D1), (-D2+sqrt(delta.poly))/(2*D1))
        phi.beta2.roots = C1*phi.beta1.roots+C2

        # model = function(x){
        #   beta1 = x[1]; beta2 = x[2]
        #   F1 = B11*(beta1-beta1.hat)^2+2*B12*(beta1-beta1.hat)*(beta2-beta2.hat)+B22*(beta2-beta2.hat)^2-R
        #   F2 = beta1^2+beta2^2-beta1.hat^2-beta2.hat^2
        #   return(c(F1 = F1, F2 = F2))
        # }
        #
        # start.points = list(c(phi.beta1.roots[1], phi.beta2.roots[1]),
        #                     c(phi.beta1.roots[2], phi.beta2.roots[2]),
        #                     c(beta1.hat/2, beta2.hat/2),
        #                     c(beta1.hat/2, beta2.hat*2),
        #                     c(beta1.hat*2, beta2.hat/2),
        #                     c(beta1.hat*2, beta2.hat*2))
        #
        # solutions.A.PlugIn = lapply(start.points, function(a){
        #   res = tryCatch(rootSolve::multiroot(f=model, start = a , maxiter = 100)$root, warning=function(cnd){NA})
        #   return(res)
        # })
        #
        # dist.solutions = as.matrix(dist(do.call(rbind, solutions.A.PlugIn)))
        # solution1.A.PlugIn = solutions.A.PlugIn[[1]]
        # solution2.A.PlugIn = solutions.A.PlugIn[[which(round(dist.solutions[, 1], 1)!=0)[1]]]
        A.hat2 = sqrt(beta1.hat^2+beta2.hat^2)

        fun = function(phi){
          B11*A.hat2^2*cos(phi)^2+B22*A.hat2^2*sin(phi)^2+2*B12*A.hat2^2*sin(phi)*cos(phi)-
            (2*B11*beta1.hat*A.hat2+2*B12*beta2.hat*A.hat2)*cos(phi)-
            (2*B12*beta1.hat*A.hat2+2*B22*beta2.hat*A.hat2)*sin(phi)+
            B11*beta1.hat^2+2*B12*beta1.hat*beta2.hat+B22*beta2.hat^2-R
        }
        all.solutions = rootSolve::uniroot.all(fun, c(0, 2*pi))

        solution1.A.PlugIn = c(A.hat2*cos(all.solutions[1]), A.hat2*sin(all.solutions[1]))
        solution2.A.PlugIn = c(A.hat2*cos(all.solutions[2]), A.hat2*sin(all.solutions[2]))


        #change to adjusted phi limits
        phi.limits = get_phaseForCI(b1.x = beta1.hat, b2.x = beta2.hat,
                                    b1.r1 = solution1.A.PlugIn[1], b2.r1 = solution1.A.PlugIn[2],
                                    b1.r2 = solution2.A.PlugIn[1], b2.r2 = solution2.A.PlugIn[2])
        phi.lower.limit = phi.limits$phi.lower.limit
        phi.upper.limit = phi.limits$phi.upper.limit
      }
    }else{
      phi.lower.limit = list(tan = 99, phase = 99)
      phi.upper.limit = list(tan = 99, phase = 99)
    }

    #calculate CI of A
    #calculate normal ellipse parameters
    ellipse.parameters = solve.ellipse.parameters(a = B11, b = B22, c = 2*B12,
                                                  d = -2*(B11*beta1.hat+B12*beta2.hat),
                                                  e = -2*(B12*beta1.hat+B22*beta2.hat),
                                                  f = B11*beta1.hat^2+B22*beta2.hat^2+2*B12*beta1.hat*beta2.hat-q*sigma2.hat*stats::qf(1-alpha, q, n-r.full))
    angle.point.newOrigin.major =
      get.angle_point.newOrigin.major(x0 = ellipse.parameters$x0,
                                      y0 = ellipse.parameters$y0,
                                      theta.rotate = ellipse.parameters$theta.rotate)
    r.point.to.center = sqrt(ellipse.parameters$x0^2+ellipse.parameters$y0^2)
    x.new = r.point.to.center*cos(angle.point.newOrigin.major)
    y.new = r.point.to.center*sin(angle.point.newOrigin.major)

    #check the position of the new origin to the ellipse
    if(x.new==0&y.new==0){
      A.limit1 = ellipse.parameters$minor
      A.limit2 = ellipse.parameters$major
    }else if(x.new==0){
      A.limit1 = abs(ellipse.parameters$minor-y.new)
      A.limit2 = ellipse.parameters$minor+y.new
    }else if(y.new==0){
      A.limit1 = abs(ellipse.parameters$major-x.new)
      A.limit2 = ellipse.parameters$major+x.new
    }else if(x.new!=0&y.new!=0){
      x.new = abs(x.new)
      y.new = abs(y.new)
      r1.xx2 = ellipse.parameters$major^2*ellipse.parameters$minor^2/(ellipse.parameters$major^2+ellipse.parameters$minor^2*ellipse.parameters$tan.rotate^2)
      r1 = sqrt(r1.xx2+r1.xx2*ellipse.parameters$tan.rotate^2)
      dmin = sqrt(x.new^2+y.new^2)-r1
      dmax = sqrt(x.new^2+y.new^2)+r1

      A.limit1 = min(dmin, dmax)
      A.limit2 = max(dmin, dmax)
    }
  }

  if(zero.in.ellipse){
    A.limit1 = 0
  }

  return(list(CI_phase = c(phi.lower.limit$phase, phi.upper.limit$phase),
              CI_A = c(A.limit1, A.limit2)))

}
get_phaseForCI = function(b1.x = beta1.hat, b2.x = beta2.hat,
                          b1.r1 = phi.beta1.roots[1], b2.r1 = phi.beta2.roots[1],
                          b1.r2 = phi.beta1.roots[2], b2.r2 = phi.beta2.roots[2]){
  # b1.x = beta1.hat; b2.x = beta2.hat;
  # b1.r1 = phi.beta1.roots[1]; b2.r1 = phi.beta2.roots[1];
  # b1.r2 = phi.beta1.roots[2]; b2.r2 = phi.beta2.roots[2]

  #b1.r1 is the beta1 estimate of root 1.
  #step 1: find which root is in the same quadrant as OLS estimate b1.x, b2.x

  x.quad = get_quad(b1.q = b1.x, b2.q = b2.x)
  phase.res = get_phase(b1.x, b2.x)
  r1.quad = get_quad(b1.q = b1.r1, b2.q = b2.r1)
  r2.quad = get_quad(b1.q = b1.r2, b2.q = b2.r2)

  ##the easy scenario: three are in the same quadrant
  if(r1.quad==x.quad&r2.quad==x.quad){
    phi.limit1 = get_phase(b1.r1, b2.r1)
    phi.limit2 = get_phase(b1.r2, b2.r2)
    if(phi.limit1$phase<phi.limit2$phase){
      phi.lower.limit = phi.limit1
      phi.upper.limit = phi.limit2
    }else{
      phi.lower.limit = phi.limit2
      phi.upper.limit = phi.limit1
    }
  }else if(r1.quad == x.quad|r2.quad==x.quad){
    #more complicated scenario: one of the root is in the other quadrant
    if(r1.quad == x.quad){
      b1.s1 = b1.r1; b2.s1 = b2.r1 #s1 stands for the root that is in the same quadrant
      b1.s2 = b1.r2; b2.s2 = b2.r2
      quad.s1 = r1.quad; quad.s2 = r2.quad
    }else{
      b1.s1 = b1.r2; b2.s1 = b2.r2 #s1 stands for the root that is in the same quadrant
      b1.s2 = b1.r1; b2.s2 = b2.r1
      quad.s1 = r2.quad; quad.s2 = r1.quad
    }
    phi.s1 = get_phase(b1.s1, b2.s1)
    if(phi.s1$phase<phase.res$phase){
      phi.lower.limit = phi.s1
      s2.type = "s1 lower, s2 upper"
      phi.upper.limit = get_phaseForCI_QuadTab(s1.quad = quad.s1, s2.quad = quad.s2, s1s2 = s2.type,
                                               b1.s2.q = b1.s2, b2.s2.q = b2.s2)
    }else{
      phi.upper.limit = phi.s1
      s2.type = "s1 upper, s1 lower"
      phi.lower.limit = get_phaseForCI_QuadTab(s1.quad = quad.s1, s2.quad = quad.s2, s1s2 = s2.type,
                                               b1.s2.q = b1.s2, b2.s2.q = b2.s2)
    }
  }else{
    #this is the scenario that the ellipse covers three quadrants
    phi.limits = get_phaseForCI_QuadTab2(x.quad, r1.quad, r2.quad,
                                         b1.t1 = b1.r1, b2.t1 = b2.r1,
                                         b1.t2 = b1.r2, b2.t2 = b2.r2)
    phi.lower.limit = phi.limits$phi.lower.limit
    phi.upper.limit = phi.limits$phi.upper.limit
  }

  return(list(phi.lower.limit = phi.lower.limit,
              phi.upper.limit = phi.upper.limit))

}

#make a table of the quadrants
get_quad = function(b1.q = b1.x, b2.q = b2.x){
  quad.table = as.data.frame(
    list(beta1.gt.0 = c(TRUE, TRUE, FALSE, FALSE),
         beta2.gt.0 = c(TRUE, FALSE, TRUE, FALSE),
         quad = c(1, 4, 2, 3))
  )
  quad = quad.table[(quad.table$beta1.gt.0==(b1.q>0))&(quad.table$beta2.gt.0==(b2.q>0)),
                    "quad"]
  return(quad)
}

# s1.quad = quad.s1; s2.quad = quad.s2; s1s2 = s2.type;
# b1.s2.q = b1.s2; b2.s2.q = b2.s2
get_phaseForCI_QuadTab = function(s1.quad = quad.s1, s2.quad = quad.s2, s1s2 = s2.type,
                                  b1.s2.q = b1.s2, b2.s2.q = b2.s2){
  # s1.quad = quad.s1; s2.quad = quad.s2; s1s2 = s2.type;
  # b1.s2.q = b1.s2; b2.s2.q = b2.s2
  quad.table2 = as.data.frame(
    list(s1.quad = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4),
         s2.quad = c(2, 2, 3, 3, 4, 4, 1, 1, 3, 3, 4, 4, 1, 1, 2, 2, 4, 4, 1, 1, 2, 2, 3, 3),
         s1s2 = rep(c("s1 lower, s2 upper", "s1 upper, s1 lower"), 12),
         s2.low = c(3, 1, 5/2, 1/2, 2, 0, 3/2, -1/2, 5/2, 1/2, 2, 0, 3/2, -1/2, 1, -1, 2, 0, 3/2, -1/2, 1, -1, 1/2, -3/2),
         s2.up  = c(7/2, 3/2, 3, 1, 5/2, 1/2, 2, 0, 3, 1, 5/2, 1/2, 2, 0, 3/2, -1/2, 5/2, 1/2, 2, 0, 3/2, -1/2, 1, -1))
  )

  s2.low = quad.table2[quad.table2$s1.quad==s1.quad&quad.table2$s2.quad==s2.quad&quad.table2$s1s2==s1s2,
                       "s2.low"]
  s2.up  = quad.table2[quad.table2$s1.quad==s1.quad&quad.table2$s2.quad==s2.quad&quad.table2$s1s2==s1s2,
                       "s2.up"]

  phi.s2 = atan(-b2.s2.q/b1.s2.q)
  if(phi.s2<s2.low*pi){
    while(phi.s2<s2.low*pi){
      phi.s2 = phi.s2+pi
    }
  }else if(phi.s2>s2.up*pi){
    while(phi.s2>s2.up*pi){
      phi.s2 = phi.s2-pi
    }
  }

  # if((phi.s2<(s2.up*pi))&(phi.s2>(s2.low*pi))){
  #   ##print("#S2 is in interval")
  # }else{
  #   ##print("S2 is not in the int!!!!!!!!!!!!!!!!")
  # }
  return(list(phase = phi.s2,
              tan = -b2.s2.q/b1.s2.q))
}

get_phaseForCI_QuadTab2 = function(x.quad, r1.quad, r2.quad,
                                   b1.t1 = b1.r1, b2.t1 = b2.r1,
                                   b1.t2 = b1.r2, b2.t2 = b2.r2){
  #t1 stands for three-quadrants, solution 1
  #This function is for when the ellipse covers three quadrants
  quad.table3 = as.data.frame(
    list(x.quad = c(1, 1, 2, 2, 3, 3, 4, 4),
         t.quad = c(2, 4, 3, 1, 4, 2, 1, 3),
         t.low = c(1, 2, 1/2, 3/2, 0, 1, -1/2, 1/2),
         t.up = c(3/2, 5/2, 1, 2, 1/2, 3/2, 0, 1))
  )
  t1.low = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r1.quad, "t.low"]
  t1.up = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r1.quad, "t.up"]
  t2.low = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r2.quad, "t.low"]
  t2.up = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r2.quad, "t.up"]

  phi.t1 = atan(-b2.t1/b1.t1)
  if(phi.t1<t1.low*pi){
    while(phi.t1<t1.low*pi){
      phi.t1 = phi.t1+pi
    }
  }else if(phi.t1>t1.up*pi){
    while(phi.t1>t1.up*pi){
      phi.t1 = phi.t1-pi
    }
  }

  # if((phi.t1<(t1.up*pi))&(phi.t1>(t1.low*pi))){
  #  # #print("#t1 is in interval")
  # }else{
  #  # #print("t1 is not in the int!!!!!!!!!!!!!!!!")
  # }

  phi.t2 = atan(-b2.t2/b1.t2)
  if(phi.t2<t2.low*pi){
    while(phi.t2<t2.low*pi){
      phi.t2 = phi.t2+pi
    }
  }else if(phi.t2>t2.up*pi){
    while(phi.t2>t2.up*pi){
      phi.t2 = phi.t2-pi
    }
  }

  # if((phi.t2<(t2.up*pi))&(phi.t2>(t2.low*pi))){
  #   # #print("#t2 is in interval")
  # }else{
  #   # #print("t2 is not in the int!!!!!!!!!!!!!!!!")
  # }

  if(phi.t1<phi.t2){
    phi.lower.limit = list(phase = phi.t1,
                           tan = -b2.t1/b1.t1)
    phi.upper.limit = list(phase = phi.t2,
                           tan = -b2.t2/b1.t2)
  }else{
    phi.lower.limit = list(phase = phi.t2,
                           tan = -b2.t2/b1.t2)
    phi.upper.limit = list(phase = phi.t1,
                           tan = -b2.t1/b1.t1)
  }
  return(list(phi.lower.limit = phi.lower.limit,
              phi.upper.limit = phi.upper.limit))
}

solve.ellipse.parameters = function(a = a, b = b, c = c, d = d, e = e, f = f){
  #a * x ^ 2 + b * y ^ 2 + c * x * y + d * x + e * y + f = 0
  delta1 = c^2 -4*a*b
  if(delta1>=0){
    return("Not a proper epplise")
  }else{
    #the transformation function is from wikepedia
    denominator.factor1 = 2*(a*e^2+b*d^2-c*d*e+delta1*f)
    major.minor.square.diff = sqrt(((a-b)^2+c^2))
    denominator.factor2 = c(a+b+major.minor.square.diff, a+b-major.minor.square.diff)
    major.minor = -sqrt(denominator.factor1*denominator.factor2)/delta1
    x0 = (2*b*d-c*e)/delta1
    y0 = (2*a*e-c*d)/delta1
    if(c==0){
      if(a<b){
        theta = 0
        tan.rotate = 0
      }else{
        theta = pi/2
        tan.rotate = Inf
      }
    }else{
      tan.rotate = (b-a-major.minor.square.diff)/c
      theta = atan(tan.rotate)
    }
    return(list(major = major.minor[1],
                minor = major.minor[2],
                x0 = x0,
                y0 = y0,
                theta.rotate = theta,
                tan.rotate = tan.rotate))
  }
}

get.angle_point.newOrigin.major = function(x0 = ellipse.parameters$x0,
                                           y0 = ellipse.parameters$y0,
                                           theta.rotate = ellipse.parameters$theta.rotate){
  #if both x0=0 and y0=0, then the distance is already there:
  #min distance = b; max distance = a
  if(x0==0&y0==0){
    return("ellipse is centered at (0, 0)")
  }else if(y0==0){
    angle = theta.rotate
    return(angle)
  }else if(x0==0){
    angle = pi/2-theta.rotate
    return(angle)
  }else{
    tan.center.0.x_axis = y0/x0
    theta.center.0.x_axis = atan(tan.center.0.x_axis)
    if(theta.center.0.x_axis*theta.rotate<0){
      angle = abs(theta.center.0.x_axis)+abs(theta.rotate)
      return(angle)
    }else if(theta.center.0.x_axis*theta.rotate>0){
      angle = abs(theta.center.0.x_axis-theta.rotate)
      return(angle)
    }else if(theta.center.0.x_axis*theta.rotate==0){
      angle = abs(theta.center.0.x_axis)
      return(angle)
    }
  }
}

# Functions for model selection -------------------------------------------
#Fisher_BS, Fisher_FS, AWFisher_BS, AWFisher_FS, Sidak_BS, Sidak_FS, Nominal_BS, VDA, AWFisher, VDA
SeqModelSel = function(action = c(1, 2), pv = c(0.01, 0.02), alpha = 0.05, method = "Fisher_FS"){
  action = as.numeric(action)
  pv = as.numeric(pv)
  if(method == "Fisher_BS"){
    p_combined = fishers.p(pv)
    stop = StopRule(action, pv, alpha, "BS")
    if(p_combined<alpha){
      stop= ifelse(stop$type=="arrhy", stop$one, stop$type)
    }else{
      stop = "arrhy"
    }
  }else if(method == "Fisher_FS"){
    p_combined = fishers.p(pv)
    stop = StopRule(action, pv, alpha, "FS")
    if(p_combined<alpha){
      stop= ifelse(stop$type=="arrhy", stop$one, stop$type)
    }else{
      stop = "arrhy"
    }
  }else if(method == "AWFisher_BS"){
    p_combined = AWFisher::AWFisher_pvalue(pv)
    p_combined = p_combined$pvalues
    stop = StopRule(action, pv, alpha, "BS")
    if(p_combined<alpha){
      stop= ifelse(stop$type=="arrhy", stop$one, stop$type)
    }else{
      stop = "arrhy"
    }
  }else if(method == "AWFisher_FS"){
    p_combined = AWFisher::AWFisher_pvalue(pv)
    p_combined = p_combined$pvalues
    stop = StopRule(action, pv, alpha, "FS")
    if(p_combined<alpha){
      stop= ifelse(stop$type=="arrhy", stop$one, stop$type)
    }else{
      stop = "arrhy"
    }
  }else if(method == "Sidak_BS"){
    stop = StopRule(action, pv, alpha.adjust(alpha, length(pv), "sidak"), "BS")$type
  }else if(method == "Sidak_FS"){
    stop = StopRule(action, pv, alpha.adjust(alpha, length(pv), "sidak"), "FS")$type
  }else if(method == "Nominal_BS"){
    stop = StopRule(action, pv, alpha, "BS")$type
  }else if(method == "Nominal_FS"){
    stop = StopRule(action, pv, alpha, "FS")$type
  }else if(method == "VDA"){
    is.rhy = pv<alpha
    if(sum(is.rhy)==2){
      stop = "both"
    }else if(sum(is.rhy)==0){
      stop = "arrhy"
    }else if(is.rhy[1]){
      stop = "rhyI"
    }else{
      stop = "rhyII"
    }
  }else if(method == "AWFisher"){
    p_combined = AWFisher::AWFisher_pvalue(as.numeric(pv))
    aw.wight = p_combined$weights
    p_combined = p_combined$pvalues
    if(p_combined>alpha){
      stop = "arrhy"
    }else if(sum(aw.wight)==2){
      stop = "both"
    }else if(aw.wight[1, 1]==1){
      stop = "rhyI"
    }else{
      stop = "rhyII"
    }
  }
  return(stop)
}

StopRule = function(action, pv, alpha, method){
  if(method == "FS"){
    a.stop = forwardStop(pv, alpha)$stop
  }else if(method == "BS"){
    if(all(pv<alpha)){
      a.stop=length(action)
    }else{
      a.stop = max(0, min(which(pv>alpha))-1)
    }
  }
  sel = action[1:a.stop]
  if(a.stop==0){
    a.type = "arrhy"
  }else if(length(sel)==1){
    if(sel==1){
      a.type = "rhyI"
    }else{
      a.type = "rhyII"
    }
  }else{
    a.type = "both"
  }
  return(list(type = a.type,
              one = ifelse(action[1]==1, "rhyI", "rhyII")))
}

forwardStop = function (pv, alpha = 0.1)
{
  if (alpha < 0 || alpha > 1)
    stop("alpha must be in [0,1]")
  if (min(pv, na.rm = T) < 0 || max(pv, na.rm = T) > 1)
    stop("pvalues must be in [0,1]")
  val = -(1/(1:length(pv))) * cumsum(log(1 - pv))
  oo = which(val <= alpha)
  if (length(oo) == 0)
    out = 0
  else out = oo[length(oo)]
  return(list(stop = out,
              p = -(1/out) * sum(log(1 - pv))))
}

fishers.p = function(ps){
  fisher.chi = -2*sum(log(ps))
  p.fisher.chi = stats::pchisq(fisher.chi, 2*length(ps), lower.tail = FALSE)
  return(p.fisher.chi)
}

alpha.adjust = function(alpha, k, method = "bonferroni"){
  if(method == "bonferroni"){
    alpha.new = alpha/k
  }else if(method == "sidak"){
    alpha.new = 1-(1-alpha)^(1/k)
  }
  return(alpha.new)
}

fishers.p = function(ps){
  fisher.chi = -2*sum(log(ps))
  p.fisher.chi = stats::pchisq(fisher.chi, 2*length(ps), lower.tail = FALSE)
  return(p.fisher.chi)
}
