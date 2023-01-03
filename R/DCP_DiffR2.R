#' Differential Rhythmicity Fitness Test
#'
#' @param x outpuf of DCP_Rhythmicity(x1, x2)
#' @param method one of "LR", "permutation", "bootstrap". "LR" is recommended.
#' @param TOJR toTOJR output. If NULL, rhythm.joint object in x will be used.
#' @param alpha cutoff of p-values for significant differential rhythm fitness change.
#' @param nSampling number of samplings if "permutation" or "bootstrap" is chosen for method
#' @param Sampling.save directory to save sampling results. If NULL, the result will not be saved.
#' @param Sampling.file.label character string. Used as file lable if sampling results are saved.
#' @param p.adjust.method input for p.adjust() in R package \code{stat}
#' @param parallel.ncores integer. Number of cores used if using parallel computing with \code{mclapply()}. Not functional for windows system.
#'
#' @return A dataframe of differential rhythm finess test results.
#' @export
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' rhythm.diffR2 = DCP_DiffPar(rhythm.res)
DCP_DiffR2 = function(x, method = "LR", TOJR = NULL,  alpha = 0.05,  nSampling=1000, Sampling.save = NULL,
                     Sampling.file.label = "gII_vs_gI",p.adjust.method = "BH", parallel.ncores = 1){
  stopifnot('method should be one of "LR", "permutation", "bootstrap". ' = method %in% c("LR", "permutation", "bootstrap", "LR_sigma2"))

  if(is.null(TOJR)){
    overlap.g = x$rhythm.joint$gname[x$rhythm.joint$TOJR!="arrhy"]
  }else{
    stopifnot('The input number of types of joint rhythmicity does not match that of overlapping genes in two groups ' =
                length(x$rhythm.joint$gname)==length(TOJR))
    overlap.g = x$gname_overlap[TOJR != "arrhy"]
  }

  x1 = x$x1
  x2 = x$x2
  x1.overlap = x1$data[match(overlap.g, x1$gname), ]
  x2.overlap = x2$data[match(overlap.g, x2$gname), ]
  t1 = x1$time
  t2 = x2$time
  x1.rhythm = x1$rhythm[match(overlap.g, x1$gname),]
  x2.rhythm = x2$rhythm[match(overlap.g, x2$gname),]
  stopifnot("x$x1$P is not equal to x$x2$P. " = x$x1$P==x$x2$P)
  period = x$x1$P

  x.list = lapply(seq_along(overlap.g), function(a){
    list(x1.time = t1,
         x2.time = t2,
         y1 = as.numeric(x1.overlap[a, ]),
         y2 = as.numeric(x2.overlap[a, ]))
  })

  if(method == "LR"){
    res.list = parallel::mclapply(seq_along(x.list), function(a){
      one.res = LR_deltaR2(x.list[[a]]$x1.time, x.list[[a]]$y1,
                                           x.list[[a]]$x2.time, x.list[[a]]$y2, period,
                                           FN = TRUE)
      one.res.tab = data.frame(gname = overlap.g[a], R2.1 = x1.rhythm$R2[a], R2.2 = x2.rhythm$R2[a],
                               delta.R2 =  x2.rhythm$R2[a]-x1.rhythm$R2[a], p.R2 = one.res)
      return(one.res.tab)
    }, mc.cores = parallel.ncores)
    diffR2.tab = do.call(rbind.data.frame, res.list)
  }else if(method == "LR_sigma2"){

    res.list = parallel::mclapply(seq_along(x.list), function(a){
      one.res = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1,
                                       x.list[[a]]$x2.time, x.list[[a]]$y2, period,
                                       FN = TRUE, type = "fit")
      # one.res.tab = data.frame(gname = overlap.g[a], R2.1 = x1.rhythm$R2[a], R2.2 = x2.rhythm$R2[a],
      #                          delta.R2 =  one.res[[2]]-one.res[[1]], p.R2 = one.res$pvalue)
      one.res.tab = data.frame(gname = overlap.g[a], R2.1 = one.res[[1]], R2.2 = one.res[[2]],
                               delta.R2 =  one.res[[2]]-one.res[[1]], p.R2 = one.res$pvalue)
      return(one.res.tab)
    }, mc.cores = parallel.ncores)
    diffR2.tab = do.call(rbind.data.frame, res.list)

  }else if(method == "permutation"){
    if(!is.null(Sampling.save)){
      if(!dir.exists(Sampling.save)){
        dir.create(file.path(Sampling.save), recursive = TRUE)
        a.message = paste0("Directory created. Permutation results will be saved in ", Sampling.save)
        message(a.message)
      }
    }
    res.tab = diff_rhythmicity_permutation(x1.overlap,x2.overlap,t1,t2,overlap.g, period, x1.rhythm, x2.rhythm, nSampling,
                                           Sampling.save, Sampling.file.label, parallel.ncores, alpha)
    diffR2.tab = res.tab[, c("gname", "delta.R2", "p.R2")]
  }else if(method == "bootstrap"){
    if(!is.null(Sampling.save)){
      if(!dir.exists(Sampling.save)){
        dir.create(file.path(Sampling.save), recursive = TRUE)
        a.message = paste0("Directory created. Bootstrap results will be saved in ", Sampling.save)
        message(a.message)
      }
    }
    res.tab = diff_rhythmicity_bootstrap(x1.overlap,x2.overlap,t1,t2,overlap.g, period, x1.rhythm, x2.rhythm, nSampling,
                                         Sampling.save, Sampling.file.label, parallel.ncores,alpha)
    diffR2.tab = res.tab[, c("gname", "delta.R2", "p.R2")]
  }

  diffR2.tab$q.R2 = stats::p.adjust(diffR2.tab$p.R2, p.adjust.method)
  return(diffR2.tab)
}

# Differential rhythmicity functions --------------------------------------
##' Likelihood method to obtain p-value for differential R2.
##'
##' Likelihood method to obtain p-value for differential R2.
##' @title Likelihood method to obtain p-value for differential R2.
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param FN Correct for finite sample.
##' @return P-value for delta R2.
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' @author Caleb (copied from Caleb's github and fixed a bug)
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24)
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' LR_deltaR2(tt1, yy1, tt2, yy2)
LR_deltaR2 <- function(tt1, yy1, tt2, yy2, period = 24, FN=TRUE){

  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  # fit1 <- fitSinCurve(tt1, yy1, period = 24)
  # fit2 <- fitSinCurve(tt2, yy2, period = 24)
  fit1 <- fitSinCurve(tt1, yy1, period)
  fit2 <- fitSinCurve(tt2, yy2, period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  E1 <- A1 * cos(w * phase1)
  F1 <- A1 * sin(w * phase1)

  E2 <- A2 * cos(w * phase2)
  F2 <- A2 * sin(w * phase2)

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(E1, F1, basal1, theta1)
  p2 <- c(E2, F2, basal2, theta2)

  x_Ha <- c(p1, p2)

  asin1 <- sin(w * tt1)
  acos1 <- cos(w * tt1)
  asin2 <- sin(w * tt2)
  acos2 <- cos(w * tt2)

  eval_f_list <- function(x,asin1,acos1,asin2,acos2) {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_E1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_F1 <- - theta1 * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_E2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_F2 <- - theta2 * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2


    return( list( "objective" = -sum(ll1) - sum(ll2),
                  "gradient"  = c(partial_E1, partial_F1, partial_C1, partial_theta1,
                                  partial_E2, partial_F2, partial_C2, partial_theta2)
    )
    )
  }

  # Equality constraints
  eval_g_eq <- function(x,asin1,acos1,asin2,acos2)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    #basel1 <- p1[3]
    theta1 <- p1[4]
    #yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    #basel2 <- p2[3]
    theta2 <- p2[4]
    #yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    A2_1 <- (E1^2 + F1^2)
    A2_2 <- (E2^2 + F2^2)
    A2_1 * theta1 - A2_2 * theta2
  }

  # Equality constraints
  eval_g_eq_jac <- function(x,asin1,acos1,asin2,acos2)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    E1 <- p1[1]
    F1 <- p1[2]
    #basel1 <- p1[3]
    theta1 <- p1[4]
    #yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    #basel2 <- p2[3]
    theta2 <- p2[4]
    #yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    A2_1 <- (E1^2 + F1^2)
    A2_2 <- (E2^2 + F2^2)
    A2_1 * theta1 - A2_2 * theta2

    c(theta1 * 2 * E1, theta1 * 2 * F1, 0, A2_1,
      - theta2 * 2 * E2, - theta2 * 2 * F2, 0, - A2_2)
  }


  # Lower and upper bounds
  lb <- c(-Inf,-Inf,-Inf,0, -Inf, -Inf,-Inf, 0)
  ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
  #initial values

  ## Error in is.nloptr(ret) :
  #  If you want to use equality constraints, then you should use one of these algorithms NLOPT_LD_AUGLAG, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_AUGLAG_EQ, NLOPT_GN_ISRES, NLOPT_LD_SLSQP

  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  "local_opts" = local_opts
  opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0
                #"check_derivatives"=TRUE
  )

  res <- nloptr::nloptr ( x0 = x_Ha,
                  eval_f = eval_f_list,
                  #eval_grad_f=eval_g,
                  lb = lb,
                  ub = ub,
                  #eval_g_ineq = eval_g_ineq,
                  eval_g_eq = eval_g_eq,
                  eval_jac_g_eq = eval_g_eq_jac,
                  opts = opts,
                  asin1=asin1,
                  acos1=acos1,
                  asin2=asin2,
                  acos2=acos2)

  #
  #x_Ha
  x_H0 <- res$solution

  l0 <- - eval_f_list(x_H0,asin1,acos1,asin2,acos2)$objective
  la <- - eval_f_list(x_Ha,asin1,acos1,asin2,acos2)$objective

  LR_stat <- -2*(l0-la)

  dfdiff <- 1
  if(!FN){
    pvalue <- stats::pchisq(LR_stat,dfdiff,lower.tail = FALSE)
  } else if(FN){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- stats::pf(Fstat,df1 = r, df2 = n-k, lower.tail = FALSE)
  } else{
    stop("FN has to be TRUE or FALSE")
  }
  pvalue
}

##' Fit sin function
##'
##' Fit a sine curve where tt is time, and yy is expression value.
##' @title Fit Data Based on Sine Curve
##' @param tt Time vector.
##' @param yy Expression vector.
##' @param period Period of the sine curve. Default is 24.
##' @param parStart Initial value for optimization purpose.
##' @return A list of amp, phase, offset, peak, A, B, SST, SSE, R2.
##' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A * sin(2\pi/period \times tt) + B * cos(2*\pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1.}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period).}
##' \item{offset}{Basal level (vertical shift) based on formula 1 or on formula 2.}
##' \item{A}{A based on formula 2.}
##' \item{B}{B based on formula 2.}
##' \item{tss}{Total sum of square.}
##' \item{rss}{Residual sum of square, SSE/n is the MLE of the variance sigma square.}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss.}
##' @author Caleb (copied from Caleb's github)
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24)
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' fitSinCurve(tt, yy)

fitSinCurve <- function(tt, yy, period = 24, parStart = list(amp=3,phase=0, offset=0)){

  getPred <- function(parS, tt) {
    parS$amp * sin(2*pi/period * (tt + parS$phase)) + parS$offset
  }

  residFun <- function(p, yy, tt) yy - getPred(p,tt)

  nls.out <- minpack.lm::nls.lm(par=parStart, fn = residFun, yy = yy,	tt = tt)

  apar <- nls.out$par

  amp0 <- apar$amp
  asign <- sign(amp0)
  ## restrict amp > 0
  amp <- amp0 * asign

  phase0 <- apar$phase
  #phase <- (round(apar$phase) + ifelse(asign==1,0,12)) %% period
  phase <- (phase0 + ifelse(asign==1,0,period/2)) %% period
  offset <- apar$offset

  peak <- (period/2 * sign(amp0) - period/4 - phase) %%period
  if(peak > period/4*3) peak = peak - period

  A <- amp0 * cos(2*pi/period * phase0)
  B <- amp0 * sin(2*pi/period * phase0)

  rss <- sum(nls.out$fvec^2)
  tss <- sum((yy - mean(yy))^2)
  R2 <- 1 - rss/tss

  if(FALSE){
    amp <- apar$amp
    phase <- apar$phase
    offset <- apar$offset
  }

  res <- list(amp=amp, phase=phase, offset=offset, peak=peak, A=A, B=B, tss=tss, rss=rss, R2=R2)
  res
}

diff_rhythmicity_permutation = function(x1.overlap,x2.overlap,t1,t2,gname, period, x1.rhythm, x2.rhythm, nPermutation = 1000,
                                        permutation.save = NULL, permutation.file.label = "gII_vs_gI", parallel.ncores = 1, alpha){

  x12.overlap = cbind(x1.overlap, x2.overlap)
  t.all = c(t1, t2)

  indexes = seq_len(length(t.all))
  n1 = length(t1)

  perm.list = parallel::mclapply(seq_along(nPermutation), function(b){
    # set.seed(b)
    index1 = sample(indexes, n1)
    index2 = setdiff(indexes, index1)
    x1.b = list(data = x12.overlap[, index1],
                time = t.all[index1],
                gname = gname)
    x2.b = list(data = x12.overlap[, index2],
                time = t.all[index2],
                gname = gname)
    perm.b = DCP_Rhythmicity(x1.b, x2.b, method = "Sidak_FS", period, alpha, CI = FALSE)
    x1.rhythm.b = perm.b$x1$rhythm
    x2.rhythm.b = perm.b$x2$rhythm
    diffparas_null = list(M_null = x2.rhythm.b$M-x1.rhythm.b$M,
                          A_null = x2.rhythm.b$A-x1.rhythm.b$A,
                          phase_null = choose.abs.min(cbind(x2.rhythm.b$phase-x1.rhythm.b$phase,
                                                            x2.rhythm.b$phase+2*pi-x1.rhythm.b$phase,
                                                            x2.rhythm.b$phase-2*pi-x1.rhythm.b$phase)),
                          peak_null = choose.abs.min(cbind(x2.rhythm.b$peak-x1.rhythm.b$peak,
                                                           x2.rhythm.b$peak+period-x1.rhythm.b$peak,
                                                           x2.rhythm.b$peak-period-x1.rhythm.b$peak)),
                          R2_null = x2.rhythm.b$R2-x1.rhythm.b$R2)
    if(!is.null(permutation.save)){
      save(x1.rhythm.b, x2.rhythm.b, index1, index2, gname,
           file = paste0(file.path(permutation.save, paste0(permutation.file.label, "_PermBetween", b, ".rData"))))
    }

    if(b%%(nPermutation/100)==0){
      cat(paste0("Differential Rhythmicity Analysis permutation step: ", round(b/nPermutation, 2)*100, "% done. \n"))
    }

    return(diffparas_null)

  }, mc.cores = parallel.ncores)

  M_null = do.call(cbind, lapply(perm.list, `[[`, "M_null"))
  A_null = do.call(cbind, lapply(perm.list, `[[`, "A_null"))
  phase_null = do.call(cbind, lapply(perm.list, `[[`, "phase_null"))
  peak_null = do.call(cbind, lapply(perm.list, `[[`, "peak_null"))
  R2_null = do.call(cbind, lapply(perm.list, `[[`, "R2_null"))

  x1.true = x1.rhythm
  x2.true = x2.rhythm

  M_obs = x2.true$M-x1.true$M
  A_obs = x2.true$A-x1.true$A
  phase_obs = choose.abs.min(cbind(x2.true$phase-x1.true$phase, x2.true$phase+2*pi-x1.true$phase, x2.true$phase-2*pi-x1.true$phase))
  peak_obs = choose.abs.min(cbind(x2.true$peak-x1.true$peak, x2.true$peak+period-x1.true$peak, x2.true$peak-period-x1.true$peak))
  R2_obs = x2.true$R2-x1.true$R2

  #ap_R2_perGene <- apply(R2_null - R2_obs,1,function(x) min(mean(x >= 0), mean(x <= 0)) * 2)

  ## permutation p-value by pooling all genes together
  ap_M_allGene0 <- rank(cbind(M_obs, M_null))[seq_along(gname)]/(length(gname)*(nPermutation+ 1))
  ap_M_allGene <- pmin(ap_M_allGene0, 1 - ap_M_allGene0) * 2

  ap_A_allGene0 <- rank(cbind(A_obs, A_null))[seq_along(gname)]/(length(gname)*(nPermutation+ 1))
  ap_A_allGene <- pmin(ap_A_allGene0, 1 - ap_A_allGene0) * 2

  ap_phase_allGene0 <- rank(cbind(phase_obs, phase_null))[seq_along(gname)]/(length(gname)*(nPermutation+ 1))
  ap_phase_allGene <- pmin(ap_phase_allGene0, 1 - ap_phase_allGene0) * 2

  ap_peak_allGene0 <- rank(cbind(peak_obs, peak_null))[seq_along(gname)]/(length(gname)*(nPermutation+ 1))
  ap_peak_allGene <- pmin(ap_peak_allGene0, 1 - ap_peak_allGene0) * 2

  ap_R2_allGene0 <- rank(cbind(R2_obs, R2_null))[seq_along(gname)]/(length(gname)*(nPermutation+ 1))
  ap_R2_allGene <- pmin(ap_R2_allGene0, 1 - ap_R2_allGene0) * 2

  out = data.frame(gname = gname,
                   delta.M = M_obs,
                   p.M = ap_M_allGene,
                   delta.A = A_obs,
                   p.A = ap_A_allGene,
                   delta.phase = phase_obs,
                   delta.peak = peak_obs,
                   p.phase = ap_phase_allGene,
                   #p.peak = ap_peak_allGene, #p for permutaiton p for peak and phase are almost the same, no need for both
                   delta.R2 = R2_obs,
                   p.R2 = ap_R2_allGene)
  if(!is.null(permutation.save)){
    utils::write.csv(out, paste0(file.path(permutation.save, paste0(permutation.file.label, "_PermBetween", "_out", ".csv"))))
  }
  return(out)
}
diff_rhythmicity_bootstrap = function(x1.overlap,x2.overlap,t1,t2,gname, period, x1.rhythm, x2.rhythm, nSampling = 1000,
                                      Sampling.save = NULL, Sampling.file.label = "gII_vs_gI", parallel.ncores = 1, alpha){

  #indexes = seq_len(length(t.all))
  index1.0 = seq_len(length(t1))
  index2.0 = seq_len(length(t2))

  boot.list = parallel::mclapply(seq_along(nSampling), function(b){
    # set.seed(b)
    index1 = sample(index1.0, replace = TRUE)
    index2 = sample(index2.0, replace = TRUE)
    x1.b = list(data = x1.overlap[, index1],
                time = t1[index1],
                gname = gname)
    x2.b = list(data = x2.overlap[, index2],
                time = t2[index2],
                gname = gname)
    boot.b = DCP_Rhythmicity(x1.b, x2.b, method = "Sidak_FS", period, alpha, CI = FALSE)
    x1.rhythm.b = boot.b$x1$rhythm
    x2.rhythm.b = boot.b$x2$rhythm
    diffparas_null = list(M_null = x2.rhythm.b$M-x1.rhythm.b$M,
                          A_null = x2.rhythm.b$A-x1.rhythm.b$A,
                          phase_null = choose.abs.min(cbind(x2.rhythm.b$phase-x1.rhythm.b$phase,
                                                            x2.rhythm.b$phase+2*pi-x1.rhythm.b$phase,
                                                            x2.rhythm.b$phase-2*pi-x1.rhythm.b$phase)),
                          peak_null = choose.abs.min(cbind(x2.rhythm.b$peak-x1.rhythm.b$peak,
                                                           x2.rhythm.b$peak+period-x1.rhythm.b$peak,
                                                           x2.rhythm.b$peak-period-x1.rhythm.b$peak)),
                          R2_null = x2.rhythm.b$R2-x1.rhythm.b$R2)
    if(!is.null(Sampling.save)){
      save(x1.rhythm.b, x2.rhythm.b, index1, index2, gname,
           file = paste0(file.path(Sampling.save, paste0(Sampling.file.label, "_Bootstap", b, ".rData"))))
    }

    if(b%%(nSampling/100)==0){
      cat(paste0("Differential Rhythmicity Analysis Bootstrap step: ", round(b/nSampling, 2)*100, "% done. \n"))
    }

    return(diffparas_null)

  }, mc.cores = parallel.ncores)


  M_null = do.call(cbind, lapply(boot.list, `[[`, "M_null")); M_null_sd = apply(M_null, 1, stats::sd)
  A_null = do.call(cbind, lapply(boot.list, `[[`, "A_null")); A_null_sd = apply(A_null, 1, stats::sd)
  phase_null = do.call(cbind, lapply(boot.list, `[[`, "phase_null")); phase_null_sd = apply(phase_null, 1, function(a){
    a.max = max(a, na.rm = TRUE)
    a[a<(a.max-2*pi)] = a[a<(a.max-2*pi)]+2*pi
    stats::sd(a)
  })
  #phase_null_sd = apply(phase_null, 1, stats::sd)
  peak_null = do.call(cbind, lapply(boot.list, `[[`, "peak_null")); peak_null_sd = apply(peak_null, 1, function(a){
    a.max = max(a, na.rm = TRUE)
    a[a<(a.max-period)] = a[a<(a.max-period)]+period
    stats::sd(a)
  })
  # peak_null_sd = apply(peak_null, 1, stats::sd)
  R2_null = do.call(cbind, lapply(boot.list, `[[`, "R2_null")); R2_null_sd = apply(R2_null, 1, stats::sd)

  x1.true = x1.rhythm
  x2.true = x2.rhythm

  M_obs = x2.true$M-x1.true$M
  A_obs = x2.true$A-x1.true$A
  phase_obs = choose.abs.min(cbind(x2.true$phase-x1.true$phase, x2.true$phase+2*pi-x1.true$phase, x2.true$phase-2*pi-x1.true$phase))
  peak_obs = choose.abs.min(cbind(x2.true$peak-x1.true$peak, x2.true$peak+period-x1.true$peak, x2.true$peak-period-x1.true$peak))
  R2_obs = x2.true$R2-x1.true$R2

  #ap_R2_perGene <- apply(R2_null - R2_obs,1,function(x) min(mean(x >= 0), mean(x <= 0)) * 2)

  ## Sampling p-value by pooling all genes together
  ap_M_allGene0 <- stats::pnorm(M_obs/M_null_sd)
  ap_M_allGene <- pmin(ap_M_allGene0, 1 - ap_M_allGene0) * 2

  ap_A_allGene0 <- stats::pnorm(A_obs/A_null_sd)
  ap_A_allGene <- pmin(ap_A_allGene0, 1 - ap_A_allGene0) * 2

  ap_phase_allGene0 <- stats::pnorm(phase_obs/phase_null_sd)
  ap_phase_allGene <- pmin(ap_phase_allGene0, 1 - ap_phase_allGene0) * 2

  ap_peak_allGene0 <- stats::pnorm(peak_obs/peak_null_sd)
  ap_peak_allGene <- pmin(ap_peak_allGene0, 1 - ap_peak_allGene0) * 2

  ap_R2_allGene0 <- stats::pnorm(R2_obs/R2_null_sd)
  ap_R2_allGene <- pmin(ap_R2_allGene0, 1 - ap_R2_allGene0) * 2

  out = data.frame(gname = gname,
                   delta.M = M_obs,
                   p.M = ap_M_allGene,
                   delta.A = A_obs,
                   p.A = ap_A_allGene,
                   delta.phase = phase_obs,
                   delta.peak = peak_obs,
                   p.phase = ap_phase_allGene,
                   #p.peak = ap_peak_allGene, #p for bootutaiton p for peak and phase are almost the same, no need for both
                   delta.R2 = R2_obs,
                   p.R2 = ap_R2_allGene)
  if(!is.null(Sampling.save)){
    utils::write.csv(out, paste0(file.path(Sampling.save, paste0(Sampling.file.label, "_Bootstrap", "_out", ".csv"))))
  }
  return(out)
}

choose.abs.min = function(mat = matrix(stats::rnorm(21), ncol = 3)){
  abs.min.vec = apply(mat, 1, function(a){
    min.abs.idx = which.min(abs(a))
    return(a[min.abs.idx])
  })
  return(abs.min.vec)
}

