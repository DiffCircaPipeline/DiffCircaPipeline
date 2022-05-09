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

  x.list = lapply(1:length(overlap.g), function(a){
    list(x1.time = t1,
         x2.time = t2,
         y1 = as.numeric(x1.overlap[a, ]),
         y2 = as.numeric(x2.overlap[a, ]))
  })

  if(method == "LR"){
    res.list = parallel::mclapply(1:length(x.list), function(a){
      one.res = differentialR2::LR_deltaR2(x.list[[a]]$x1.time, x.list[[a]]$y1,
                                           x.list[[a]]$x2.time, x.list[[a]]$y2,
                                           FN = TRUE)
      one.res.tab = data.frame(gname = overlap.g[a], R2.1 = x1.rhythm$R2[a], R2.2 = x2.rhythm$R2[a],
                               delta.R2 =  x2.rhythm$R2[a]-x1.rhythm$R2[a], p.R2 = one.res)
      return(one.res.tab)
    }, mc.cores = parallel.ncores)
    diffR2.tab = do.call(rbind.data.frame, res.list)
  }else if(method == "LR_sigma2"){

    res.list = parallel::mclapply(1:length(x.list), function(a){
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
        message(paste0("Directory created. Permutation results will be saved in ", Sampling.save))
      }
    }
    res.tab = diff_rhythmicity_permutation(x1.overlap,x2.overlap,t1,t2,overlap.g, period, x1.rhythm, x2.rhythm, nSampling,
                                           Sampling.save, Sampling.file.label, parallel.ncores, alpha)
    diffR2.tab = res.tab[, c("gname", "delta.R2", "p.R2")]
  }else if(method == "bootstrap"){
    if(!is.null(Sampling.save)){
      if(!dir.exists(Sampling.save)){
        dir.create(file.path(Sampling.save), recursive = TRUE)
        message(paste0("Directory created. Bootstrap results will be saved in ", Sampling.save))
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

diff_rhythmicity_permutation = function(x1.overlap,x2.overlap,t1,t2,gname, period, x1.rhythm, x2.rhythm, nPermutation = 1000,
                                        permutation.save = NULL, permutation.file.label = "gII_vs_gI", parallel.ncores = 1, alpha){

  x12.overlap = cbind(x1.overlap, x2.overlap)
  t.all = c(t1, t2)

  indexes = seq_len(length(t.all))
  n1 = length(t1)

  perm.list = parallel::mclapply(1:nPermutation, function(b){
    set.seed(b)
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
  ap_M_allGene0 <- rank(cbind(M_obs, M_null))[1:length(gname)]/(length(gname)*(nPermutation+ 1))
  ap_M_allGene <- pmin(ap_M_allGene0, 1 - ap_M_allGene0) * 2

  ap_A_allGene0 <- rank(cbind(A_obs, A_null))[1:length(gname)]/(length(gname)*(nPermutation+ 1))
  ap_A_allGene <- pmin(ap_A_allGene0, 1 - ap_A_allGene0) * 2

  ap_phase_allGene0 <- rank(cbind(phase_obs, phase_null))[1:length(gname)]/(length(gname)*(nPermutation+ 1))
  ap_phase_allGene <- pmin(ap_phase_allGene0, 1 - ap_phase_allGene0) * 2

  ap_peak_allGene0 <- rank(cbind(peak_obs, peak_null))[1:length(gname)]/(length(gname)*(nPermutation+ 1))
  ap_peak_allGene <- pmin(ap_peak_allGene0, 1 - ap_peak_allGene0) * 2

  ap_R2_allGene0 <- rank(cbind(R2_obs, R2_null))[1:length(gname)]/(length(gname)*(nPermutation+ 1))
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

  boot.list = parallel::mclapply(1:nSampling, function(b){
    set.seed(b)
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

