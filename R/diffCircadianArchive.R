LR_diff <- function(tt1,yy1,tt2,yy2,period=24,method="LR",FN=TRUE,type="all"){
  # Likelihood-based tests for differential circadian pattern detection.
  # LR
  if (method=="LR"){
    # LR + all
    if(type=="all"){
      list(LRTest_diff_amp(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
           LRTest_diff_phase(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
           LRTest_diff_offset(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
           LRTest_diff_sigma2(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN))
    }
    # LR + amp
    else if(type=="amplitude"){
      LRTest_diff_amp(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    else if(type=="phase"){
      LRTest_diff_phase(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    else if(type=="basal"){
      LRTest_diff_offset(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    else if(type=="fit"){
      LRTest_diff_sigma2(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    else{
      stop("Please check your input! type = 'all','amplitude','phase','offset' or 'fit' and test = 'LR' or 'Wald'")
    }
  }





  # # Wald
  # else if (method=="Wald"){
  #   # Wald + all
  #   if(type=="all"){
  #     list(WaldTest_diff_amp(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
  #          WaldTest_diff_phase(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
  #          WaldTest_diff_offset(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
  #          WaldTest_diff_sigma2(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN))
  #   }
  #   # Wald + amp
  #   else if(type=="amplitude"){
  #     WaldTest_diff_amp(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
  #   }
  #   # Wald + phase
  #   else if(type=="phase"){
  #     WaldTest_diff_phase(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
  #   }
  #   # Wald + offset
  #   else if(type=="basal"){
  #     WaldTest_diff_offset(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
  #   }
  #   # Wald + sigma2
  #   else if(type=="fit"){
  #     WaldTest_diff_sigma2(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
  #   }
  #   else("Please check your input! type = 'all','amplitude','phase','offset' or 'rhythmicity' and test = 'LR' or 'Wald'")
  # }
  # else("Please check your input! type = 'all','amplitude','phase','basal' or 'fit' and test = 'LR' or 'Wald'")

}

LRTest_diff_amp<- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){

  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

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
    A2_1 - A2_2
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

    c(2 * E1, 2 * F1, 0, 0,
      - 2 * E2, - 2 * F2, 0, 0)
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

  amp_c <- sqrt(x_H0[1]^2 + x_H0[2]^2)
  amp_c2 <- sqrt(x_H0[5]^2 + x_H0[6]^2)

  res <- list(amp_1=A1, amp_2=A2, amp_c=amp_c,
              l0=l0,
              la=la,
              #df = dfdiff,
              stat=-2*(l0-la),
              pvalue=pvalue)
  return(res)
}

LRTest_diff_offset <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

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
    basel1 <- p1[3]
    theta1 <- p1[4]
    #yhat1 <- E1 * asin1 + F1 * acos1 + basel1

    E2 <- p2[1]
    F2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    #yhat2 <- E2 * asin2 + F2 * acos2 + basel2

    basel1 - basel2
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

    c(0, 0, 1, 0,
      0, 0, -1, 0)
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

  offset_c <- x_H0[3]
  offset_c2 <- x_H0[7]


  res <- list(offset_1=basal1, offset_2=basal2, offset_c=offset_c,
              l0=l0,
              la=la,
              #df = dfdiff,
              stat=-2*(l0-la),
              pvalue=pvalue)
  return(res)
}

LRTest_diff_phase <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){

  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  #period <- 24
  w <- 2*pi/period

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)

  A1 <- fit1$amp
  A2 <- fit2$amp

  phase1 <- fit1$phase
  phase2 <- fit2$phase

  if(phase2 - phase1 > period/2){
    phase2 <- phase2 - period
  } else if(phase1 - phase2 > period/2){
    phase1 <- phase1 - period
  }

  basal1 <- fit1$offset
  basal2 <- fit2$offset

  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2

  p1 <- c(A1, phase1, basal1, theta1)
  p2 <- c(A2, phase2, basal2, theta2)

  x_Ha <- c(p1, p2)

  eval_f_list <- function(x) {
    p1 <- x[1:4]
    p2 <- x[5:8]

    A1 <- p1[1]
    phase1 <- p1[2]
    basel1 <- p1[3]
    theta1 <- p1[4]
    asin1 <- sin(w * (tt1 + phase1) )
    acos1 <- cos(w * (tt1 + phase1) )
    yhat1 <- A1 * asin1 + basel1

    A2 <- p2[1]
    phase2 <- p2[2]
    basel2 <- p2[3]
    theta2 <- p2[4]
    asin2 <- sin(w * (tt2 + phase2) )
    acos2 <- cos(w * (tt2 + phase2) )
    yhat2 <- A2 * asin2 + basel2

    ll1_a <- log(theta1)/2
    ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
    ll1 <- ll1_a - ll1_b

    ll2_a <- log(theta2)/2
    ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
    ll2 <- ll2_a - ll2_b

    partial_A1 <- - theta1 * sum((yy1 - yhat1) * asin1)
    partial_phase1 <- - theta1 * A1 * w * sum((yy1 - yhat1) * acos1)
    partial_C1 <- - theta1 * sum(yy1 - yhat1)
    partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

    partial_A2 <- - theta2 * sum((yy2 - yhat2) * asin2)
    partial_phase2 <- - theta2 * A2 * w * sum((yy2 - yhat2) * acos2)
    partial_C2 <- - theta2 * sum(yy2 - yhat2)
    partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2


    return( list( "objective" = -sum(ll1) - sum(ll2),
                  "gradient"  = c(partial_A1, partial_phase1, partial_C1, partial_theta1,
                                  partial_A2, partial_phase2, partial_C2, partial_theta2)
    )
    )
  }

  # Equality constraints
  eval_g_eq <- function(x)
  {
    p1 <- x[1:4]
    p2 <- x[5:8]

    phase1 <- p1[2]
    phase2 <- p2[2]

    phase1 - phase2
  }

  # Equality constraints
  eval_g_eq_jac <- function(x)
  {
    c(0, 1, 0, 0,
      0, -1, 0, 0)
  }


  # Lower and upper bounds
  lb <- c(0,-Inf,-Inf,0, 0, -Inf,-Inf, 0)
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
                  opts = opts)

  #
  #x_Ha
  x_H0 <- res$solution

  l0 <- - eval_f_list(x_H0)$objective
  la <- - eval_f_list(x_Ha)$objective

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

  phase_c <- x_H0[2]
  phase_c2 <- x_H0[6]

  res <- list(phase_1=phase1, phase_2=phase2, phase_c=phase_c,
              l0=l0,
              la=la,
              #df = dfdiff,
              stat=-2*(l0-la),
              pvalue=pvalue)
  return(res)
}

LRTest_diff_sigma2 <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  fit1 <- fitSinCurve(tt1, yy1, period = period)
  fit2 <- fitSinCurve(tt2, yy2, period = period)
  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss
  sigma2_C <- 1/(n1 + n2) * (sigma2_1 * n1 + sigma2_2 * n2)

  ## H0 (sigma2_C)
  l0_1 <- -n1/2*log(2*pi*sigma2_C) # - 1/(2*sigma2_C)*sum(residual_1^2)
  l0_2 <- -n2/2*log(2*pi*sigma2_C) # - 1/(2*sigma2_C)*sum(residual_2^2)
  l0 <- l0_1 + l0_2

  ## H1 (sigma2_1, sigma2_2)
  la_1 <- -n1/2*log(2*pi*sigma2_1) # - 1/(2*sigma2_1)*sum(residual_1^2)
  la_2 <- -n2/2*log(2*pi*sigma2_2) # - 1/(2*sigma2_2)*sum(residual_2^2)
  la <- la_1 + la_2

  dfdiff <- 1
  if(FN==FALSE){
    pvalue <- stats::pchisq(-2*(l0-la),dfdiff,lower.tail = FALSE)
  } else if(FN==TRUE){
    LR_stat <- -2*(l0-la)
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- stats::pf(Fstat,df1 = r, df2 = n-k, lower.tail = FALSE)
  }


  res <- list(sigma2_1=sigma2_1, sigma2_2=sigma2_2, sigma2_c=sigma2_C,
              l0=l0,
              la=la,
              #df = dfdiff,
              stat=-2*(l0-la),
              pvalue=pvalue)
  return(res)
}
