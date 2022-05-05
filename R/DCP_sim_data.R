#' Simulating Data for Rhythmicity Comparison
#'
#' \code{DCP_sim_data} simulates a simple two-group study. The simulated data is just used for illustration purposes and is not intended to mimic real data.
#' @param ngene integer. Number of features. Shared by both groups
#' @param nsample integer. Number of samples Shared by both groups
#' @param A1 vector of length 2. c(min, max) amplitude for group I
#' @param A2 vector of length 2. c(min, max) amplitude for group II
#' @param phase1 vector of length 2. c(min, max) phase for group I (values not in the 0 to \eqn{2\pi} interval will be converted)
#' @param phase2 vector of length 2. c(min, max) phase for group II
#' @param M1 vector of length 2. c(min, max) MESOR for group I
#' @param M2 vector of length 2. c(min, max) MESOR for group I
#' @param sigma1 numeric. Noise level for group I.
#' @param sigma2 numeric. Noise level for group II.
#'
#' @return A list of length two, each is data for one group.
#' @export
#'
#' @examples
#' #simulate data with shifted phase
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
DCP_sim_data = function(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3), phase1=c(0, 2*pi), phase2=c(0, 2*pi), M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1){
  a.A1 = stats::runif(ngene, A1[1], A1[2]); a.A2 = stats::runif(ngene, A2[1], A2[2]);
  a.m1 = stats::runif(ngene, M1[1], M1[2]); a.m2 = stats::runif(ngene, M2[1], M2[2]);
  a.phase1 = stats::runif(ngene, phase1[1], phase1[2]); a.phase2 = stats::runif(ngene, phase2[1], phase2[2]);
  a.sigma1 = sigma1; a.sigma2 = sigma2
  x1.time = stats::runif(nsample, min = 0, max = 24)
  x2.time = stats::runif(nsample, min = 0, max = 24)

  noise.mat1 = matrix(stats::rnorm(ngene*nsample, 0, a.sigma1), ncol = nsample, nrow = ngene)
  signal.mat1 = t(sapply(1:ngene, function(a){a.m1[a]+a.A1[a]*cos(2*pi/24*x1.time+a.phase1[a])}))

  noise.mat2 = matrix(stats::rnorm(ngene*nsample, 0, a.sigma2), ncol = nsample, nrow = ngene)
  signal.mat2 = t(sapply(1:ngene, function(a){a.m2[a]+a.A2[a]*cos(2*pi/24*x2.time+a.phase2[a])}))

  x1.data = as.data.frame(noise.mat1 + signal.mat1)
  x2.data = as.data.frame(noise.mat2 + signal.mat2)

  x1 = list(data = x1.data,
            time = x1.time,
            gname = paste("gene", seq_len(ngene)))

  x2 = list(data = x2.data,
            time = x2.time,
            gname = paste("gene", seq_len(ngene)))
  return(list(x1, x2))
}
