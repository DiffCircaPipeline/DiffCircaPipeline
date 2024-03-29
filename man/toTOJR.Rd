% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCP_Rhythmicity.R
\name{toTOJR}
\alias{toTOJR}
\title{Title Types of Joint Rhythmicity (TOJR)}
\usage{
toTOJR(
  x,
  method = "Sidak_FS",
  amp.cutoff = 0,
  alpha = 0.05,
  adjustP = TRUE,
  p.adjust.method = "BH",
  parallel.ncores = 1
)
}
\arguments{
\item{x}{one of the following two: \itemize{
\item output of DCP_rhythmicity(x1, x2), both x1 and x2 are not NULL. (A list with the rhythm.joint component).
\item A list of with two outputs from DCP_rhythmicity(x1, x2 = NULL), each using data from a different group. }}

\item{method}{character string specifying the algorithm used for joint rhythmicity categorization. Should be one of "Sidak_FS", "Sidak_BS", "VDA", "AWFisher".}

\item{amp.cutoff}{Only genes with amplitude greater than amp.cutoff are consirdered rhythmic}

\item{alpha}{integer. Threshold for rhythmicity p-value in joint rhythmicity categorization.}

\item{adjustP}{logic. Should joint rhythmicity categorization be based on adjusted p-value?}

\item{p.adjust.method}{input for p.adjust() in R package \code{stat}}

\item{parallel.ncores}{integer. Number of cores used if using parallel computing with \code{mclapply()}. Not functional for windows system.}
}
\value{
Joint rhythmicity categories of genes
}
\description{
Categorize genes to four types of joint rhythmicity (TOJR).
}
\examples{
#Re-calculate TOJR for DCP_rhythmicity(x1, x2) output with q-value cutoff 0.1
x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
TOJR.new = toTOJR(rhythm.res, alpha = 0.1, adjustP = TRUE)

#Calculate TOJR for two DCP_rhythmicity(x1, x2 = NULL) outputs with p-value cutoff 0.05
x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
rhythm.res1 = DCP_Rhythmicity(x1 = x[[1]])
rhythm.res2 = DCP_Rhythmicity(x1 = x[[2]])
TOJR = toTOJR(x = list(x1 = rhythm.res1, x2 = rhythm.res2),
alpha = 0.05, adjustP = FALSE)
}
