#' Summary of Differential Rhythmicity Tests
#' Summarize differential rhythmicity Tests results into a table.
#' @param result a list of vectors, usually p-values. Each vector will be summarized with a different cutoff.
#' @param test a vector of strings indicating the name of tests.
#' @param type a vector of strings indicating the type of the results.
#' @param val The cutoff used to dichotomize the result vectors.
#' @param out Should the output table be "long" or "wide"?
#'
#' @return a table of summarized results
#' @export
#'
#' @examples
#' x = DCP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#' rhythm.res = DCP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' rhythm.diffR2 = DCP_DiffPar(rhythm.res)
#' SummarizeDR(result = list(rhythm.diffR2$p.R2, rhythm.diffR2$p.R2, rhythm.diffR2$q.R2),
#' test = c(rep("DRF", 2), "DRF"),
#' type = c(rep("p-value", 2), "q-value"),
#' val = c(0.1, 0.05, 0.4), "long")
SummarizeDR = function(result = DCP_dR2$p.R2, test = "DRF", type = "p-value", val = 0.05, out = "long"){
  if(class(result)!="list"){
    result = list(result)
  }
  l = length(result)
  stopifnot("Arguments `result`, `test`, `type`, `val` must have the same length" = (length(test)==l&length(type)==l&length(val)==l))
  tab = lapply(1:l, function(a){
    a.result = result[[a]]
    a.test = test[a]
    a.type = type[a]
    a.val = val[a]
    a.tab0 = table(a.result<a.val)
    a.tab = data.frame(nTRUE = unname(ifelse(is.na(a.tab0["TRUE"]), 0, a.tab0["TRUE"])),
                       nFALSE = unname(ifelse(is.na(a.tab0["FALSE"]), 0, a.tab0["FALSE"])),
                       test = a.test,
                       cutoff = paste0(a.type, "<", a.val))
    return(a.tab)
  })
  tab = do.call(rbind.data.frame, tab)
  tab$nTRUE = paste0(tab$nTRUE, "(/", tab$nTRUE+tab$nFALSE, ")")
  if(out=="long"){
    return(tab[, c(1, 3,4)])
  }else if(out =="wide"){
    tab.list = lapply(1:nrow(tab[, c(1, 3,4)]), function(i){
      a.cell = data.frame(x = c(tab[i, "nTRUE"]))
      rownames(a.cell) = NULL
      colnames(a.cell) = paste0(tab[i, "test"], " ", tab[i, "cutoff"])
      return(a.cell)
    })
    tab2 = do.call(cbind.data.frame,
                   tab.list)
    return(tab2)
  }
}
