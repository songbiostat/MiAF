#' Associated taxa selected by combining lower- and upper-tail p-values
#'
#' @description Associated taxa selected by combining lower- and upper-tail p-values of a test
#'              using one of abundance representations.
#'
#' @param lowerlist The index of selected taxa by combining lower-tail p-values.
#' @param upperlist The index of selected taxa by combining upper-tail p-values.
#' @param com The index of selected one-sided test after combining two one-sided community-level tests.
#'            1: the one-sided test using lower-tail p-values.
#'            2: the one-sided test using upper-tail p-values.
#' @param label The lables of all non-rooted nodes.
#'
#' @return A list of selected associated taxa after combining two one-sided community-level tests:
#'         \item{lower}{Under-presented taxa based on lower-tail p-values.}
#'         \item{upper}{Over-presented taxa based on upper-tail p-values.}
#'
#' @export
#'
select.l.u <- function(lowerlist, upperlist, com, label){
  sel.lower <- c()
  sel.upper <- c()
  if(1 %in% com) sel.lower <- label[lowerlist]
  if(2 %in% com) sel.upper <- label[upperlist]
  return(list(lower = sel.lower, upper = sel.upper))
}


#' Associated taxa selected by combining p-values of four tests using different abundance representations
#'
#' @param sel.unweight The list of selected associated taxa by unweighted UniFrac-like test.
#' @param sel.weight The list of selected associated taxa by weighted UniFrac-like test.
#' @param sel.5 The list of selected associated taxa by generalized UniFrac-like test.
#' @param sel.tip The list of selected associated taxa by leaf-nodes-only test.
#' @param com The index of selected tests using different abundance representations.
#'            1: unweighted UniFrac-like test.
#'            2: weighted UniFrac-like test.
#'            3: generalized UniFrac-like test.
#'            4: leaf-nodes-only test.
#'
#' @return A list of selected associated taxa after combining four using different abundance representations:
#'         \item{lower}{Under-presented taxa based on lower-tail p-values.}
#'         \item{upper}{Over-presented taxa based on upper-tail p-values.}
#'
#' @export
#'
select.com <- function(sel.unweight, sel.weight, sel.5, sel.tip, com){
  sel.lower <- c()
  sel.upper <- c()
  if(1 %in% com){
    sel.lower <- sel.unweight$lower
    sel.upper <- sel.unweight$upper
  }
  if(2 %in% com){
    sel.lower <- c(sel.lower, sel.weight$lower)
    sel.upper <- c(sel.upper, sel.weight$upper)
  }
  if(3 %in% com){
    sel.lower <- c(sel.lower, sel.5$lower)
    sel.upper <- c(sel.upper, sel.5$upper)
  }
  if(4 %in% com){
    sel.lower <- c(sel.lower, sel.tip$lower)
    sel.upper <- c(sel.upper, sel.tip$upper)
  }
  sel.lower <- unique(sel.lower)
  sel.upper <- unique(sel.upper)
  return(list(lower = sel.lower, upper = sel.upper))
}
