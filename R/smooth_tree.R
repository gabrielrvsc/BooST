#' Estimate BooST
#'
#' Estimates Boosting of Smooth Trees (BooST)
#'
#'
#' @param x Design matrix with explanatory variables.
#' @param y Response variable.
#' @param p Proportion of variables tested in each node split (default 1).
#' @param d_max Number of splits in each tree (default 4).
#' @param gamma Transiction function intensity. Bigger numbers makes the transition less smoth. The default is a sequence of values (0.5:5) to be randomized in each new node. Multiple values may be supplied in a vector to increase the model randomness.
#'
#' @return An object with S3 class "SmoothTree".
#' \item{Model}{A list with all trees.}
#' \item{fitted.values}{Final model fitted values.}
#' \item{nvar}{Number of variables in x.}
#' \item{varnames}{colnames of x to be used in other functions.}
#' \item{call}{The matched call.}
#' @keywords BooST, Boosting, Smooth Tree, Partial Effects
#' @export
#' @examples
#' ## == to be made == ##
#'
#' @references
#' blablabla
#'
# @seealso \code{\link{BooST}}, \code{\link{predict.SmoothTree}},  \code{\link{estimate_derivative}}


smooth_tree=function(x, y, p = 1, d_max = 4, gamma = seq(0.5,5,0.01),node_obs=nrow(x)/200){
  tree=grow_tree(x,y, p = p, d_max = d_max, gamma = gamma, node_obs=node_obs)
  fitted.values=tree$fitted.values
  result=list(tree=tree$tree, fitted.values=fitted.values, nvar = ncol(x) , varnames=colnames(x) ,call=match.call())
  class(result)="SmoothTree"
  return(result)
}
