#' Creates analytical expression for partial effects
#'
#' Creates analytical expression for partial effects estimation of a BoosT or a SmoothTree object.
#'
#'
#' @param object A BooST or a SmoothTree object.
#' @keywords BooST, Boosting, Smooth Tree, Partial Effects, Predict
#' @export
#' @examples
#' ## == to be made == ##
#'
# @seealso \code{\link{BooST}}, \code{\link{smooth_tree}}, \code{\link{estimate_derivative}}

derivative_expression <- function(object) {

  model = object$Model
  n = length(model)
  ncol=object$nvar
  variable=1:ncol

  trees=lapply(model,tree_derivative,variable=variable)

  return('Analytical' = trees)
}
