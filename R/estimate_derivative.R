#' Creates analytical expression for partial effects
#'
#' Creates analytical expression for partial effects estimation of a BoosT or a SmoothTree object.
#'
#'
#' @param object A BooST or a SmoothTree object.
#' @param x Matrix of obsevations where the partial derivatives are desired.
#' @param de_exp Derivative expressions from the derivative_expression function.
#' @keywords BooST, Boosting, Smooth Tree, Partial Effects
#' @export
#' @examples
#' ## == to be made == ##
#'
# @seealso \code{\link{BooST}}, \code{\link{smooth_tree}}, \code{\link{derivative_expression}}

estimate_derivatives=function(object,x,de_exp){
  if(class(object)=="BooST"){
    v=object$v
    rho=object$rho
  }
  if(class(object)=="SmoothTree"){
    v=1
    rho=1
  }
  if(is.vector(x)){x=matrix(x,nrow=1)}
  partial_effects=apply(x,1,function(x){
    eval_gradient(de_exp, x ,v=v,rho=rho)
  })
  partial_effects=t(partial_effects)
  colnames(partial_effects)=object$varnames
  return(partial_effects)
}
