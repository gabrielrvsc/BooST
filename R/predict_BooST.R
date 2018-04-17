#' Predict Method for BooST Fits
#'
#' Obtains predictions for a BooST model and a given design matrix.
#'
#'
#' @param object A BooST or a SmoothTree object.
#' @param newx Matrix of new values of x at wich predictions are to be made.
#' @param ... Additional arguments for other methods.
#' @keywords BooST, Boosting, Smooth Tree, Partial Effects, Predict
#' @export
#' @examples
#' ## == to be made == ##
#'
# @seealso \code{\link{BooST}}, \code{\link{smooth_tree}}

predict.BooST=function(object,newx=NULL,...){

  if(is.null(newx)){
    return(stats::fitted(object))
  }

  if(is.vector(newx)){newx=matrix(newx,nrow=1)}

  v=object$v
  y0=object$ybar
  rho=object$rho
  model=object$Model
  rhov=rho*v
  fitaux=t(t(Reduce("cbind",lapply(model,function(t)eval_tree(newx,t$tree))))*rhov)
  fitted.values=y0+rowSums(fitaux)
  return(fitted.values)
}
