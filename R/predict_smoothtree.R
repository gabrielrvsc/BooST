#' Predict Method for Smooth Tree Fits
#'
#' Obtains predictions for a Smooth Tree model and a given design matrix.
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

predict.SmoothTree=function(object,newx=NULL,...){

  if(is.null(newx)){
    return(stats::fitted(object))
  }

  if(is.vector(newx)){newx=matrix(newx,nrow=1)}

  model=object$tree
  fitted.values=eval_tree(newx,model)

  return(fitted.values)
}
