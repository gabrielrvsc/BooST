#' Creates analytical expression for partial effects
#'
#' Creates analytical expression for partial effects estimation of a BoosT or a SmoothTree object.
#'
#'
#' @param object A BooST or a SmoothTree object.
#' @param x Matrix of obsevations where the partial derivatives are desired.
#' @param variable Index of the variable the partial effects should be estimated.
#' @keywords BooST, Boosting, Smooth Tree, Partial Effects
#' @export
#' @examples
#' ## == to be made == ##
#'
# @seealso \code{\link{BooST}}, \code{\link{smooth_tree}}

estimate_derivatives=function(object,x,variable){
  if(is.vector(x)){x=matrix(x,nrow=1)}
  if(class(object)=="BooST"){
    v=object$v
    rho=object$rho
    vrho=v*rho
    partialeffects=lapply(object$Model,function(z) gradient_st(z,x,variable) )
    for(i in 1:length(partialeffects)){
      partialeffects[[i]]=partialeffects[[i]]*vrho[i]
    }
    partialeffects=Reduce('+',partialeffects)
  }
  if(class(object)=="SmoothTree"){
    partialeffects=gradient_st(object,x,variable)
  }
  colnames(partialeffects)=colnames(x)[variable]
  return(partialeffects)
}
