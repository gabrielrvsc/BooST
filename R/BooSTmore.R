#' Adds more trees to an estimated BooST object
#'
#' Adds more trees to an estimated BooST object
#'
#' @param x Design matrix with explanatory variables.
#' @param y Response variable.
#' @param M Number of trees.
#' @param display If TRUE, displays iteration counter.
#'
#' @return An object with S3 class "Boost".
#' \item{Model}{A list with all trees.}
#' \item{fitted.values}{Final model fitted values.}
#' \item{brmse}{Boost rmse in each iteratiob.}
#' \item{Model}{A list with all trees.}
#' \item{ybar}{Average value of y used in the first iteration.}
#' \item{v}{Chosen learning rate.}
#' \item{rho}{Vector of gradient estimates for each iteration.}
#' \item{nvar}{Numver of variables in x}
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
# @seealso \code{\link{predict.BooST}}, \code{\link{smooth_tree}}, \code{\link{estimate_derivative}}


BooST.more = function(x, y, object, M, display = FALSE) {

  if(class(object)!="BooST"){
    stop("Object must be of class BooST")
  }

  params=as.list(object$call)

  d_max = ifelse(is.null(params$d_max),4,params$d_max)
  gamma = if(is.null(params$gamma)){seq(0.5,5,0.01)}else{eval(params$gamma)}
  Mold = ifelse(is.null(params$M),300,params$M)
  stochastic = ifelse(is.null(params$stochastic),FALSE,params$stochastic)
  s_prop = ifelse(is.null(params$s_prop),0.5,params$s_prop)
  node_obs = ifelse(is.null(params$node_obs),nrow(x)/200,params$node_obs)
  p = ifelse(is.null(params$p),2/3,params$p)
  v = ifelse(is.null(params$v),0.2,params$v)

  save_rho=object$rho
  ybar = object$ybar

  d_max=d_max-1
  N=length(y)
  phi=predict(object,x)

  brmse=object$brmse
  savetree=object$Model
  sq=seq(Mold+1,Mold+M,1)

  if(stochastic==TRUE){
    for(i in sq){
      s=sample(1:N,floor(N*s_prop),replace = FALSE)
      u=y-phi

      step=grow_tree(x=x[s,],y=u[s],p=p,d_max=d_max,gamma=gamma,node_obs=node_obs)
      fitstep=eval_tree(x,step$tree)
      rho=stats::coef(stats::lm(y[s]-phi[s]~-1+fitstep[s]))

      phitest=phi+v*rho*fitstep
      savetree[[i]]=step
      brmse[i]=sqrt(mean((y-phitest)^2))

      if(i>1){
        if(brmse[i]/brmse[i-1]>1.02){
          rho=0
          phitest=phi+v*rho*fitstep
          savetree[[i]]=step
          brmse[i]=sqrt(mean((y-phitest)^2))
          cat("stag")
        }
      }
      phi=phitest
      save_rho[i]=rho
      if(display==TRUE){
        cat(i," RMSE = ",brmse[i],"\n")
      }

    }

  }else{

    for(i in sq){
      u=y-phi
      step=grow_tree(x=x,y=u,p=p,d_max=d_max,gamma=gamma,node_obs=node_obs)
      fitstep=stats::fitted(step)
      rho=stats::coef(stats::lm(y-phi~-1+fitstep))
      phitest=phi+v*rho*fitstep
      savetree[[i]]=step
      brmse[i]=sqrt(mean((y-phitest)^2))

      # if(i>1){
      #   if(brmse[i]>brmse[i-1]){
      #     rho=0
      #     phitest=phi+v*rho*fitstep
      #     savetree[[i]]=step
      #     brmse[i]=sqrt(mean((y-phitest)^2))
      #     cat("stag")
      #   }
      # }
      phi=phitest
      save_rho[i]=rho
      if(display==TRUE){
        cat(i," RMSE = ",brmse[i],"\n")
      }
    }

  }

  result=list(Model=savetree,fitted.values=phi,brmse=brmse,ybar=ybar,v=v,rho=save_rho,nvar=ncol(x),varnames=colnames(x),call=match.call())
  class(result)="BooST"
  return(result)
}
