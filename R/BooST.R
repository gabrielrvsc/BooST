#' Estimate BooST
#'
#' Estimates Boosting of Smooth Trees (BooST)
#'
#'
#' @param x Design matrix with explanatory variables.
#' @param y Response variable.
#' @param v Learning rate (default 0.2).
#' @param p Proportion of variables tested in each node split (default 2/3).
#' @param d_max Number of splits in each tree (default 4).
#' @param gamma Transiction function intensity. Bigger numbers makes the transition less smoth. The default is a sequence of values (0.5:5) to be randomized in each new node. Multiple values may be supplied in a vector to increase the model randomness.
#' @param M Number of trees.
#' @param display If TRUE, displays iteration counter.
#' @param stochastic If TRUE the model will be estimated using Stochasting Gradient Boosting.
#' @param s_prop Used only if stochastic=TRUE. Determines the proportion of data used in each tree.
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


BooST = function(x, y, v=0.2, p = 2/3, d_max = 4, gamma = seq(0.5,5,0.01),
                 M = 300, display=FALSE,
                 stochastic=FALSE,s_prop=0.5) {

  d_max=d_max-1
  N=length(y)
  phi=rep(mean(y),length(y))

  brmse=rep(NA,M)
  savetree=list()
  save_rho=rep(NA,M)


  if(stochastic==TRUE){
    for(i in 1:M){
      s=sample(1:N,floor(N*s_prop),replace = FALSE)
      u=y-phi

      step=grow_tree(x=x[s,],y=u[s],p=p,d_max=d_max,gamma=gamma)
      fitstep=eval_tree(x,step$tree)
      rho=stats::coef(stats::lm(y[s]-phi[s]~-1+fitstep[s]))

      phitest=phi+v*rho*fitstep
      savetree[[i]]=step
      brmse[i]=sqrt(mean((y-phitest)^2))

      if(i>1){
        if(brmse[i]/brmse[i-1]>1.05){
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
        cat(i,"\n")
      }

    }

  }else{

    for(i in 1:M){
      u=y-phi
      step=grow_tree(x=x,y=u,p=p,d_max=d_max,gamma=gamma)
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
        cat(i,"\n")
      }
    }

  }

  result=list(Model=savetree,fitted.values=phi,brmse=brmse,ybar=mean(y),v=v,rho=save_rho,nvar=ncol(x),varnames=colnames(x),call=match.call())
  class(result)="BooST"
  return(result)
}
