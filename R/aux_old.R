grow_tree=function(x, y, p, d_max, gamma, node_obs){
  bf=0
  variables = sample(ncol(x), round(p*ncol(x)))#sort(sample(1:ncol(x), max(round(p*ncol(x)),2)  ))
  gammai=gamma[sample(1:length(gamma),1)]
  N=length(y)
  fit=list()
  for(i in 1:length(variables)){
    xtest=x[,variables[i]]
    #xtest=stats::runif(20,min(xtest)-0.1*stats::sd(xtest),max(xtest)+0.1*stats::sd(xtest))
    xtest=sample(xtest,min(20,N))

    #gammascale=max(stats::IQR(x[,variables[i]]),0.5)
    gammascale=max(stats::sd(x[,variables[i]]),0.1)
    ssr=sapply(xtest,initial_node_var_test,x=x[,variables[i]],y=y,gamma=gammai/gammascale,node_obs=node_obs)
    ssr=t(ssr)
    best=which.min(ssr[,1])
    res0=c(xtest[best],gammai/gammascale,ssr[best,])
    names(res0)=c("c0","gamma","val","b0","b1")
    fit[[i]]=res0

  }
  best=which.min(lapply(fit,function(x)x["val"]))
  node=fit[[best]]

  nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
                      , "parent" = 0, "terminal" = "yes", variable = variables[best], id = 1,deep=1)
  noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
                       , "parent" = 0, "terminal" = "yes", variable = variables[best], id = 2,deep=1)
  tree=rbind(nodeleft,noderight)
  tree$terminal=as.character(tree$terminal)

  Pmat=1/(1+exp(-node["gamma"]*(x[,variables[best]]-node["c0"])))
  Pmat=cbind(Pmat,1-Pmat)

  ################
  iter=1
  while(iter<=d_max){

    gammai=gamma[sample(1:length(gamma),1)]
    terminal=which(tree$terminal=="yes")
    variables=sample(ncol(x), round(p*ncol(x)))#sort(sample(1:ncol(x), max(round(p*ncol(x)),2)))
    test=expand.grid(variables,terminal)
    colnames(test)=c("variable","terminal")

    t1=Sys.time()
    fit=list()
    for(i in 1:nrow(test)){
      xt=x[,test[i,"variable"]]
      #xtest=stats::runif(10,min(xt)-0.1*stats::sd(xt),max(xt)+0.1*stats::sd(xt))
      xtest=sample(xt,min(20,N),prob = Pmat[,test$terminal[i]]+0.01)

      gammascale=max(stats::sd(xt),0.1)
      middlenodes=which(is.na(colSums(Pmat)))
      if(length(xtest)<=1){
        fit[[i]]=c(val=Inf)
      }else{
        ssr=sapply(xtest,node_var_test,x=x[,test[i,"variable"]],y=y,gamma=gammai/gammascale,Pmat=Pmat,
                   terminal=test$terminal[i],middlenodes=middlenodes,deep=tree$deep[test$terminal[i]]+1,node_obs=node_obs)
        ssr=t(ssr)
        ssr[is.nan(ssr)]=Inf
        best=which.min(ssr[,1])
        res=c(xtest[best],gammai/gammascale,ssr[best,])
        names(res)=c("c0","gamma","val","b0","b1")
        fit[[i]]=res
      }


    }
    Sys.time()-t1
    best=which.min(lapply(fit,function(x)x["val"]))
    node=fit[[best]]
    if(bf==5){
      iter=d_max+1
      next
    }
    if(is.infinite(node["val"])){
      bf=bf+1
      next
    }

    nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
                        , "parent" = tree$id[test$terminal[best]], "terminal" = "yes", variable = test[best,"variable"],id=nrow(tree)+1,
                        deep=tree$deep[test$terminal[best]]+1)
    noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
                         , "parent" = tree$id[test$terminal[best]], "terminal" = "yes", variable = test[best,"variable"],id=nrow(tree)+2,
                         deep=tree$deep[test$terminal[best]]+1)


    tree$terminal[test$terminal[best]]="no"

    tree=rbind(tree,nodeleft,noderight)

    p0=(1/(1+exp(-node["gamma"]*(x[,nodeleft$variable]-node["c0"]))))*Pmat[,test$terminal[best]]
    p1=(1-(1/(1+exp(-node["gamma"]*(x[,noderight$variable]-node["c0"])))))*Pmat[,test$terminal[best]]
    Pmat=cbind(Pmat,p0,p1)
    Pmat[,test$terminal[best]]=NA
    tree$b[tree$terminal=="yes"] =  node[ c(6:length(node),4,5)]
    iter=iter+1
  }
  Pmat[is.na(Pmat)]=0
  fitted=Pmat%*%tree$b
  result=list(tree=tree,fitted.values=fitted)
  return(result)
}

initial_node_var_test=function(c0,x,y,gamma,node_obs){
  logit=1/(1+exp(-gamma*(x-c0)))
  b0=logit;b1=1-logit
  X=cbind(b0,b1)
  l0=length(which(b0>=0.5))
  l1=length(which(b1>=0.5))
  if(l0<node_obs | l1<node_obs){
    return(c(Inf,rep(NA,ncol(X))))
  }
  b=tryCatch(stats::coef(stats::.lm.fit(X,y)),error=function(e)Inf)
  if(is.infinite(b[1])){
    return(c(b[1],rep(NA,ncol(X))))
  }
  c(sum((y-X%*%b)^2),b)
}

node_var_test=function(c0,x,y,gamma,Pmat,terminal,middlenodes,deep,node_obs){
  logit=1/(1+exp(-gamma*(x-c0)))
  b0=logit*Pmat[,terminal]
  b1=(1-logit)*Pmat[,terminal]
  # b0=(1/(1+exp(-gamma*(x-c0))))*Pmat[,terminal]
  # b1=(1-b0)*Pmat[,terminal]
  X=cbind(b0,b1,Pmat[,-c(terminal,middlenodes)])
  l0=length(which(b0>=0.5^deep))
  l1=length(which(b1>=0.5^deep))
  if(l0<node_obs | l1<node_obs){
    return(c(Inf,rep(NA,ncol(X))))
  }
  b=tryCatch(stats::coef(stats::.lm.fit(X,y)),error=function(e)Inf)
  if(is.infinite(b[1])){
    return(c(b[1],rep(NA,ncol(X))))
  }
  e=sum((y-X%*%b)^2)
  c(e,b)
}

gradient_st=function(object,x,variable){
  tree=object$tree
  terminal=tree[which(tree$terminal=="yes"),]
  logimat=matrix(NA,nrow(x),nrow(terminal))

  for(i in 1:nrow(terminal)){
    node=terminal[i,]
    logit=1/(1+exp(-node$gamma*(x[,node$variable]-node$c0)))
    if(node$variable == variable){
      dlogit=node$gamma*logit*(1-logit)
    }else{
      dlogit=0
    }
    if(node$side==2){
      logit=1-logit
      dlogit=-dlogit
    }
    parent=node$parent
    while(parent!=0){
      node=tree[parent,]

      logitaux=1/(1+exp(-node$gamma*(x[,node$variable]-node$c0)))
      if(node$variable == variable){
        dlogitaux=node$gamma*logitaux*(1-logitaux)
      }else{
        dlogitaux=0
      }
      if(node$side==2){
        logitaux=1-logitaux
        dlogitaux=-dlogitaux
      }

      dlogit=logit*dlogitaux+dlogit*logitaux
      logit=logit*logitaux

      parent=node$parent
    }
    logimat[,i]=dlogit
  }

  fitted=logimat%*%terminal$b
  return(fitted)
}

eval_tree=function(x,tree){
  terminal=tree[which(tree$terminal=="yes"),]
  logimat=matrix(NA,nrow(x),nrow(terminal))
  for(i in 1:nrow(terminal)){
    node=terminal[i,]
    logit=1/(1+exp(-node$gamma*(x[,node$variable]-node$c0)))
    if(node$side==2){logit=1-logit}
    parent=node$parent
    while(parent!=0){
      node=tree[parent,]
      logitaux=1/(1+exp(-node$gamma*(x[,node$variable]-node$c0)))
      if(node$side==2){logitaux=1-logitaux}
      logit=logit*logitaux
      parent=node$parent
    }
    logimat[,i]=logit
  }

  fitted=logimat%*%terminal$b
  return(fitted)
}

