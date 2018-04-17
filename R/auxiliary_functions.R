
grow_tree=function(x, y, p, d_max, gamma){

  variables = sample(ncol(x), round(p*ncol(x)))#sort(sample(1:ncol(x), max(round(p*ncol(x)),2)  ))
  gammai=gamma[sample(1:length(gamma),1)]

  fit=list()
  for(i in 1:length(variables)){
    xtest=x[,variables[i]]
    xtest=stats::runif(10,min(xtest)-0.1*stats::sd(xtest),max(xtest)+0.1*stats::sd(xtest))

    gammascale=max(stats::IQR(x[,variables[i]]),0.5)
    ssr=sapply(xtest,initial_node_var_test,x=x[,variables[i]],y=y,gamma=gammai*gammascale)
    ssr=t(ssr)
    best=which.min(ssr[,1])
    res=c(xtest[best],gammai*gammascale,ssr[best,])
    names(res)=c("c0","gamma","val","b0","b1")
    fit[[i]]=res

  }
  best=which.min(lapply(fit,function(x)x["val"]))
  node=fit[[best]]

  nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
                      , "parent" = 0, "terminal" = "yes", variable = variables[best], id = 1)
  noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
                       , "parent" = 0, "terminal" = "yes", variable = variables[best], id = 2)
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
    i=1
    for(i in 1:nrow(test)){
      xt=x[,test[i,"variable"]]
      xtest=stats::runif(10,min(xt)-0.1*stats::sd(xt),max(xt)+0.1*stats::sd(xt))
      gammascale=max(stats::IQR(xt),0.1)
      middlenodes=which(is.na(colSums(Pmat)))
      if(length(xtest)<=1){
        fit[[i]]=c(val=Inf)
      }else{
        ssr=sapply(xtest,node_var_test,x=x[,test[i,"variable"]],y=y,gamma=gammai*gammascale,Pmat=Pmat,
                   terminal=test$terminal[i],middlenodes=middlenodes)
        ssr=t(ssr)
        best=which.min(ssr[,1])
        res=c(xtest[best],gammai*gammascale,ssr[best,])
        names(res)=c("c0","gamma","val","b0","b1")
        fit[[i]]=res
      }


    }
    Sys.time()-t1
    best=which.min(lapply(fit,function(x)x["val"]))
    node=fit[[best]]
    if(is.infinite(node["val"])){
      next
    }

    nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
                        , "parent" = tree$id[test$terminal[best]], "terminal" = "yes", variable = test[best,"variable"],id=nrow(tree)+1)
    noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
                         , "parent" = tree$id[test$terminal[best]], "terminal" = "yes", variable = test[best,"variable"],id=nrow(tree)+2)


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

initial_node_var_test=function(c0,x,y,gamma){
  logit=1/(1+exp(-gamma*(x-c0)))
  b0=logit;b1=1-logit
  X=cbind(b0,b1)
  #b=solve(t(X)%*%X)%*%t(X)%*%y
  b=tryCatch(stats::coef(stats::.lm.fit(X,y)),error=function(e)Inf)
  if(is.infinite(b[1])){
    return(c(b[1],rep(NA,ncol(X))))
  }
  c(sum((y-X%*%b)^2),b)
}

node_var_test=function(c0,x,y,gamma,Pmat,terminal,middlenodes){
  b0=(1/(1+exp(-gamma*(x-c0))))*Pmat[,terminal]
  b1=(1-b0)*Pmat[,terminal]
  X=cbind(b0,b1,Pmat[,-c(terminal,middlenodes)])
  b=tryCatch(stats::coef(stats::.lm.fit(X,y)),error=function(e)Inf)
  if(is.infinite(b[1])){
    return(c(b[1],rep(NA,ncol(X))))
  }
  e=sum((y-X%*%b)^2)
  c(e,b)
}

