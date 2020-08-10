
grow0 = function(x,y,p,N,gamma,node_obs){
  variables = sample(ncol(x), round(p*ncol(x)))#sort(sample(1:ncol(x), max(round(p*ncol(x)),2)  ))
  gammai=gamma[sample(1:length(gamma),1)]
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
  return(list(tree = tree, Pmat=Pmat))
}

grow = function(x,y,p,N,gamma,node_obs,tree,Pmat,d_max,iter,bf){
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

    fit[[i]]=c(val=Inf)
    ssr=sapply(xtest,node_var_test,x=x[,test[i,"variable"]],y=y,gamma=gammai/gammascale,Pmat=Pmat,
               terminal=test$terminal[i],middlenodes=middlenodes,deep=tree$deep[test$terminal[i]]+1,node_obs=node_obs)
    ssr=t(ssr)
    ssr[is.nan(ssr)]=Inf
    best=which.min(ssr[,1])
    res=c(xtest[best],gammai/gammascale,ssr[best,])
    names(res)=c("c0","gamma","val","b0","b1")
    fit[[i]]=res

  }
  Sys.time()-t1
  best=which.min(lapply(fit,function(x)x["val"]))
  node=fit[[best]]
  if(bf==5){
    iter=d_max+1
    return(list(tree = tree, Pmat = Pmat, iter = iter, bf = bf))
  }
  if(is.infinite(node["val"])){
    bf=bf+1
    return(list(tree = tree, Pmat = Pmat, iter = iter, bf = bf))
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
  return(list(tree = tree, Pmat = Pmat, iter = iter, bf = bf))
}

grow_tree=function(x, y, p, d_max, gamma, node_obs){
  bf=0

  N=length(y)

  t0 = grow0(x,y,p,N,gamma,node_obs)
  tree = t0$tree
  Pmat = t0$Pmat

  ################
  iter=1
  while(iter<=d_max){
    titer=grow(x,y,p,N,gamma,node_obs,tree,Pmat,d_max,iter,bf)
    tree = titer$tree
    Pmat = titer$Pmat
    iter = titer$iter
    bf = titer$bf
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

grow0_random = function(x,y,p,N,gamma,node_obs){
  gammai=gamma[sample(1:length(gamma),1)]

  variable = sample(ncol(x), 1)#sort(sample(1:ncol(x), max(round(p*ncol(x)),2)  ))


  xtest=x[,variable]

  good_tree = FALSE

  while(good_tree==FALSE){

    c0 = sample(xtest,1)

    gammascale=max(stats::sd(xtest),0.1)

    #ssr = initial_node_var_test(c0,xtest,y,gammai/gammascale,node_obs = node_obs)

    res=c(c0,gammai/gammascale,NA,NA,NA)
    names(res)=c("c0","gamma","val","b0","b1")

    node = res

    nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
                        , "parent" = 0, "terminal" = "yes", variable = variable, id = 1,deep=1)
    noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
                         , "parent" = 0, "terminal" = "yes", variable = variable, id = 2,deep=1)
    tree=rbind(nodeleft,noderight)
    tree$terminal=as.character(tree$terminal)

    Pmat=1/(1+exp(-node["gamma"]*(x[,variable]-node["c0"])))

    l0=length(which(Pmat>=0.5^nodeleft$deep))
    if(l0<node_obs){
      good_tree = FALSE
    }else{
      good_tree = TRUE
    }

  }

  Pmat=cbind(Pmat,1-Pmat)

  return(list(tree = tree, Pmat=Pmat))
}

grow_random = function(x,y,p,N,gamma,node_obs,tree,Pmat,d_max,iter){
  gammai=gamma[sample(1:length(gamma),1)]
  terminal=which(tree$terminal=="yes")

  variable = sample(ncol(x), 1)
  test=expand.grid(variable,terminal)
  colnames(test)=c("variable","terminal")
  xtest = x[,variable]

  good_tree = FALSE

  while(good_tree==FALSE){
    splitnode = sample(terminal,1)
    gammascale=max(stats::sd(xtest),0.1)
    middlenodes=which(is.na(colSums(Pmat)))

    c0 = sample(xtest,1)

    res=c(c0,gammai/gammascale,NA,NA,NA)
    names(res)=c("c0","gamma","val","b0","b1")


    node=res


    nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
                        , "parent" = splitnode, "terminal" = "yes", variable = variable,id=nrow(tree)+1,
                        deep=tree$deep[splitnode]+1)
    noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
                         , "parent" = splitnode, "terminal" = "yes", variable = variable,id=nrow(tree)+2,
                         deep=tree$deep[splitnode]+1)



    p0=(1/(1+exp(-node["gamma"]*(x[,nodeleft$variable]-node["c0"]))))*Pmat[,splitnode]
    p1=(1-(1/(1+exp(-node["gamma"]*(x[,noderight$variable]-node["c0"])))))*Pmat[,splitnode]

    l0=length(which(p0>=0.5^nodeleft$deep))
    l1=length(which(p1>=0.5^noderight$deep))
    if(l0<node_obs | l1<node_obs){
      good_tree = FALSE
    }else{
      good_tree = TRUE
    }

  }
  tree$terminal[splitnode]="no"

  tree=rbind(tree,nodeleft,noderight)

  Pmat=cbind(Pmat,p0,p1)
  Pmat[,splitnode]=NA
  iter=iter+1
  return(list(tree = tree, Pmat = Pmat, iter = iter))
}

grow_tree_random=function(x, y, p, d_max, gamma, node_obs){

  N=length(y)

  t0 = grow0_random(x,y,p,N,gamma,node_obs)
  tree = t0$tree
  Pmat = t0$Pmat

  ################
  iter=1
  while(iter<=d_max){
    titer=grow_random(x,y,p,N,gamma,node_obs,tree,Pmat,d_max,iter)
    tree = titer$tree
    Pmat = titer$Pmat
    iter = titer$iter
  }
  Pmat[is.na(Pmat)]=0

  ##

  betas = coef(lm.fit(Pmat,y, tol = 1e-3))
  betas[is.na(betas)] = 0
  tree$b = betas

  fitted=Pmat%*%tree$b
  result=list(tree=tree,fitted.values=fitted)
  return(result)
}




grow_tree_random_aux=function(x, y, p, d_max, gamma, node_obs){

  N=length(y)

  crit = 1
  while(length(crit)>0){
    t0 = grow0_random(x,y,p,N,gamma,node_obs)
    tree = t0$tree
    Pmat = t0$Pmat

    ################
    iter=1
    while(iter<=d_max){
      titer=grow_random(x,y,p,N,gamma,node_obs,tree,Pmat,d_max,iter)
      tree = titer$tree
      Pmat = titer$Pmat
      iter = titer$iter
    }
    Pmat[is.na(Pmat)]=0

    deep_crit = tree$deep
    for(kk in 1:length(deep_crit)){
      l=length(which(Pmat[,kk]>=0.5^deep_crit[kk]))
      if(l<node_obs){
        Pmat[,kk]=0
      }
    }


    aux = rowSums(Pmat)
    crit = which(aux==0)
    d_max = d_max-1
  }
  ##
  Pmat = Pmat/rowSums(Pmat)

  betas = coef(lm(y~-1 + Pmat))
  betas[is.na(betas)] = 0
  tree$b = betas

  fitted=Pmat%*%tree$b
  result=list(tree=tree,fitted.values=fitted)
  return(result)
}
