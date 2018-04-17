# initial_node_var_test_op=function(c0,x,y,gamma){
#   logit=1/(1+exp(-gamma*(x-c0)))
#   b0=logit;b1=1-logit
#   X=cbind(b0,b1)
#   #b=solve(t(X)%*%X)%*%t(X)%*%y
#   b=tryCatch(coef(.lm.fit(X,y)),error=function(e)Inf)
#   if(is.infinite(b[1])){
#     return(c(b[1]))
#   }
#   c(sum((y-X%*%b)^2))
# }
#
# node_var_test_op=function(c0,x,y,gamma,Pmat,terminal,middlenodes){
#   b0=(1/(1+exp(-gamma*(x-c0))))*Pmat[,terminal]
#   b1=(1-b0)*Pmat[,terminal]
#   X=cbind(b0,b1,Pmat[,-c(terminal,middlenodes)])
#   b=tryCatch(coef(.lm.fit(X,y)),error=function(e)Inf)
#   if(is.infinite(b[1])){
#     return(c(b[1]))
#   }
#   e=sum((y-X%*%b)^2)
#   c(e)
# }
#

# grow_treeDE=function(x, y, p, d_max, gamma){
#
#   variables = sample(ncol(x), round(p*ncol(x)))#sort(sample(1:ncol(x), max(round(p*ncol(x)),2)  ))
#   gammai=gamma[sample(1:length(gamma),1)]
#
#   fit=list()
#   for(i in 1:length(variables)){
#     xtest=x[,variables[i]]
#     gammascale=max(IQR(x[,variables[i]]),0.5)
#     opt=DEoptim(initial_node_var_test_op,lower=min(xtest)-0.1*sd(xtest),upper=max(xtest)+0.1*sd(xtest),x=x[,variables[i]],y=y,gamma=gammai*gammascale,
#                 DEoptim.control(itermax = 100,reltol = 1e-4,steptol = 5,trace=FALSE))
#     res=initial_node_var_test(c0=opt$optim$bestmem,x=x[,variables[i]],y=y,gamma=gammai*gammascale)
#     res=c(opt$optim$bestmem,gammai*gammascale,res)
#     names(res)=c("c0","gamma","val","b0","b1")
#     fit[[i]]=res
#
#   }
#
#   best=which.min(lapply(fit,function(x)x["val"]))
#   node=fit[[best]]
#
#   nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
#                       , "parent" = 0, "terminal" = "yes", variable = variables[best], id = 1)
#   noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
#                        , "parent" = 0, "terminal" = "yes", variable = variables[best], id = 2)
#   tree=rbind(nodeleft,noderight)
#   tree$terminal=as.character(tree$terminal)
#
#   Pmat=1/(1+exp(-node["gamma"]*(x[,variables[best]]-node["c0"])))
#   Pmat=cbind(Pmat,1-Pmat)
#
#   ################
#   iter=1
#
#   while(iter<=d_max){
#
#     gammai=gamma[sample(1:length(gamma),1)]
#     terminal=which(tree$terminal=="yes")
#     variables=sample(ncol(x), round(p*ncol(x)))#sort(sample(1:ncol(x), max(round(p*ncol(x)),2)))
#     test=expand.grid(variables,terminal)
#     colnames(test)=c("variable","terminal")
#
#     fit=list()
#     for(i in 1:nrow(test)){
#       xtest=x[,test[i,"variable"]]
#
#       gammascale=max(IQR(xtest),0.5)
#       middlenodes=which(is.na(colSums(Pmat)))
#       if(length(xtest)==0){
#         fit[[i]]=c(val=Inf)
#       }else{
#
#         opt=DEoptim(node_var_test_op,lower=min(xtest)-0.1*sd(xtest),upper=max(xtest)+0.1*sd(xtest),
#                     x=x[,test[i,"variable"]],y=y,gamma=gammai*gammascale,Pmat=Pmat,
#                     terminal=test$terminal[i],middlenodes=middlenodes,
#                     DEoptim.control(itermax = 100,reltol = 1e-4,steptol = 5,trace=FALSE)
#         )
#         res=node_var_test(opt$optim$bestmem,x=x[,test[i,"variable"]],y=y,gamma=gammai*gammascale,Pmat=Pmat,
#                           terminal=test$terminal[i],middlenodes=middlenodes)
#
#         res=c(opt$optim$bestmem,gammai*gammascale,res)
#         names(res)=c("c0","gamma","val","b0","b1")
#         fit[[i]]=res
#       }
#
#     }
#
#     best=which.min(lapply(fit,function(x)x["val"]))
#     node=fit[[best]]
#     if(is.infinite(node["val"])){
#       next
#     }
#
#     nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
#                         , "parent" = tree$id[test$terminal[best]], "terminal" = "yes", variable = test[best,"variable"],id=nrow(tree)+1)
#     noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
#                          , "parent" = tree$id[test$terminal[best]], "terminal" = "yes", variable = test[best,"variable"],id=nrow(tree)+2)
#
#
#     tree$terminal[test$terminal[best]]="no"
#
#     tree=rbind(tree,nodeleft,noderight)
#
#     p0=(1/(1+exp(-node["gamma"]*(x[,nodeleft$variable]-node["c0"]))))*Pmat[,test$terminal[best]]
#     p1=(1-(1/(1+exp(-node["gamma"]*(x[,noderight$variable]-node["c0"])))))*Pmat[,test$terminal[best]]
#     Pmat=cbind(Pmat,p0,p1)
#     Pmat[,test$terminal[best]]=NA
#     tree$b[tree$terminal=="yes"] =  node[ c(6:length(node),4,5)]
#     iter=iter+1
#   }
#   Pmat[is.na(Pmat)]=0
#   fitted=Pmat%*%tree$b
#   result=list(tree=tree,fitted.values=fitted)
#   return(result)
# }
#
#
# grow_tree_full=function(x, y, p, d_max, gamma){
#
#   variables = sort(sample(1:ncol(x), max(round(p*ncol(x)),2)  ))
#   gammai=gamma[sample(1:length(gamma),1)]
#
#   fit=list()
#   for(i in 1:length(variables)){
#     xtest=unique(sort(x[,variables[i]]))
#     if(length(xtest)>10){
#       xtest=xtest[quantile(1:length(xtest),probs = seq(0,1,0.1))]
#     }
#
#     gammascale=max(IQR(x[,variables[i]]),0.1)
#     ssr=sapply(xtest,initial_node_var_test,x=x[,variables[i]],y=y,gamma=gammai*gammascale)
#     ssr=t(ssr)
#     best=which.min(ssr[,1])
#     res=c(xtest[best],gammai*gammascale,ssr[best,])
#     names(res)=c("c0","gamma","val","b0","b1")
#     fit[[i]]=res
#   }
#   best=which.min(lapply(fit,function(x)x["val"]))
#   node=fit[[best]]
#
#   nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
#                       , "parent" = 0, "terminal" = "yes", variable = variables[best], id = 1, layer=1)
#   noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
#                        , "parent" = 0, "terminal" = "yes", variable = variables[best], id = 2,layer=1)
#   tree=rbind(nodeleft,noderight)
#   tree$terminal=as.character(tree$terminal)
#
#   Pmat=1/(1+exp(-node["gamma"]*(x[,variables[best]]-node["c0"])))
#   Pmat=cbind(Pmat,1-Pmat)
#
#   ################
#   for(iter in 1:d_max){
#     if(max(tree$layer) < iter){
#       break
#     }
#     terminal=which(tree$terminal=="yes" & tree$layer==max(tree$layer))
#
#     for(t in 1:length(terminal)){
#       gammai=gamma[sample(1:length(gamma),1)]
#       variables=sort(sample(1:ncol(x), max(round(p*ncol(x)),2)))
#       test=data.frame(variables,terminal[t])
#       colnames(test)=c("variable","terminal")
#
#       t1=Sys.time()
#       fit=list()
#       for(i in 1:nrow(test)){
#         xt=x[,test[i,"variable"]]
#         nodeobs=which(apply(Pmat,1,which.max)==test$terminal[i])
#         xtest=unique(sort(xt[nodeobs]))
#
#         nodep=tree[tree$id==test$terminal[i],]
#         used=c()
#         if(nodep$variable==test[i,"variable"]) used=c(used,nodep$c0)
#
#         parent=ifelse(length(nodep$parent)!=0,nodep$parent,-1)
#         while(length(parent)!=0){
#           if(nodep$variable==test[i,"variable"]){
#             used=c(used,nodep$c0)
#           }
#           nodep=tree[tree$id==parent,]
#           parent=nodep$parent
#         }
#         xtest=setdiff(xtest,used)
#         if(length(xtest)>10){
#           xtest=xtest[quantile(1:length(xtest),probs = seq(0,1,0.1))]
#         }
#
#         gammascale=max(IQR(xt),0.1)
#         middlenodes=which(is.na(colSums(Pmat)))
#
#         if(length(xtest)<=1){
#           fit[[i]]=c(val=Inf)
#         }else{
#           ssr=sapply(xtest,node_var_test,x=x[,test[i,"variable"]],y=y,gamma=gammai*gammascale,Pmat=Pmat,
#                      terminal=test$terminal[i],middlenodes=middlenodes)
#           ssr=t(ssr)
#           best=which.min(ssr[,1])
#           res=c(xtest[best],gammai*gammascale,ssr[best,])
#           names(res)=c("c0","gamma","val","b0","b1")
#           fit[[i]]=res
#         }
#
#
#       }
#       Sys.time()-t1
#       best=which.min(lapply(fit,function(x)x["val"]))
#       node=fit[[best]]
#       if(is.infinite(node["val"])){
#         next
#       }
#
#       nodeleft=data.frame("side" = 1,"b" = node["b0"],"c0" = node["c0"],gamma = node["gamma"]
#                           , "parent" = tree$id[test$terminal[best]], "terminal" = "yes", variable = test[best,"variable"],id=nrow(tree)+1,layer=tree$layer[test$terminal[best]]+1)
#       noderight=data.frame(side = 2,  "b" = node["b1"],"c0" = node["c0"],gamma = node["gamma"]
#                            , "parent" = tree$id[test$terminal[best]], "terminal" = "yes", variable = test[best,"variable"],id=nrow(tree)+2,layer=tree$layer[test$terminal[best]]+1)
#
#
#       tree$terminal[test$terminal[best]]="no"
#
#       tree=rbind(tree,nodeleft,noderight)
#
#       p0=(1/(1+exp(-node["gamma"]*(x[,nodeleft$variable]-node["c0"]))))*Pmat[,test$terminal[best]]
#       p1=(1-(1/(1+exp(-node["gamma"]*(x[,noderight$variable]-node["c0"])))))*Pmat[,test$terminal[best]]
#       Pmat=cbind(Pmat,p0,p1)
#       Pmat[,test$terminal[best]]=NA
#       tree$b[tree$terminal=="yes"] =  node[ c(6:length(node),4,5)]
#
#     }
#
#   }
#   Pmat[is.na(Pmat)]=0
#   fitted=Pmat%*%tree$b
#   result=list(tree=tree,fitted.values=fitted)
#   return(result)
# }
