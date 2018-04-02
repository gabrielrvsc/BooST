#' @importFrom magrittr set_colnames %>%

new_node=function(tree, variable, left = NA, right = NA, parent = NA,
         gamma, location, beta0, beta1, side = NA, variance = NA){

  d = length(tree) + 1

  node = list('variable' = variable,
              'left' = left, 'right' = right,
              'parent' = parent, 'gamma' = gamma,
              'location' = location,
              'beta0' = beta0, 'beta1' = beta1,
              'side' = side, 'variance' = variance, "id"=d)



  tree[[d]] <- node

  if (is.na(parent) == FALSE) {
    if (side == 'right') {
      tree[[parent]]$right = d
      tree[[parent]]$beta1 = NA
    } else {
      tree[[parent]]$left = d
      tree[[parent]]$beta0 = NA
    }
  }
  return(tree)
}

grow_tree=function(x, y, p, d_max, gamma,tol,analitical_grad){

  if(analitical_grad==FALSE){
    grr=NULL
  }else{
    grr=objective_gradient
  }

  if(length(d_max)==0){
    d_max=sample(3:7,1)
  }

  tree = list()

  #### CRIAR PRIMEIRO NODE

  variables = sample(ncol(x), round(p*ncol(x)))
  first_gamma=sample(gamma,1)
  first_node = list()

  #par.start=c(0,0,0)
  for (i in 1:length(variables)) {
    par.start = c(sample(x[,variables[i]],1), sample(y,2))

    first_node[[i]] = stats::optim(par = par.start, fn = objective, gr=grr, method = 'BFGS',
                            x = x, y = y, tree = tree, variable = variables[i], parent = NA,
                            gamma = first_gamma*max(IQR(x[,variables[i]]),0.1), side = NA
                            ,control=list(reltol=tol))
  }
  best_fit = which.min(lapply(first_node, function(x) x$value))


  new_tree = new_node(tree = tree, variable = variables[best_fit], gamma = first_gamma,
                      location = first_node[[best_fit]]$par[1], beta0 = first_node[[best_fit]]$par[2],
                      beta1 = first_node[[best_fit]]$par[3])

  var_explained = 1 - stats::var(y - eval_tree(x, new_tree[[1]], new_tree))/stats::var(y)
  new_tree[[length(new_tree)]]$variance = var_explained

  #------------

  #### GROWTH and STOP criteria
  tree = new_tree
  for(i in 1:d_max){
    tree <- grow_node(x = x, y = y, tree = tree, p = p, gamma = gamma, tol=tol,analitical_grad)
  }

  return(tree)
}

objective=function(par, x, y, tree, variable, side, gamma, parent){

  new_tree = new_node(tree = tree, variable = variable, side = side,
                      gamma = gamma, location = par[1], beta0 = par[2],
                      beta1 = par[3], parent = parent)

  response = eval_tree(x, new_tree[[1]], new_tree)
  residual = (y - response)^2

  return(sum(residual))
}

objective_gradient=function(par, x, y, tree, variable, side, gamma, parent){

  grr = c()
  gradient = list()

  new_tree = new_node(tree = tree, variable = variable, side = side,
                      gamma = gamma, location = par[1], beta0 = par[2],
                      beta1 = par[3], parent = parent)

  response = eval_tree(x, new_tree[[1]], new_tree)

  chain_rule = eval_gradient_optim(par, x, new_tree[[1]], new_tree, gradient)

  grr[1] = sum(-2*(y - response)*chain_rule$location)
  grr[2] = sum(-2*(y - response)*chain_rule$beta0)
  grr[3] = sum(-2*(y - response)*chain_rule$beta1)

  return(grr)
}

grow_node=function(x, y, tree, gamma, p, tol,analitical_grad){

  if(analitical_grad==FALSE){
    grr=NULL
  }else{
    grr=objective_gradient
  }

  gamma=sample(gamma,1)

  pos = data.frame('son' = rep(0, length(tree) + 1), 'parent' = rep(0, length(tree) + 1))
  colnames(pos) = c('son', 'parent')

  j = 1
  for (i in 1:length(tree)) {
    if (is.na(tree[[i]]$left) == TRUE) {
      pos[j,2] = i
      pos[j,1] = 'left'
      j = j + 1
    }

    if (is.na(tree[[i]]$right) == TRUE) {
      pos[j,2] = i
      pos[j,1] = 'right'
      j = j + 1
    }
  }

  fit = list()
  variables = sample(ncol(x), round(p*ncol(x)))

  comb = pos %>% merge(variables, by=NULL) %>%
    set_colnames(c('side', 'parent', 'variable'))


  #par=c(0,0,0)
  for (i in 1:nrow(comb)) {
    par = c(sample(x[,comb[i,'variable']],1), sample(y,2))
    fit[[i]] = stats::optim(par = par, fn = objective, gr=grr, method = 'BFGS',
                     tree = tree, side = as.character(comb[i,'side']), gamma = gamma*max(IQR(x[, comb[i,'variable']]),0.1),
                     variable = comb[i,'variable'], parent = comb[i,'parent'], x = x, y = y
                     ,control=list(reltol=tol))
  }

  best_fit = which.min(lapply(fit, function(x) x$value))

  old_fit = eval_tree(x, tree[[1]], tree)

  new_tree = new_node(tree = tree, variable = comb[best_fit, 'variable'],
                      parent = comb[best_fit, 'parent'], gamma = gamma,
                      location = fit[[best_fit]]$par[1], beta0 = fit[[best_fit]]$par[2],
                      beta1 = fit[[best_fit]]$par[3], side = comb[best_fit, 'side'],
                      variance = NA)

  new_fit = eval_tree(x, new_tree[[1]], new_tree)

  var_explained = (sum((y - old_fit)^2) - sum((y - new_fit)^2))/sum((y - mean(y))^2)
  new_tree[[length(new_tree)]]$variance = var_explained

  return(new_tree)
}

logit_func=function(x, node, side){

  variable = node$variable
  x = x[,variable]

  gamma = node$gamma
  location = node$location

  fuzzy = 1/(1 + exp(-gamma*(x - location)))

  if (side == 'left') fuzzy = 1 - fuzzy

  return(fuzzy)
}

eval_tree=function(x, node, tree){

  if (is.na(node$right) == TRUE) {
    fit_right = logit_func(x = x, node = node, side = 'right')*node$beta1
  } else {
    aux = eval_tree(x, node = tree[[node$right]], tree)
    fit_right = logit_func(x = x, node = node, side = 'right')*aux
  }


  if (is.na(node$left) == TRUE) {
    fit_left = logit_func(x = x, node = node, side = 'left')*node$beta0
  } else {
    aux = eval_tree(x, node = tree[[node$left]], tree)
    fit_left = logit_func(x = x, node = node, side = 'left')*aux
  }

  final_fit = fit_right + fit_left

  return(final_fit)
}

eval_gradient_optim=function(par, x, node, tree, gradient){

  if (node$id == length(tree)) {

    gradient$beta1 = logit_func(x = x, node = node, side = 'right')
    gradient$beta0 = logit_func(x = x, node = node, side = 'left')
    func = exp(-node$gamma*(x-par[1]))
    gradient$location = node$gamma*func*(par[2] - par[3])/((1+func)^2)

  }

  if (is.na(node$right) != TRUE) {

    aux = eval_gradient_optim(par, x, node = tree[[node$right]], tree, gradient)
    gradient$beta1 = logit_func(x = x, node = node, side = 'right')*aux$beta1
    gradient$beta0 = logit_func(x = x, node = node, side = 'right')*aux$beta0
    gradient$location = logit_func(x = x, node = node, side = 'right')*aux$location

  }

  if (is.na(node$left) != TRUE) {

    aux = eval_gradient_optim(par, x, node = tree[[node$left]], tree, gradient)
    gradient$beta1 = logit_func(x = x, node = node, side = 'left')*aux$beta1
    gradient$beta0 = logit_func(x = x, node = node, side = 'left')*aux$beta0
    gradient$location = logit_func(x = x, node = node, side = 'left')*aux$location

  }

  return(gradient)
}

