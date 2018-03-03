

st_expression <- function(node, side){

  parse(text = paste0('x', node$variable))
  gamma = node$gamma
  location = node$location

  if (side == 'right') {
    fuzzy = paste0('1/(1 + ','exp(-',gamma,'*(',parse(text = paste0('x', node$variable)),' - ',location,')))')
  } else {
    fuzzy = paste0('(1 - 1/(1 + ','exp(-',gamma,'*(',parse(text = paste0('x', node$variable)),' - ',location,'))))')
  }

  return(fuzzy)
}

tree_expression <- function(node, tree){

  if (is.na(node$right) == TRUE) {
    fit_right = paste0(st_expression(node = node, side = 'right'),'*',eval(node$beta1))
  } else {
    aux = tree_expression(node = tree[[node$right]], tree)
    fit_right = paste0(st_expression(node = node, side = 'right'),'*(',aux,')')
  }


  if (is.na(node$left) == TRUE) {
    fit_left = paste0(st_expression(node = node, side = 'left'),'*',eval(node$beta0))
  } else {
    aux = tree_expression(node = tree[[node$left]], tree)
    fit_left = paste0(st_expression(node = node, side = 'left'),'*(',aux,')')
  }

  final_fit = paste0('(',fit_right,'+',fit_left,')')
}

tree_derivative <- function(tree, variables) {

  tree_expr = tree_expression(tree[[1]], tree)

  idx = variables

  partial_d <- sapply(idx, function(i, expr) stats::D(expr, paste0('x', i)), expr = parse(text = tree_expr))

  return(list('Tree' = tree_expr, 'Gradient' = partial_d))
}

eval_gradient <- function(object, x, v=1, rho=1){

  x=as.numeric(x)
  fit = matrix(0, nrow = length(x), ncol = length(object))


  for (i in 1:length(object)) {
    for(j in 1:length(x)){
      fit[j,i] = eval({for (u in 1:length(x)) assign(paste0('x',u), x[u]); object[[i]]$Gradient[[j]]})
      #fit[2,i] = eval({x1 <- x[1]; x2 <- x[2]; object[[i]]$Gradient[[2]]})
    }
  }
  vrho=v*rho
  fit=t(t(fit)*vrho)
  return(rowSums(fit,na.rm = TRUE))
}
