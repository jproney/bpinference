#Auto-generates a Stan file for inferring a model where rates have specified functional dependencies

#functional_deps is an array of strings encoding functions of  variable 'x1, x2, ..., x9' and constant parameters 'c1','c2',... 'c9'
#priors is a m-dimensional list of with the prior objects for all parameters associated with a given rate
#example prior object: p = alist(name="gamma",params=c(1,1))
#filename is name of generated Stan file
generate = function(functional_deps, priors, filename){
  
  template =  readChar("stan_template.txt", file.info("stan_template.txt")$size) #load the template file
  nFunc = length(functional_deps)
  funcs = rep(0,nFunc)
  for(i in 1:nFunc){
    # convert the expression
    exprn = functional_deps[[i]]
    stanstr = exp_to_stan(exprn)
    funcs[i] = sprintf("\t\t\tif(func_type[k] == %d){\n\t\t\t\tR_prime[k, P[k]] = %s; \n\t\t\t}\n",i,stanstr) 
  }
  
  for(dist in 1:length(priors)){
    
  }
  write(sprintf(template, paste(subs, collapse="")), filename)
}

# takes simple mathematical R expressions and turns them into a specific type of Stan code to fill in the template
exp_to_stan = function(exprn){
  if(is.atomic(exprn) || is.name(exprn)){
    first = exprn
  }
  else{
    first = exprn[[1]]
  }
  
  if(is.atomic(first)){
    return(first)
  }
  if(deparse(first) %in% c('+','-','/','*','^')){ #allowed operations
    return(paste("(",exp_to_stan(exprn[[2]]),")", deparse(first), "(", exp_to_stan(exprn[[3]]), ")", sep=" "))
  }
  if(deparse(first) == '('){
    return(exp_to_stan(exprn[[2]]))
  }
  if(!is.na(str_extract(deparse(first),"^c[1-9]$"))){ # function parameters
    num = strtoi(substr(str_extract(deparse(first),"^c[1-9]$"),2,2))-1
    if(num == 0){
      return(sprintf("R[event_idx[k]]", num))
    }
    return(sprintf("R[event_idx[k]+%d]", num))
  }
  if(!is.na(str_extract(deparse(first),"^x[1-9]$"))){ # variables
    num = strtoi(substr(str_extract(deparse(first),"^x[1-9]$"),2,2))
    return(sprintf("function_var[l, %d]", num))
  }
  else{
    stop("Invalid expression!")
  }
}
