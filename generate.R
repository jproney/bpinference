#Auto-generates a Stan file for inferring a model where rates have specified functional dependencies
library(stringr)

#functional_deps is an array of strings encoding functions of  variable 'x1, x2, ..., x9' and constant parameters 'c1','c2',... 'c9'
#priors is list of with prior objects for all parameters in order of the rate they contribute to
#example prior object: p = list(name="gamma",params=c(1,1))
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
  dists = readLines("allowed_distributions.txt") #load distributions file
  nPri = length(priors)
  priorStrs = rep(0, nPri)
  for(p in 1:nPri){
    ind = which(!is.na(str_extract(dists, paste("^",priors[[p]]$name, " ",sep=""))))
    if(ind){
      nparm = strtoi(str_extract(dists[ind], "[1-9]+"))
      if(nparm != length(priors[[p]]$params)){
        stop("Incorrect number of parameters")
      }
      constraints = str_extract(dists[ind], "[_+]+")
      for(i in 1:nparm){
        if(priors[[p]]$params[i] <= 0 && substr(constraints,i,i) == "+"){
          stop("Negative parameter not allowed here!")
        }
      }
    }
    else{
      stop("Invalid prior name!")
    }
    priorStrs[p] = sprintf("\tR[%d] ~ %s(%s)\n", p, priors[[p]]$name, paste(priors[[p]]$params, collapse=","))
  }
  write(sprintf(template, paste(funcs, collapse=""),  paste(priorStrs,collapse="")), filename)
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
