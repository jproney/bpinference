#Auto-generates a Stan file for inferring a model where rates have specified functional dependencies
library(stringr)

#functional_deps is an array of strings encoding functions of  variable 'x1, x2, ...' and constant parameters 'c1','c2',...
#priors is list of with prior objects for all parameters in order of the rate they contribute to
#example prior object: p = list(name="gamma",params=c(1,1), init = c(0,3), bounds=c(0,3))
#filename is name of generated Stan file
generate = function(functional_deps, priors, filename){
  
  template =  readChar("stan_template.txt", file.info("stan_template.txt")$size) #load the template file
  nFunc = length(functional_deps)
  funcs = rep(0,nFunc)
  
  dists = readLines("allowed_distributions.txt") #load distributions file
  nPri = length(priors)
  declStrs = rep(0,nPri)
  priorStrs = rep(0, nPri)
  
  for(i in 1:nFunc){
    # convert the expression
    exprn = functional_deps[[i]]
    stanstr = exp_to_stan(exprn, nPri)
    funcs[i] = sprintf("\t\tR[%d, P[%d]] = %s;\n",i,i,stanstr) 
  }

  for(p in 1:nPri){
    ind = which(!is.na(str_extract(dists, paste("^",priors[[p]]$name, " ",sep=""))))
    if(length(ind) > 0){
      nparm = strtoi(str_extract(dists[ind], "[1-9]+"))
      if(nparm != length(priors[[p]]$params)){
        stop("Incorrect number of parameters for prior distribution")
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
    priorStrs[p] = sprintf("\tTheta%d ~ %s(%s);\n", p, priors[[p]]$name, paste(priors[[p]]$params, collapse=","))
    if(!is.na(priors[[p]]$bounds)){
      declStrs[p] = sprintf("\treal<lower=%d, upper=%d> Theta%d;\n", priors[[p]]$bounds[1], priors[[p]]$bounds[2], p)
    }
  }
  write(sprintf(template, paste(declStrs, collapse=""), paste(funcs, collapse=""),  paste(priorStrs,collapse="")), filename)
}

# takes simple mathematical R expressions and turns them into a specific type of Stan code to fill in the template
exp_to_stan = function(exprn, maxParams){
  if(is.atomic(exprn) || is.name(exprn)){
    first = exprn
  }
  else{
    first = exprn[[1]]
  }
  
  if(is.atomic(first)){
    return(first)
  }
  if(deparse(first) %in% c('+','-','/','*','^','exp', 'log')){ #allowed operations
    if(length(exprn) > 2){
      return(paste("(",exp_to_stan(exprn[[2]], maxParams),")", deparse(first), "(", exp_to_stan(exprn[[3]], maxParams), ")", sep=" "))
    }
    else{
      return(paste(deparse(first),'(', exp_to_stan(exprn[[2]], maxParams), ')', sep=""))
    }
  }
  if(deparse(first) == '('){
    return(exp_to_stan(exprn[[2]], maxParams))
  }
  if(!is.na(str_extract(deparse(first),"^c[1-9]$"))){ # function parameters
    s = str_extract(deparse(first),"^c[1-9]+$")
    num = strtoi(substr(s,2,nchar(s)))
    if(num > maxParams){
      stop(paste("There is no prior for parameter",s,sep=" "))
    }
    return(sprintf("Theta%d", num))
  }
  if(!is.null(str_extract(deparse(first),"^x[1-9]$"))){ # variables
    s = str_extract(deparse(first),"^x[1-9]+$")
    num = strtoi(substr(s,2,nchar(s)))
    return(sprintf("function_var[i, %d]", num))
  }
  else{
    stop("Invalid expression!")
  }
}
