
transf_funcdefpar <- function(
  idparsfuncdefpar,
  functions_defining_params,
  idfactorsopt,
  trparsfix,
  trparsopt,
  idparsfix,
  idparsopt
){
  trparfuncdefpar <- NULL
  ids_all <- c(idparsfix, idparsopt)
  
  values_all <- c(trparsfix / (1 - trparsfix), trparsopt / (1 - trparsopt))
  a_new_envir <- new.env()
  x <- as.list(values_all) ## To declare all the ids as variables
  
  if (is.null(idfactorsopt)) {
    names(x) <- paste0("par_", ids_all)
  } else {
    names(x) <-
      c(paste0("par_", ids_all),
        paste0("factor_", idfactorsopt))
  }
  list2env(x, envir = a_new_envir)
  
  for (jj in 1:length(functions_defining_params)) {
    myfunc <- functions_defining_params[[jj]]
    environment(myfunc) <- a_new_envir
    #   local(myfunc,envir = a_new_envir)
    value_func_defining_parm <- local(myfunc(),envir = a_new_envir)

    ## Now, declare the variable that is just calculated, so it is
    ## available for the next calculation if needed
    y <- as.list(value_func_defining_parm)
    names(y) <- paste0("par_", idparsfuncdefpar[jj])
    list2env(y, envir = a_new_envir)
    
    if (is.numeric(value_func_defining_parm) == FALSE) {
      stop("something went wrong with the calculation of parameters in 'functions_param_struct'")
    }
    trparfuncdefpar <- c(trparfuncdefpar, value_func_defining_parm)
  }
  trparfuncdefpar <- trparfuncdefpar / (1 + trparfuncdefpar)
  rm(a_new_envir)
  return(trparfuncdefpar)
}


lemad_transform_parameters <-
  function(trparsopt,
           trparsfix,
           idparsopt,
           idparsfix,
           idparslist,
           structure_func
  ) {
    if(is.null(structure_func)==FALSE){
      idparsfuncdefpar <- structure_func[[1]] 
      functions_defining_params <- structure_func[[2]] 
    #  idfactorsopt <- structure_func[[3]] <- idfactorsopt
      
      if(length(structure_func[[3]])>1){
        idfactorsopt <- structure_func[[3]] 
      } else {
        if(structure_func[[3]] =="noFactor"){
          idfactorsopt <- NULL
        } else {
          idfactorsopt <- structure_func[[3]] 
        }
      }
      
      trparfuncdefpar <- transf_funcdefpar(idparsfuncdefpar = idparsfuncdefpar,
                                         functions_defining_params = functions_defining_params,
                                         idfactorsopt = idfactorsopt,
                                         trparsfix = trparsfix,
                                         trparsopt = trparsopt,
                                         idparsfix=idparsfix,
                                         idparsopt = idparsopt)
    }
    
    if(class(idparslist[[1]]) == "list"){ # when the ml function is called from cla_lemad
      trpars1 <- idparslist
      
      for(j in 1:nrow(trpars1[[3]])){
        trpars1[[1]][[j]][,] <- NA
      }
      
      for(j in 2:3){
        trpars1[[j]][] <- NA
      }
      
      
      if(length(idparsfix) != 0){
        
        for(i in 1:length(idparsfix)){
          
          for(j in 1:nrow(trpars1[[3]])){
            id <- which(idparslist[[1]][[j]] == idparsfix[i])
            trpars1[[1]][[j]][id] <- trparsfix[i]
          }
          for(j in 2:3) {
            id <- which(idparslist[[j]] == idparsfix[i])
            trpars1[[j]][id] <- trparsfix[i]
          }
        }
      }
      
      for(i in 1:length(idparsopt)){
        for(j in 1:nrow(trpars1[[3]])){
          id <- which(idparslist[[1]][[j]] == idparsopt[i])
          trpars1[[1]][[j]][id] <- trparsopt[i]
        }
        
        
        for(j in 2:3){
          id <- which(idparslist[[j]] == idparsopt[i])
          trpars1[[j]][id] <- trparsopt[i]
        }
      }
      
      ## structure_func part
      
      if(is.null(structure_func) == FALSE){
        for (i in 1:length(idparsfuncdefpar)) {
          for(j in 1:nrow(trpars1[[3]])){
            id <- which(idparslist[[1]][[j]] == idparsfuncdefpar[i])
            trpars1[[1]][[j]][id] <- trparfuncdefpar[i]
          }
          
          for (j in 2:3)
          {
            id <- which(idparslist[[j]] == idparsfuncdefpar[i])
            trpars1[[j]][id] <- trparfuncdefpar[i]
          }
        }
      }  
      
      pre_pars1<-list()
      pars1 <- list()
      
      for(j in 1:nrow(trpars1[[3]])){
        
        pre_pars1[[j]] <- trpars1[[1]][[j]][,]/(1 - trpars1[[1]][[j]][,])
      }
      
      pars1[[1]] <- pre_pars1
      for(j in 2:3){
        
        pars1[[j]] <- trpars1[[j]]/(1 - trpars1[[j]])
      }
      
    } else { #### when non-cla option is called
      
      trpars1 <- idparslist
      for (j in 1:3) {
        trpars1[[j]][] = NA
      }
      if (length(idparsfix) != 0) {
        for (i in 1:length(idparsfix)) {
          for (j in 1:3) {
            id <- which(idparslist[[j]] == idparsfix[i])
            trpars1[[j]][id] <- trparsfix[i]
            
          }
        }
      }
      for (i in 1:length(idparsopt)) {
        for (j in 1:3)
        {
          id <- which(idparslist[[j]] == idparsopt[i])
          trpars1[[j]][id] <- trparsopt[i]
          
        }
      }
      
      ## if structure_func part
      
      if(is.null(structure_func) == FALSE){
        for (i in 1:length(idparsfuncdefpar)) {
          
          for (j in 1:3)
          {
            id <- which(idparslist[[j]] == idparsfuncdefpar[i])
            trpars1[[j]][id] <- trparfuncdefpar[i]
            
          }
        }
      }  
      pars1 <- list()
      for (j in 1:3) {
        pars1[[j]] <- trpars1[[j]] / (1 - trpars1[[j]])
      }
    }
    return(pars1)
  }

lemad_loglik_choosepar <-
  function(trparsopt,
           trparsfix,
           idparsopt,
           idparsfix,
           idparslist,
           structure_func = structure_func,
           phy = phy,
           traits = traits,
           num_max_multiregion = num_max_multiregion,
           cond = cond,
                                    root_state_weight = root_state_weight, 
                                    sampling_fraction = sampling_fraction, 
                                    setting_calculation = setting_calculation,
                                    see_ancestral_states = see_ancestral_states, 
                                    loglik_penalty = loglik_penalty,
                                    is_complete_tree = is_complete_tree, 
                                    verbose = verbose,
                                    num_threads = num_threads,
                                    atol = atol,
                                    rtol = rtol,
                                    method = method) {
    alltrpars <- c(trparsopt, trparsfix)
    if (max(alltrpars) > 1 | min(alltrpars) < 0) {
      loglik = -Inf
    } else {
      pars1 <-
        lemad_transform_parameters(
          trparsopt,
          trparsfix,
          idparsopt,
          idparsfix,
          idparslist,
          structure_func
        )
      
      if(class(pars1[[1]]) == "list"){ # is the cla_ used?
        loglik <-
          cla_lemad_loglik(
            parameter = pars1,
            phy = phy,
            traits = traits,
            num_max_multiregion = num_max_multiregion,
                                        cond = cond,
                                        root_state_weight = root_state_weight,
                                        sampling_fraction = sampling_fraction,
                                        setting_calculation = setting_calculation, 
                                        see_ancestral_states = see_ancestral_states,
                                        loglik_penalty = loglik_penalty,
                                        is_complete_tree = is_complete_tree,
                                        num_threads = num_threads,
                                        atol = atol,
                                        rtol = rtol,
                                        method = method)
      } else {
        loglik <-
          lemad_loglik(
            parameter = pars1,
            phy = phy,
            traits = traits,
            use_fortran = use_fortran,
            methode = methode,
            cond = cond,
            root_state_weight = root_state_weight,
            sampling_fraction = sampling_fraction,
            run_parallel = run_parallel,
            setting_calculation = setting_calculation,
            setting_parallel = setting_parallel,
            see_ancestral_states = see_ancestral_states,
            loglik_penalty = loglik_penalty
          )
      }
      if (is.nan(loglik) || is.na(loglik)) {
        #print(trparsopt) ## new thing
        cat("There are parameter values used which cause numerical problems.\n")
        loglik <- -Inf
      }
    }
    if(verbose){
      cat(c(trparsopt / (1 - trparsopt),loglik),"\n")
    }
    return(loglik)
  }
# lemad_transform_parameters <- function(trparsopt,trparsfix,idparsopt,idparsfix,idparslist){
#   trpars1 <- idparslist
#   for(j in 1:3){
#     trpars1[[j]][] = NA
#   }
#   if(length(idparsfix) != 0){
#     
#     for(i in 1:length(idparsfix)){
#       for(j in 1:3) {
#         id <- which(idparslist[[j]] == idparsfix[i])
#         trpars1[[j]][id] <- trparsfix[i]
#         
#       }
#     }
#   }
#   for(i in 1:length(idparsopt))
#   {
#     for(j in 1:3)
#     {
#       id <- which(idparslist[[j]] == idparsopt[i])
#       trpars1[[j]][id] <- trparsopt[i]
#       
#     }
#   }
#   pars1 <- list()
#   for(j in 1:3)
#   {
#     pars1[[j]] <- trpars1[[j]]/(1 - trpars1[[j]])
#   }
#   return(pars1)
# }
# 
# lemad_loglik_choosepar <- function(trparsopt,trparsfix,idparsopt,idparsfix,idparslist,phy=phy,traits=traits,num_concealed_states=num_concealed_states,use_fortran=use_fortran,methode,cond=cond,root_state_weight=root_state_weight,sampling_fraction=sampling_fraction,setting_calculation=setting_calculation,run_parallel=run_parallel,setting_parallel=setting_parallel,see_ancestral_states=see_ancestral_states){
#   alltrpars <- c(trparsopt,trparsfix)
#   if(max(alltrpars) > 1 | min(alltrpars) < 0) {
#     loglik = -Inf
#   } else {
#     pars1 <- lemad_transform_parameters(trparsopt,trparsfix,idparsopt,idparsfix,idparslist)
#      
#       loglik <- lemad_loglik(parameter=pars1,phy=phy,traits=traits,num_concealed_states=num_concealed_states,use_fortran=use_fortran,methode = methode,cond=cond,root_state_weight=root_state_weight,sampling_fraction=sampling_fraction,run_parallel=run_parallel,setting_calculation=setting_calculation,setting_parallel=setting_parallel,see_ancestral_states=see_ancestral_states)
#           
#     if(is.nan(loglik) || is.na(loglik)){
#       print(trparsopt) ## new thing
#       cat("There are parameter values used which cause numerical problems.\n")
#       loglik <- -Inf
#     }
#   }
#   return(loglik)
# }



lemad_ml <- function(
  phy,
  traits,
  num_concealed_states,
  idparslist,
  idparsopt,
  initparsopt,
  idparsfix,
  parsfix,
  cond = "proper_cond",
  root_state_weight = "proper_weights",
  sampling_fraction,
  tol = c(1E-4, 1E-5, 1E-7),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  use_fortran = TRUE,
  methode = "ode45",
  optimmethod = 'simplex',
  num_cycles = 1,
  run_parallel = FALSE,
  loglik_penalty = 0
){
  structure_func<-NULL
  check_input(traits,phy,sampling_fraction,root_state_weight)
  
  if(is.matrix(traits)){
    cat("you are setting a model where some species had more than one trait state \n")
  }
  
  if(length(initparsopt)!=length(idparsopt)){
    stop("initparsopt must be the same length as idparsopt. Number of parameters to optimize does not match the number of initial values for the search")
  }
  
  if(length(idparsfix)!=length(parsfix)){
    stop("idparsfix and parsfix must be the same length.Number of fixed elements does not match the fixed figures")
  }
  
  if(anyDuplicated(c(idparsopt,idparsfix))!=0){
    stop("at least one element was asked to be both fixed and estimated ")
  }
  
  if(identical(as.numeric(sort(c(idparsopt,idparsfix))),as.numeric(sort(unique(unlist(idparslist)))))==FALSE){
    stop("All elements in idparslist must be included in either idparsopt or idparsfix ")
  }
  
  if(anyDuplicated(c(unique(sort(as.vector(idparslist[[3]]))),idparsfix[which(parsfix==0)]))!=0){
    cat("You set some transitions as impossible to happen","\n")
  }
  
  see_ancestral_states <- FALSE 
  
  #options(warn=-1)
  cat("Calculating the likelihood for the initial parameters.","\n")
  utils::flush.console()
  trparsopt <- initparsopt/(1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] = 1
  trparsfix <- parsfix/(1 + parsfix)
  trparsfix[which(parsfix == Inf)] = 1
  optimpars <- c(tol,maxiter)
  
  
  if(.Platform$OS.type=="windows" && run_parallel==TRUE){
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    setting_calculation <- build_initStates_time_bigtree(phy, traits, num_concealed_states, sampling_fraction)
    setting_parallel<-1
    on.exit(parallel::stopCluster(cl))
    }
  
  if(.Platform$OS.type=="unix" && run_parallel==TRUE){
    doMC::registerDoMC(2)
    setting_calculation <- build_initStates_time_bigtree(phy, traits, num_concealed_states, sampling_fraction)
    setting_parallel <- 1
  } 
  
  if(run_parallel == FALSE){
    setting_calculation <- build_initStates_time(phy,traits,num_concealed_states,sampling_fraction)
    setting_parallel <- NULL
  }
  if(optimmethod == 'subplex') {verbose <- TRUE} else {verbose <- FALSE}
  initloglik <- lemad_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparslist = idparslist, structure_func = structure_func, phy = phy, traits = traits,cond=cond,root_state_weight=root_state_weight,sampling_fraction=sampling_fraction,setting_calculation=setting_calculation,see_ancestral_states = see_ancestral_states, loglik_penalty = loglik_penalty,verbose = verbose)
  cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
  if(initloglik == -Inf)
  {
    stop("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.")
  } else {
    cat("Optimizing the likelihood - this may take a while.","\n")
    utils::flush.console()
    cat(setting_parallel,"\n")
    out <- DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = lemad_loglik_choosepar,trparsopt = trparsopt,idparsopt = idparsopt,trparsfix = trparsfix,idparsfix = idparsfix,idparslist = idparslist, structure_func = structure_func,phy = phy, traits = traits,num_concealed_states=num_concealed_states,use_fortran=use_fortran,methode = methode,cond=cond,root_state_weight=root_state_weight,sampling_fraction=sampling_fraction,setting_calculation=setting_calculation,run_parallel=run_parallel,setting_parallel=setting_parallel,see_ancestral_states=see_ancestral_states, num_cycles = num_cycles, loglik_penalty = loglik_penalty,verbose = verbose)
    if(out$conv != 0)
    {
      stop("Optimization has not converged. Try again with different initial values.\n")
    } else {
      MLpars1 <- lemad_transform_parameters(as.numeric(unlist(out$par)),trparsfix,idparsopt,idparsfix,idparslist, structure_func)
      out2 <- list(MLpars = MLpars1,ML = as.numeric(unlist(out$fvalues)),conv = out$conv)
    }
  }
  return(out2)
}
