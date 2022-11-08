lemad_ml_func_def_pars <- function(phy,
                                    traits,
                                    num_concealed_states,
                                    idparslist,
                                    idparsopt,
                                    initparsopt,
                                    idfactorsopt,
                                    initfactors,
                                    idparsfix,
                                    parsfix,
                                    idparsfuncdefpar,
                                    functions_defining_params,
                                    cond = "proper_cond",
                                    root_state_weight = "proper_weights",
                                    sampling_fraction,
                                    tol = c(1E-4, 1E-5, 1E-7),
                                    maxiter = 1000 * round((1.25) ^ length(idparsopt)),
                                    use_fortran = TRUE,
                                    methode = "ode45",
                                    optimmethod = 'simplex',
                                    num_cycles = 1,
                                    run_parallel = FALSE,
                                    loglik_penalty = 0) {
  
  structure_func <- list()
  structure_func[[1]] <- idparsfuncdefpar
  structure_func[[2]] <- functions_defining_params
  if(is.null(idfactorsopt)){
    structure_func[[3]] <- "noFactor"
  } else {
    structure_func[[3]] <- idfactorsopt
  }
  
  see_ancestral_states <- FALSE
  if (is.null(idfactorsopt) == FALSE) {
    if (length(initfactors) != length(idfactorsopt)) {
      stop("idfactorsopt should have the same length as initfactors.")
    }
  }
  
  if (is.list(functions_defining_params) == FALSE) {
    stop(
      "The argument functions_defining_params should be a list of functions. See example and vignette"
    )
  }
  
  if (length(functions_defining_params) != length(idparsfuncdefpar)) {
    stop(
      "The argument functions_defining_params should have the same length than idparsfuncdefpar"
    )
  }
  
  if (is.matrix(traits)) {
    cat("You are setting a model where some species had more than one trait state \n")
  }
  
  if (length(initparsopt) != length(idparsopt)) {
    stop(
      "initparsopt must be the same length as idparsopt. Number of parameters to optimize does not match the number of initial values for the search"
    )
  }
  
  if (length(idparsfix) != length(parsfix)) {
    stop(
      "idparsfix and parsfix must be the same length.Number of fixed elements does not match the fixed figures"
    )
  }
  
  if (anyDuplicated(c(idparsopt, idparsfix, idparsfuncdefpar)) != 0) {
    stop("At least one element was asked to be fixed, estimated or a function at the same time")
  }
  
  if (identical(as.numeric(sort(
    c(idparsopt, idparsfix, idparsfuncdefpar)
  )), as.numeric(sort(unique(
    unlist(idparslist)
  )))) == FALSE) {
    stop(
      "All elements in idparslist must be included in either idparsopt or idparsfix or idparsfuncdefpar "
    )
  }
  
  if (anyDuplicated(c(unique(sort(
    as.vector(idparslist[[3]])
  )), idparsfix[which(parsfix == 0)])) != 0) {
    cat("You set some transitions as impossible to happen", "\n")
  }
  
  
  #options(warn = -1)
  cat("Calculating the likelihood for the initial parameters.", "\n")
  utils::flush.console()
  
  initparsopt2 <- c(initparsopt, initfactors)
  
  trparsopt <- initparsopt2 / (1 + initparsopt2)
  trparsopt[which(initparsopt2 == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  
  
  optimpars <- c(tol, maxiter)
  
  
  
  if (.Platform$OS.type == "windows" && run_parallel == TRUE) {
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    setting_calculation <-
      build_initStates_time_bigtree(phy, traits, num_concealed_states, sampling_fraction)
    setting_parallel <- 1
    on.exit(parallel::stopCluster(cl))
  }
  
  if (.Platform$OS.type == "unix" && run_parallel == TRUE) {
    doMC::registerDoMC(2)
    setting_calculation <-
      build_initStates_time_bigtree(phy, traits, num_concealed_states, sampling_fraction)
    setting_parallel <- 1
  }
  
  if (run_parallel == FALSE) {
    setting_calculation <-
      build_initStates_time(phy, traits, num_max_multiregion, sampling_fraction)
    setting_parallel <- NULL
  }
  
  if(optimmethod == 'subplex') {verbose <- TRUE} else {verbose <- FALSE}
  initloglik <-
    lemad_loglik_choosepar(
      trparsopt = trparsopt,
      trparsfix = trparsfix,
      idparsopt = idparsopt,
      idparsfix = idparsfix,
      idparslist = idparslist,
      structure_func = structure_func,
      phy = phy,
      traits = traits,
      num_concealed_states =
        num_concealed_states,
      use_fortran = use_fortran,
      methode = methode,
      cond = cond,
      root_state_weight = root_state_weight,
      sampling_fraction = sampling_fraction,
      setting_calculation =
        setting_calculation,
      run_parallel = run_parallel,
      setting_parallel = setting_parallel,
      see_ancestral_states = see_ancestral_states,
      loglik_penalty = loglik_penalty,
      verbose = verbose
    )
  cat("The loglikelihood for the initial parameter values is",
      initloglik,
      "\n")
  if (initloglik == -Inf)
  {
    stop(
      "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values."
    )
  } else {
    cat("Optimizing the likelihood - this may take a while.", "\n")
    utils::flush.console()
    cat(setting_parallel, "\n")
    out <-
      DDD::optimizer(
        optimmethod = optimmethod,
        optimpars = optimpars,
        fun = lemad_loglik_choosepar,
        trparsopt = trparsopt,
        idparsopt = idparsopt,
        trparsfix = trparsfix,
        idparsfix = idparsfix,
        idparslist = idparslist,
        structure_func = structure_func,
        phy = phy,
        traits = traits,
        num_concealed_states = num_concealed_states,
        use_fortran = use_fortran,
        methode = methode,
        cond = cond,
        root_state_weight = root_state_weight,
        sampling_fraction = sampling_fraction,
        setting_calculation = setting_calculation,
        run_parallel = run_parallel,
        setting_parallel = setting_parallel,
        see_ancestral_states = see_ancestral_states,
        num_cycles = num_cycles,
        loglik_penalty = loglik_penalty,
        verbose = verbose
      )
    if (out$conv != 0)
    {
      stop("Optimization has not converged. Try again with different initial values.\n")
    } else {
      MLpars1 <-
        lemad_transform_parameters(
          as.numeric(unlist(out$par)),
          trparsfix,
          idparsopt,
          idparsfix,
          idparslist,
          structure_func
        )
      out2 <-
        list(MLpars = MLpars1, ML = as.numeric(unlist(out$fvalues)),conv = out$conv)
    }
  }
  return(out2)
}
