#' Maximum likehood search under Lineage Extinction Model of Ancestral Distribution (LEMAD) 
#' @title Maximum likehood optimization for (lemad) and probabilities of ancestral lineages at regions
#' @param phylotree_recons phylogenetic tree of class phylo, ultrametric, rooted and with branch lengths.
#' @param species_presence a vector (of class character) with regions for each species. Coded in letters: "C"  "A"  "C"  "C" . Same order than tree tips
#' @param areas a vector that includes all the regions. Do not include any combination of locations
#' @param num_max_multiregion indicates the maximum number of regions a lineage can possible take at a given point in time. It can go from 2 to length(areas). Alternatively, it can be a vector of length 2: first element to be the maximum number of regions and the second a letter that represents "without any modern species" e.g.,num_max_multiregion = c(3,"D")  
#' @param condition_on_origin If NULL, the conditioning on the crown origin follows Herrera-Alsina et al. 2019 (Syst. Biol.). But the user can provide the hypothesized root state: condition_on_origin <- "A"
#' @param DEC_events DIVA and DEC models differ sligthly in the type of events that are allowed e.g., vacariance narrow 
#' @param missing_spp_areas a list informing the proportion of missing species (missing from the tree) per location. Only include incomplete locations in this list. See example. If tree is complete, use NULL.
#' @param lineage_extinction set "free" to have it estimated or provide a rate to fix it at.
#' @param initial_lambda Vector of length 2 for initial in-situ and vicariance rates to start the ML search. If NULL, the starting lambda will be taken from a Birth-death process.
#' @param initial_disperextirpation Intial values for dispersal-local extinction during the ML search. If NULL, it will be equal to lambda/5.
#' @return List with model's loglik, number of free parameters, estimates of rates, and ancestral locations probabilities.
#' @examples
#'# Example of how to set the arguments for a Maximum Likelihood search.
#'library(lemad)
#'library(DDD)
#'areas <- c("A", "B", "C")
#'set.seed(4)
#'phylotree_recons <- ape::rcoal(4, tip.label = 1:4)
#'species_presence <- c("B","C","A","AC")
#' DEC_events <- FALSE
#' # In case the tree has  80% of species inhabiting in A. And 95% of species in AC.
#' # Notice that regions in vector 1 and proportions in vector 2 have to be in the same order
#'missing_spp_areas <- list()
#'missing_spp_areas[[1]] <- c("A","AC")
#'missing_spp_areas[[2]] <- c(0.8,0.95)
#'
#'output <- lemad_analysis(
#'phylotree_recons,
#'species_presence,
#'areas,
#'num_max_multiregion = 3,
#'condition_on_origin = NULL,
#'DEC_events,
#'missing_spp_areas,
#'lineage_extinction =  0.005,
#'initial_lambda = NULL,
#'initial_disperextirpation = NULL
#')
#' 
#' output$model_ml #  -9.893469 the loglikelihood for the model
#' @export

lemad_analysis <- function(phylotree_recons,species_presence,areas,num_max_multiregion,
                           condition_on_origin = NULL,
                           DEC_events,
                           missing_spp_areas,
                           lineage_extinction = "free",
                           initial_lambda = NULL,
                           initial_disperextirpation = NULL){
  
  most_widespread_spp <- max(nchar(species_presence),na.rm = TRUE)
  
  unique_areas_tocheck <- sort(unique(unlist(strsplit(species_presence, split = ""))))
  
  if(class(num_max_multiregion[2]) == "character"){
    if(any(c(areas,unique_areas_tocheck) == num_max_multiregion[2]) == FALSE){ # is num_max_multiregion[2] in the areas or the presence vector?
      cat( "you are doing an analysis with -areas without any modern species-   \n")
    } else {
      stop ("if you are trying an analysis with areas with no modern species, double check your areas and presence vector: num_max_multiregion[2] should not be included in either vector ")
    }
    
  } 
  
  if(length(unique_areas_tocheck) != length(areas)){   
    stop("your vector of species presences has to be consistent with vector of areas")
  }
  
  if(num_max_multiregion[1] == 1){
    stop("num_max_multiregion should be higher than 1 so you can have multi-region species")
  }
  if(num_max_multiregion[1] < most_widespread_spp){
    
    stop("some of your species are in more areas than num_max_multiregion ")
  }
  
  if(class(lineage_extinction) != "character" && class(lineage_extinction) != "numeric"){
    stop(" lineage_extinction should be either free in quouting marks or a number")    
  }
  
  
  for(ik in 1:length(species_presence)){
    string_to_sort <- species_presence[ik]
    species_presence[ik] <- paste(sort(unlist(strsplit(string_to_sort, ""))), collapse = "")
    
  }
  
  id_rate_insitu <- 1
  id_rate_vicariance <- 2
  id_the_mu <- c(3,3)
  id_q_expansion <- 4
  id_q_contraction <- 4
if(class(num_max_multiregion[2]) == "character"){ # area with no modern species
    areas <- c(areas,num_max_multiregion[2]) 
}
 
  matrices <- prepare_full_lambdas_vicariance(areas,num_max_multiregion[1],id_rate_insitu,
                                              id_rate_vicariance,DEC_events)
  
  lambdas <- matrices$all_matrices
  all_area_combination <- matrices$all_area_combination
  matrices_names <- matrices$matrices_names
  
  
  # here, the mus will be sorted in such a way that species present in more than one location, cannot go extinct
  mus <- NULL
  for(ui in 1:length(lambdas)){
    if(ui <= length(areas)){
      mus <- c(mus,id_the_mu[1])
    } else {
      mus <- c(mus,id_the_mu[2])
    }
  }
  
  qs <- lemad_prepare_q_matrix(all_area_combination,matrices_names,
                               id_q_expansion,id_q_contraction)
  idparslist <- list() 
  idparslist[[1]] <- lambdas
  idparslist[[2]] <- mus
  idparslist[[3]] <- qs
  
  if(class(missing_spp_areas) == "list"){
    sampling_fraction <- rep(1,length(lambdas))
    
    if(length(missing_spp_areas ) != 2){
      stop ("missing_spp_areas should have to vectors: missing_spp_areas[[1]] and missing_spp_areas[[2]]; see ?lemad_analysis")
      
    }
    
    if(length(missing_spp_areas[[1]] ) != length(missing_spp_areas[[2]] )){
      stop("the two vectos in missing_spp_areas should be same length")
    }
    
    
    for(II in 1:length(missing_spp_areas[[1]])){
      
      replace_sampling <- which(missing_spp_areas[[1]][II] == matrices_names)
      sampling_fraction[replace_sampling] <- missing_spp_areas[[2]][II]
      
    }
  }
  if(is.null(missing_spp_areas)){
    sampling_fraction <- rep(1,length(lambdas))
  }  
  
  
  startingpoint <- bd_ML(brts = ape::branching.times(phylotree_recons))
  # intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  ######
  
  
  
  ####### this bit will use a set of initial values for ALL the parameters to be estimated
  ####### This is intended as a routine to continue an optimization that did not finish.
  #######  The idea is to use the last parameters shown at the R console
  if(length(initial_lambda) > 2 & class(initial_disperextirpation) == "character"){
    
    if(initial_disperextirpation != "continuing"){
      stop("check your initial_disperextirpation, it should be numerical")
    }
    
    if(class(lineage_extinction) == "character"){
      idparsopt <- c(1,2,3,4)
      initparsopt <- initial_lambda # initial_lambda is the vector of all the parameter estimates.
      idparsfix <- c(0)
      parsfix <- c(0)
      
    } else {
      idparsopt <- c(1,2,4)
      initparsopt <- c(initial_lambda[1],initial_lambda[2],initial_lambda[3]) # initial_lambda[3] is the dispersal rate when extinction is "free"
      idparsfix <- c(0,3)
      parsfix <- c(0,lineage_extinction)
      
    }
    
    # end of it. Now, the regular setup of the initial points for the optimization from scratch
    
  } else {
    
    if(is.null(initial_lambda )){
      intGuessLamba <- c(startingpoint$lambda0,startingpoint$lambda0)
      # intGuessMu <- startingpoint$mu0
      
    } else {
      if(class(initial_lambda) != "numeric"){
        stop ("intial_lambda should be a number")
      }
      if(length(initial_lambda) != 2){
        stop ("it needs an intitial value for each type of speciation insitu and vicariance e.g. c(0.01,0.01)")
      }
      
      
      intGuessLamba <- initial_lambda
    }
    
    if(is.null(initial_disperextirpation )){
      initial_disperextirpation_rate <- startingpoint$mu0/5
      # intGuessMu <- startingpoint$mu0
      
    } else {
      if(class(initial_disperextirpation) != "numeric"){
        stop ("intial_lambda should be a number")
      }
      
      initial_disperextirpation_rate <- initial_disperextirpation
    } 
    
    if(class(lineage_extinction) == "character"){
      idparsopt <- c(1,2,3,4)
      initparsopt <- c(intGuessLamba[1],intGuessLamba[2],intGuessMu,initial_disperextirpation_rate)
      idparsfix <- c(0)
      parsfix <- c(0)
      
    } else {
      idparsopt <- c(1,2,4)
      initparsopt <- c(intGuessLamba[1],intGuessLamba[2],initial_disperextirpation_rate)
      idparsfix <- c(0,3)
      parsfix <- c(0,lineage_extinction)
      
    }
    
  }
  
  tol <- c(1e-04, 1e-05, 1e-07)
  maxiter <- 1000 * round((1.25) ^ length(idparsopt))
  methode <- "ode45"
  optimmethod <- "subplex" # "simplex""   subplex"
  cond <- "proper_cond"
  
  if(is.null(condition_on_origin)){
    root_state_weight <- "proper_weights" 
  } else {
    if(any(condition_on_origin == matrices_names)){
      root_state_weight <- rep(0,length(matrices_names))
      root_state_weight[which(condition_on_origin == matrices_names)] <- 1
      
    } else {
      stop("your condition_on_origin should contain a valid region")
    }
    
  }
  is_complete_tree <- FALSE
  model <- cla_lemad_ml(
    phylotree_recons,
    species_presence,
    num_max_multiregion,
    idparslist,
    idparsopt,
    initparsopt,
    idparsfix,
    parsfix,
    cond,
    root_state_weight,
    sampling_fraction,
    tol = c(1e-04, 1e-05, 1e-07),
    maxiter = 1000 * round((1.25)^length(idparsopt)),
    optimmethod = "simplex",
    num_cycles = 1, 
    loglik_penalty = 0, 
    is_complete_tree = is_complete_tree, 
    verbose = (optimmethod == "subplex"),
    num_threads = 1,
    atol = 1e-12,
    rtol = 1e-12,
    method = "odeint::bulirsch_stoer")
  
  parameter <- model$MLpars
  
  
  estimates_anagenesis <- unique(unlist(model$MLpars[[1]]))[1]
  estimates_cladogenesis <- unique(unlist(model$MLpars[[1]]))[3]
  estimates_dispersalextirpation <- unique(as.vector(model$MLpars[[3]]))[3]
  estimates_or_fixes_extinction <- unique(model$MLpars[[2]])
  
  estimated_rates <- c(estimates_anagenesis,estimates_cladogenesis,
                       estimates_dispersalextirpation,
                       estimates_or_fixes_extinction)
  
  if(class(lineage_extinction) == "character"){
    names(estimated_rates) <- c("in-situ","vicariance","disper-extirp","extinction")
  } else {
    names(estimated_rates) <- c("in-situ","vicariance","disper-extirp","Fixed_extinction")
    
  }
  
  num_threads <- 1
  output <- cla_lemad_loglik(parameter, phylotree_recons, species_presence, num_max_multiregion,
                             cond = "proper_cond",
                             root_state_weight = "proper_weights",
                             sampling_fraction,
                             setting_calculation = NULL,
                             see_ancestral_states = TRUE,
                             loglik_penalty = 0,
                             is_complete_tree = FALSE,
                             num_threads = 1,
                             method = ifelse(num_threads == 1,
                                             "odeint::bulirsch_stoer",
                                             "odeint::runge_kutta_fehlberg78"),
                             atol = 1e-16,
                             rtol = 1e-16)
  
  
  ancestral_states <- output$ancestral_states
  colnames(ancestral_states) <- give_me_states_combination(areas,num_max_multiregion[1])
  
  return(list(model_ml = model$ML,
              number_free_pars = length(idparsopt),
              estimated_rates = estimated_rates,
              ancestral_states = ancestral_states,
              matrices = lambdas))
}

cla_lemad_ml <- function(
    phy,
    traits,
    num_max_multiregion,
    idparslist,
    idparsopt,
    initparsopt,
    idparsfix,
    parsfix,
    cond = "proper_cond",
    root_state_weight = "proper_weights",
    sampling_fraction,
    tol = c(1e-04, 1e-05, 1e-07),
    maxiter = 1000 * round((1.25)^length(idparsopt)),
    optimmethod = "simplex",
    num_cycles = 1, 
    loglik_penalty = 0, 
    is_complete_tree = FALSE, 
    verbose = (optimmethod == "subplex"),
    num_threads = 1,
    atol = 1e-12,
    rtol = 1e-12,
    method = "odeint::bulirsch_stoer") {
  cat("I am checking your input right now \n")
  check_input(traits,phy,num_max_multiregion,sampling_fraction,root_state_weight)
  structure_func <- NULL
  if(is.matrix(traits)){
    stop("matrix is not implemented yet \n")
  }
  
  if(length(initparsopt) != length(idparsopt)){
    stop("initparsopt must be the same length as idparsopt. Number of parameters to optimize does not match the number of initial values for the search")
  }
  
  if(length(idparsfix) != length(parsfix)){
    stop("idparsfix and parsfix must be the same length.Number of fixed elements does not match the fixed figures")
  }
  
  if(anyDuplicated(c(idparsopt,idparsfix)) != 0){
    stop("at least one element was asked to be both fixed and estimated ")
  }
  
  if(identical(as.numeric(sort(c(idparsopt,idparsfix))),as.numeric(sort(unique(unlist(idparslist)))))==FALSE){
    stop("All elements in idparslist must be included in either idparsopt or idparsfix ")
  }
  
  if(anyDuplicated(c(unique(sort(as.vector(idparslist[[3]]))),idparsfix[which(parsfix==0)]))!=0){
    cat("You set some transitions as impossible to happen","\n")
  }
  
  if(class(idparslist[[1]]) == "matrix"){ ## it is a tailor case otherwise
    idparslist[[1]] <- lemad_prepare_full_lambdas(traits,num_concealed_states,idparslist[[1]])
  }
  
  see_ancestral_states <- FALSE 
  
  #options(warn=-1)
  cat("Calculating the likelihood for the initial parameters.","\n")
  utils::flush.console()
  trparsopt <- initparsopt/(1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix/(1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  mus <- calc_mus(is_complete_tree,
                  idparslist,
                  idparsfix,
                  parsfix,
                  idparsopt,
                  initparsopt)
  optimpars <- c(tol, maxiter)
  
  setting_calculation <- build_initStates_time(phy, 
                                               traits,
                                               num_max_multiregion,
                                               sampling_fraction,
                                               is_complete_tree, 
                                               mus)
  setting_parallel <- NULL
  
  initloglik <- lemad_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparslist = idparslist, structure_func = structure_func,phy = phy, traits = traits,num_max_multiregion = num_max_multiregion ,cond=cond,root_state_weight = root_state_weight,sampling_fraction = sampling_fraction, 
                                       setting_calculation = setting_calculation,
                                       see_ancestral_states = see_ancestral_states, 
                                       loglik_penalty = loglik_penalty,
                                       is_complete_tree = is_complete_tree, 
                                       verbose = verbose,
                                       num_threads = num_threads,
                                       atol = atol,
                                       rtol = rtol,
                                       method = method)
  cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
  if(initloglik == -Inf)
  {
    stop("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.")
  } else {
    cat("Optimizing the likelihood - this may take a while.","\n")
    utils::flush.console()
    cat(setting_parallel,"\n")
    out <- DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = lemad_loglik_choosepar,trparsopt = trparsopt,idparsopt = idparsopt,trparsfix = trparsfix,idparsfix = idparsfix,idparslist = idparslist,structure_func = structure_func,phy = phy, traits = traits,num_max_multiregion = num_max_multiregion, cond = cond,
                          root_state_weight = root_state_weight,
                          sampling_fraction = sampling_fraction,
                          setting_calculation = setting_calculation,
                          see_ancestral_states = see_ancestral_states,
                          num_cycles = num_cycles,
                          loglik_penalty = loglik_penalty,
                          is_complete_tree = is_complete_tree, 
                          verbose = verbose,
                          num_threads = num_threads,
                          atol = atol,
                          rtol = rtol,
                          method = method)
    if(out$conv != 0)
    {
      stop("Optimization has not converged. Try again with different initial values.\n")
    } else {
      MLpars1 <- lemad_transform_parameters(as.numeric(unlist(out$par)),trparsfix,idparsopt,idparsfix,idparslist,structure_func)
      out2 <- list(MLpars = MLpars1,ML = as.numeric(unlist(out$fvalues)),conv = out$conv)
    }
  }
  return(out2)
}
