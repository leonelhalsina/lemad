cla_lemad_loglik_rhs <- function(t,y,parameter){
  ly <- length(y)
  d <- ly/2
  Es <- y[1:d]
  Ds <- y[(d+1):ly]
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  Q <- parameter[[3]]
  diag(Q) <- 0
  
  all_states <- cbind(Ds,Es)
  a <- cbind(all_states[,2],all_states[,1])
  b <- t(all_states)
  cross_D_E <- a%*%b
  
  dD <- -((unlist(lapply(lambdas,sum))) + mus + Q %*% (rep(1,d)))  * Ds +( Q %*% Ds ) + unlist(lapply(lapply(lambdas,"*",cross_D_E),sum))
  dE <- -((unlist(lapply(lambdas,sum))) + mus + Q %*% (rep(1,d)))  * Es + ( Q %*% Es ) + mus + unlist(lapply(lapply(lambdas,"*",Es%*%t(Es)),sum))
  
  return(list(c(dE,dD)))
}

cla_doParalThing <- function(take_ancesSub,
                             states,
                             loglik,
                             forTime,
                             parameter,
                             use_fortran,
                             methode,
                             phy
                             
){
  #cl <- makeCluster(2)
  #registerDoParallel(cl)
  
  
  
  #.packages=c("foreach"),
  ii <- NULL
  rm(ii)
  statesNEW <- foreach::foreach (ii = 1:2,
                                 .packages = c(
                                   "lemad",
                                   #"diversitree",
                                   "deSolve",
                                   "phylobase",
                                   "foreach",
                                   "doParallel",
                                   "geiger",
                                   "apTreeshape"),
                                 .export = c(
                                   "cla_lemad_loglik",
                                   "ode_FORTRAN",
                                   #"phy",
                                   #"methode",
                                   "cla_calThruNodes"
                                   #,"use_fortran")) %dopar% {
                                 )) %dopar% { 
                                   ancesSub <- take_ancesSub[[ii]]
                                   for(i in 1:length(ancesSub)){
                                     calcul <- 
                                       cla_calThruNodes(ancesSub[i],
                                                        states,
                                                        loglik,
                                                        forTime,
                                                        parameter,
                                                        use_fortran = use_fortran,
                                                        methode = methode,
                                                        phy = phy)
                                     loglik <- calcul$loglik
                                     states <- calcul$states
                                   }
                                   list(states, loglik)
                                 }
  return(statesNEW)
}


cla_lemad_loglik <- function(parameter,
                             phy,
                             traits,
                             num_max_multiregion,
                             cond = "proper_cond",
                             root_state_weight = "proper_weights",
                             sampling_fraction,
                             setting_calculation = NULL,
                             see_ancestral_states = FALSE,
                             loglik_penalty = 0,
                             is_complete_tree = FALSE,
                             num_threads = 1,
                             method = ifelse(num_threads == 1,
                                             "odeint::bulirsch_stoer",
                                             "odeint::runge_kutta_fehlberg78"),
                             atol = 1e-16,
                             rtol = 1e-16){
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  
  
  if (is.null(setting_calculation)) {
    check_input(traits,phy,num_max_multiregion,sampling_fraction,root_state_weight)
    setting_calculation <- build_initStates_time(phy,
                                                 traits,
                                                 num_max_multiregion,
                                                 sampling_fraction,
                                                 is_complete_tree,
                                                 mus)
  }
  states <- setting_calculation$states
  forTime <- setting_calculation$forTime  # nolint
  ances <- setting_calculation$ances
  
  
  ##    if (num_concealed_states != round(num_concealed_states)) {
  ##   # for testing
  ##   d <- ncol(states) / 2
  ##   new_states <- states[, c(1:sqrt(d), (d + 1):((d + 1) + sqrt(d) - 1))]
  ##   new_states <- states[, c(1, 2, 3, 10, 11, 12)]
  ##   states <- new_states
  ## }
  
  loglik <- 0
  d <- ncol(states) / 2
  
  calcul <- c()
  if (num_threads == 1) {
    ancescpp <- ances - 1
    forTimecpp <- forTime # nolint
    forTimecpp[, c(1, 2)] <- forTimecpp[, c(1, 2)] - 1 # nolint
    calcul <- cla_calThruNodes_cpp(ancescpp,
                                   states,
                                   forTimecpp,
                                   lambdas,
                                   mus,
                                   Q,
                                   method,
                                   atol,
                                   rtol,
                                   is_complete_tree)
  } else {
    # because C++ indexes from 0, we need to adjust the indexing:
    ancescpp <- ances - 1
    forTimecpp <- forTime # nolint
    forTimecpp[, c(1, 2)] <- forTimecpp[, c(1, 2)] - 1 # nolint
    
    if (num_threads == -2) {
      calcul <- calc_cla_ll_threaded(ancescpp,
                                     states,
                                     forTimecpp,
                                     lambdas,
                                     mus,
                                     Q,
                                     1,
                                     method,
                                     is_complete_tree)
    } else {
      calcul <- calc_cla_ll_threaded(ancescpp,
                                     states,
                                     forTimecpp,
                                     lambdas,
                                     mus,
                                     Q,
                                     num_threads,
                                     method,
                                     is_complete_tree)
    }
  }
  
  mergeBranch <- calcul$mergeBranch # nolint
  nodeM <- calcul$nodeM  # nolint
  loglik <- calcul$loglik
  
  ## At the root
  
  ## At the root
  mergeBranch2 <- mergeBranch # nolint
  lmb <- length(mergeBranch2)
  if (is.numeric(root_state_weight)) {
    weightStates <- rep(root_state_weight / num_concealed_states, # nolint
                        num_concealed_states)
  } else {
    if (root_state_weight == "maddison_weights") {
      weightStates <- mergeBranch / sum(mergeBranch2)
    }
    if (root_state_weight == "proper_weights") {
      numerator <- rep(NA, lmb)
      for (j in 1:lmb) {
        numerator[j] <- mergeBranch2[j] / sum(lambdas[[j]] *
                                                ((1 - nodeM[1:d]) %o% (1 - nodeM[1:d])))
      }
      weightStates <- numerator / sum(numerator) # nolint
    }
    if (root_state_weight == "equal_weights") {
      weightStates <- rep(1 / lmb, lmb) # nolint
    }
  }
  
  if (cond == "maddison_cond") {
    preCond <- rep(NA, lmb) # nolint
    for (j in 1:lmb) {
      preCond[j] <- sum(weightStates[j] *
                          lambdas[[j]] *
                          (1 - nodeM[1:d][j]) ^ 2)
    }
    mergeBranch2 <- mergeBranch2 / sum(preCond) # nolint
  }
  
  if (is_complete_tree) {
    timeInte <- max(abs(ape::branching.times(phy))) # nolint
    y <- rep(0, lmb)
    
    nodeM <- ct_condition_cla(y, # nolint
                              timeInte,
                              lambdas,
                              mus,
                              Q,
                              "odeint::bulirsch_stoer",
                              1e-16,
                              1e-12)
    nodeM <- c(nodeM, y) # nolint
  }
  
  if (cond == "proper_cond") {
    preCond <- rep(NA, lmb) # nolint
    for (j in 1:lmb) {
      preCond[j] <- sum(lambdas[[j]] * ((1 - nodeM[1:d]) %o% (1 - nodeM[1:d])))
    }
    mergeBranch2 <- mergeBranch2 / preCond # nolint
  }
  
  wholeLike_atRoot <- sum(mergeBranch2 * weightStates) # nolint
  LL <- log(wholeLike_atRoot) + # nolint
    loglik -
    penalty(pars = parameter,
            loglik_penalty = loglik_penalty) 
  if(see_ancestral_states == TRUE){
    num_tips <- ape::Ntip(phy)
    
    
    states2 <- calcul$states
    ancestral_states <- states2[(num_tips+1):(nrow(states2)-1),]
    ancestral_states <- ancestral_states[,-(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- sort(ances)
    return(list(ancestral_states=ancestral_states,LL=LL))
  } else {
    return(LL)
  }
}
