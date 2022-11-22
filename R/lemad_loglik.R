check_input <- function(traits,phy,num_max_multiregion,sampling_fraction,root_state_weight){
  traitStates <- sort(unique(unlist(strsplit(traits, split = ""))))
  
  #traitStates <- give_me_states_combination(traitStates,num_max_multiregion)
  
  
  if(class(root_state_weight) == "numeric"){
    # if(length(root_state_weight) != length(sort(unique(traits)))){
    #   stop("you need to have as many elements in root_state_weight as traits")
    # }
    if(length(which(root_state_weight == 1)) != 1){
      stop("your root_state_weight needs only one 1")
    }
  } else {
    if(any(root_state_weight == "maddison_weights" | root_state_weight == "equal_weights"|
           root_state_weight == "proper_weights") == FALSE){
      stop("check root_state_weight for a typo or so")
    }
  }
  
  if(ape::is.rooted(phy) == FALSE){
    stop("the tree needs to be rooted")
  }
  
  if(ape::is.binary(phy) == FALSE){
    stop("the tree needs to be fully resolved")
  }
  if(ape::is.ultrametric(phy) == FALSE){
    stop("the tree needs to be ultrametric")
  }
  
  if(is.matrix(traits)){
    if(length(sampling_fraction) != length(sort(unique(traits[, 1])))){
      stop("sampling_fraction must have as many elements as traits you have")
    }
    
    if(all(sort(unique(as.vector(traits))) == sort(unique(traits[, 1]))) == 
       FALSE){
      stop(
        "Check your trait argument, if you have more than one column, make sure all your states are included in the first column"
      )
    }
  } 
  
  if(length(sort(unique(as.vector(traitStates)))) < 2)
  {
    stop("the trait has only one state")
  }
}


lemad_loglik_rhs <- function(t,y,parameter){
  ly <- length(y)
  d <- ly/2
  Es <- y[1:d]
  Ds <- y[(d+1):ly]
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  Q <- parameter[[3]]
  diag(Q) <- 0
  dE <- mus-(lambdas + mus + Q %*% (rep(1,d))) * Es + lambdas * Es * Es + ( Q %*% Es )
  dD <- -(lambdas + mus + Q %*% (rep(1,d)))  * Ds + 2 * lambdas * Es * Ds +( Q %*% Ds )
  return(list(c(dE,dD)))
}


build_initStates_time <- function(phy,
                                  traits,
                                  num_max_multiregion,
                                  sampling_fraction,
                                  is_complete_tree = FALSE,
                                  mus = NULL){ 
  if(is.matrix(traits)){
    stop("matrix option in traits not implemented yet")## trait might be a matrix when a species has more than one state
    traitStates <- sort(unique())
  }else{
    traitStates <- sort(unique(unlist(strsplit(traits, split = ""))))
  }
  
  
  
  #traits <- c("a","b","c","d")
  if(class(num_max_multiregion[2])=="character"){ # areas with no modern species
    traitStates <- c(traitStates,num_max_multiregion[2])
    traitStates <- give_me_states_combination(traitStates,num_max_multiregion[1])
  } else {
     traitStates <- give_me_states_combination(traitStates,num_max_multiregion) 
  }

  
  
  
  nb_tip <- ape::Ntip(phy)
  nb_node <- phy$Nnode
  states <- matrix(ncol = length(traitStates)*2,nrow=(nb_tip+nb_node))
  ly <- ncol(states)
  d <- ncol(states) / 2
  ## In a example of 3 states, the names of the colums would be like:
  ##
  ## colnames(states) <- c("E0A","E1A","E2A","E0B","E1B","E2B",
  ##                   "D0A","D1A","D2B","D0B","D1B","D2B")
  states[1:nb_tip,] <- 0
  if(is.matrix(traits)){
    stop("matrix option in traits not implemented yet")## I repeat the process of state assignation as many times as columns I have
    for(iv in 1:ncol(traits)){
      usetraits <- traits[,iv]
      if(anyNA(usetraits)){
        nas <- which( is.na(traits))
        for(iii in 1:length(nas) ){
          states[nas[iii],] <- c(1 - rep(sampling_fraction,num_concealed_states),
                                 rep(sampling_fraction,num_concealed_states))
        }
      }
      
      for(iii in 1:length(traitStates)){ # Initial state probabilities
        StatesPresents <- d + iii
        toPlaceOnes <- NULL
        for(jj in 1:(num_concealed_states-1)){
          toPlaceOnes <- c(toPlaceOnes, StatesPresents + (length(traitStates)* jj))
        }
        toPlaceOnes <- c(StatesPresents,toPlaceOnes)
        tipSampling <- 1 * sampling_fraction
        states[which(usetraits == traitStates[iii]),toPlaceOnes] <- tipSampling[iii]
      }
      
      for(iii in 1:nb_tip){
        states[iii,1:d] <- rep(1 - sampling_fraction,num_concealed_states)
      }
    }
  } else {
    if(anyNA(traits)){
      nas <- which( is.na(traits))
      for(iii in 1:length(nas) ){
        states[nas[iii],] <- c(1 - rep(sampling_fraction,num_concealed_states),
                               rep(sampling_fraction,num_concealed_states))
      }
    }
    ####
    for(iii in 1:length(traitStates)){ # Initial state probabilities
      
      tipSampling <- 1 * sampling_fraction
      states[which(traits == traitStates[iii]),length(traitStates) +iii] <- tipSampling[iii]
    }
    
    for(iii in 1:nb_tip){
      states[iii,1:d] <- rep(1 - sampling_fraction)
    }
  }
  #####
  phy$node.label <- NULL
  split_times <- sort(ape::branching.times(phy), decreasing = F)
  ances <- as.numeric(names(split_times))
  forTime <- matrix(NA,ncol = 3,nrow = nrow(phy$edge))
  forTime[,1:2] <- phy$edge
  
  for(ab in 1:length(ances)){
    focalTime <- ances[ab]
    desRows <- which(phy$edge[, 1] == focalTime)
    desNodes <- phy$edge[desRows, 2]
    rootward_age <- split_times[which(names(split_times) == focalTime)]
    for(desIndex in 1:2){
      ## to Find the time which the integration will be done in
      if(any(desNodes[desIndex] == names(split_times))){
        tipward_age <- split_times[which(names(split_times) == desNodes[desIndex])]
        timeInterv <- c(tipward_age, rootward_age)
      } else {
        timeInterv <- c(0, rootward_age)
      }  
      forTime[ which(forTime[,2] == desNodes[desIndex]),3] <- timeInterv[2]-timeInterv[1]
    }
  }
  return(list(
    states = states,
    ances = ances,
    forTime = forTime
  ))
}

calThruNodes <- function(
    ances,
    states,
    loglik,
    forTime,
    parameter,
    use_fortran,
    methode,
    phy
){
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  nb_node <- phy$Nnode
  reltol <- 1e-12
  abstol <- 1e-16
  ly <- ncol(states)
  d <- ncol(states) / 2
  focal <- ances
  desRows <- which(phy$edge[, 1] == focal)
  desNodes <- phy$edge[desRows, 2]
  
  nodeM <- numeric()
  nodeN <- numeric()
  
  for(desIndex in 1:2){
    y <- states[desNodes[desIndex],]
    timeInte <- forTime[which(forTime[,2] == desNodes[desIndex]),3]
    ##  To make the calculation in both lineages
    
    if(use_fortran == FALSE){
      nodeMN <- deSolve::ode(y = y, func = lemad_loglik_rhs,
                             times = c(0,timeInte), parms = parameter,  rtol = reltol, atol = abstol,
                             method = methode)
    } else {
      nodeMN <- ode_FORTRAN(y = y, func = "lemad_runmod",
                            times = c(0,timeInte), parms = parameter,  rtol = reltol, atol = abstol,
                            method = methode)
    }
    if(desIndex == 1){
      nodeN <- nodeMN
    }
    if(desIndex == 2){
      nodeM <- nodeMN 
    }
  }
  ## At the node
  nodeM <- as.numeric(nodeM[2,-1])
  nodeN <- as.numeric(nodeN[2,-1])
  ff <- lemad_normalize_loglik(nodeM[(1:d) + d],loglik); nodeM[(1:d) + d] <- ff$probs; loglik <- ff$loglik
  ff <- lemad_normalize_loglik(nodeN[(1:d) + d],loglik); nodeN[(1:d) + d] <- ff$probs; loglik <- ff$loglik
  #mergeBranch <- NULL
  #for(iii in 1:d){
  #  combProb <-  nodeM[iii + d] * nodeN[iii + d] * lambdas[iii] ## multiplication of probabilities of both branches
  #  mergeBranch <- c(mergeBranch,combProb)
  #}
  mergeBranch <- nodeM[(1:d) + d] * nodeN[(1:d) + d] * lambdas[(1:d)]
  
  ff <- lemad_normalize_loglik(mergeBranch,loglik); mergeBranch <- ff$probs; loglik <- ff$loglik
  #sumD <- sum(mergeBranch)
  #mergeBranch <- mergeBranch/sumD
  #loglik <- loglik + log(sumD)
  
  newstate <- nodeM[1:d] ## extinction probabilities
  newstate <- c(newstate,mergeBranch)
  states[focal,] <- newstate
  #print(parameter); print(loglik)
  return(list(states = states,loglik = loglik,mergeBranch = mergeBranch,nodeM = nodeM))
}
# calThruNodes <- function(ances,
#                          states,
#                          loglik,
#                          forTime,
#                          parameter,
#                          use_fortran,
#                          methode,
#                          phy){
#   
#   
#   lambdas <- parameter[[1]]
#   mus <- parameter[[2]]
#   
#   parameter[[3]][is.na(parameter[[3]])] <- 0
#   Q <- parameter[[3]]
#   
#   nb_node <- phy$Nnode
#   reltol <- 1e-10
#   abstol <- 1e-16
#   hmax <- NULL
#   ly <- ncol(states)
#   d <- ncol(states) / 2
#   
#   focal <- ances
#   desRows <- which(phy$edge[, 1] == focal)
#   desNodes <- phy$edge[desRows, 2]
#   
#   nodeM <- numeric()
#   nodeN <- numeric()
#   
#   for(desIndex in 1:2){
#     y <- states[desNodes[desIndex], ]
#     timeInte <- forTime[which(forTime[, 2] == desNodes[desIndex]), 3]
#     ##  To make the calculation in both lineages
#     
#     if(use_fortran == "Rsolver"){
#       nodeMN <- deSolve::ode(
#         y = y,
#         func = lemad_loglik_rhs,
#         times = c(0, timeInte),
#         parms = parameter,
#         rtol = reltol,
#         atol = abstol,
#         hmax = hmax,
#         method = methode
#       )
#       if(desIndex == 1){
#         nodeN <- nodeMN
#         
#       }
#       
#       if(desIndex == 2){
#         nodeM <- nodeMN
#       }
#       
#     }
#     
#     
#     if(use_fortran == "Fortran"){
#       nodeMN <- ode_FORTRAN(
#         y = y,
#         func = "lemad_loglik_rhs",
#         times = c(0, timeInte),
#         parms = parameter,
#         rtol = reltol,
#         atol = abstol,
#         hmax = hmax,
#         method = methode
#       )
#       if(desIndex == 1){
#         nodeN <- nodeMN
#         
#       }
#       
#       if(desIndex == 2){
#         nodeM <- nodeMN
#       }
#       
#     }
#     
#   }
#   ## At the node
#   
#   nodeM <- as.numeric(nodeM[2, -1])
#   nodeN <- as.numeric(nodeN[2, -1])
#   mergeBranch <- NULL
#   
#   for(iii in 1:d){
#     combProb <- 
#       nodeM[iii + d] * nodeN[iii + d] * lambdas[iii] ## multiplication of probabilities of both branches
#     mergeBranch <- c(mergeBranch, combProb)
#   }
#   
#   sumD <- sum(mergeBranch)
#   mergeBranch <- mergeBranch / sumD
#   
#   newstate <- nodeM[1:d] ## extinction probabilities
#   newstate <- c(newstate, mergeBranch)
#   
#   states[focal, ] <- newstate
#   loglik <- loglik + log(sumD)
#   return(list(
#     states = states,
#     loglik = loglik,
#     mergeBranch = mergeBranch,
#     nodeM = nodeM
#   ))
# }

doParalThing <- function(take_ancesSub,
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
  statesNEW <- foreach::foreach(ii = 1:2,
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
                                  "lemad_loglik",
                                  "ode_FORTRAN",
                                  #"phy",
                                  #"methode",
                                  "calThruNodes"
                                  #,"use_fortran")) %dopar% {
                                )) %dopar% { 
                                  ancesSub <- take_ancesSub[[ii]]
                                  for(i in 1:length(ancesSub)){
                                    calcul <- 
                                      calThruNodes(ancesSub[i], states, loglik, forTime, parameter, use_fortran = use_fortran,methode = methode, phy = phy)
                                    loglik <- calcul$loglik
                                    states <- calcul$states
                                  }
                                  list(states, loglik)
                                }
  return(statesNEW)
}

build_initStates_time_bigtree <- 
  function(phy, traits,num_max_multiregion,  sampling_fraction){
    if(is.matrix(traits)){
      stop("matrix option in traits not implemented yet")## trait might be a matrix when a species has more than one state
      ## trait might be a matrix when a species has more than one state
      traitStates <- sort(unique(traits[, 1]))
    } else{
      traitStates <- sort(unique(unlist(strsplit(traits, split = ""))))
    }
    
    
    traitStates <- give_me_states_combination(traitStates,num_max_multiregion)
    
    
    
    nb_tip <- ape::Ntip(phy)
    nb_node <- phy$Nnode
    states <- 
      matrix(ncol = (length(traitStates) * 2),
             nrow = (nb_tip + nb_node))
    
    
    ####
    ly <- ncol(states)
    d <- ncol(states) / 2
    ## In a example of 3 states, the names of the colums would be like:
    ##
    ## colnames(states) <- c("E0A","E1A","E2A","E0B","E1B","E2B",
    ##                   "D0A","D1A","D2B","D0B","D1B","D2B")
    
    states[1:nb_tip,] <- 0
    
    if(is.matrix(traits)){
      ## I repeat the process of state assignation as many times as columns I have
      for(iv in 1:ncol(traits)){
        usetraits <- traits[, iv]
        if(anyNA(usetraits)){
          nas <- which(is.na(traits))
          for(iii in 1:length(nas)){
            states[nas[iii],] <- c(1 - rep(sampling_fraction, num_concealed_states),
                                   rep(sampling_fraction, num_concealed_states))
          }
        }
        
        for(iii in 1:length(traitStates)){
          # Initial state probabilities
          StatesPresents <- d + iii
          toPlaceOnes <- NULL
          for(jj in 1:(num_concealed_states - 1)){
            toPlaceOnes <- 
              c(toPlaceOnes, StatesPresents + (length(traitStates) * jj))
          }
          toPlaceOnes <- c(StatesPresents, toPlaceOnes)
          tipSampling <- 1 * sampling_fraction
          states[which(usetraits == traitStates[iii]), toPlaceOnes] <- 
            tipSampling[iii]
        }
        
        for(iii in 1:nb_tip){
          states[iii, 1:d] <- rep(1 - sampling_fraction, num_concealed_states)
        }
        
      }
    } else {
      if(anyNA(traits)){
        nas <- which(is.na(traits))
        for(iii in 1:length(nas)){
          states[nas[iii],] <- c(1 - rep(sampling_fraction, num_concealed_states),
                                 rep(sampling_fraction, num_concealed_states))
        }
      }
      
      # for(iii in 1:length(traitStates)){
      #   # Initial state probabilities
      #   StatesPresents <- d + iii
      #   toPlaceOnes <- NULL
      #   for(jj in 1:(num_concealed_states - 1)){
      #     toPlaceOnes <- 
      #       c(toPlaceOnes, StatesPresents + (length(traitStates) * jj))
      #   }
      #   toPlaceOnes <- c(StatesPresents, toPlaceOnes)
      #   tipSampling <- 1 * sampling_fraction
      #   states[which(traits == traitStates[iii]), toPlaceOnes] <- 
      #     tipSampling[iii]
      # }
      # 
      # for(iii in 1:nb_tip){
      #   states[iii, 1:d] <- rep(1 - sampling_fraction, num_concealed_states)
      # }
      
      
      
      for(iii in 1:length(traitStates)){ # Initial state probabilities
        
        tipSampling <- 1 * sampling_fraction
        states[which(traits == traitStates[iii]),length(traitStates) +iii] <- tipSampling[iii]
      }
      
      for(iii in 1:nb_tip){
        states[iii,1:d] <- rep(1 - sampling_fraction)
      }
      
    }
    
    phy$node.label <- NULL
    split_times <- sort(ape::branching.times(phy), decreasing = F)
    ances <- as.numeric(names(split_times))
    
    forTime <- matrix(NA, ncol = 3, nrow = nrow(phy$edge))
    forTime[, 1:2] <- phy$edge
    
    for(ab in 1:length(ances)){
      focalTime <- ances[ab]
      desRows <- which(phy$edge[, 1] == focalTime)
      desNodes <- phy$edge[desRows, 2]
      
      rootward_age <- 
        split_times[which(names(split_times) == focalTime)]
      
      for(desIndex in 1:2){
        ## to Find the time which the integration will be done in
        if(any(desNodes[desIndex] == names(split_times))){
          tipward_age <- 
            split_times[which(names(split_times) == desNodes[desIndex])]
          timeInterv <- c(tipward_age, rootward_age)
        } else {
          timeInterv <- c(0, rootward_age)
        }
        
        forTime [which(forTime[, 2] == desNodes[desIndex]), 3] <- 
          timeInterv[2] - timeInterv[1]
      }
    }
    reltol <- 1e-12
    abstol <- 1e-16
    
    phySplit <- phy
    phySplit$node.label <- NULL
    nspp <- length(phySplit$tip.label)
    
    split_times <- 
      sort(ape::branching.times(phySplit), decreasing = T)
    interNodes <- as.numeric(names(split_times))
    
    formattedtree <- apTreeshape::as.treeshape(phySplit)
    smaller <- apTreeshape::smaller.clade.spectrum(formattedtree)
    smaller <- cbind(smaller, (nspp + 1):(nrow(phySplit$edge) + 1))
    
    optSplit <- NULL
    for(I in 1:length(interNodes)){
      optSplit <- 
        rbind(optSplit, smaller[which(smaller[, 3] == interNodes[I]),])
    }
    
    optSplit <- optSplit[-1,] # To prevent the root to be chosen
    NodeOptimSplit <- 
      optSplit[which(optSplit[, 2] == max(optSplit[, 2]))[1], 3]
    
    phy2 <- phylobase::phylo4(phy)
    
    sub1 <- as.numeric(phylobase::children(phy2,NodeOptimSplit)[1])
    sub2 <- as.numeric(phylobase::children(phy2,NodeOptimSplit)[2])
    jointSubs <- NodeOptimSplit
    
    cat("The best split is a subtree which has branches of size:",
        length(geiger::tips(phy, sub1)),"and",length(geiger::tips(phy, sub2)),"tips \n")    
    
    descenSub1 <- sort(phylobase::descendants(phy2, sub1, type = "all"))
    descenSub1 <- 
      as.numeric(descenSub1[(length(geiger::tips(phy, sub1)) + 1):length(descenSub1)])
    
    descenSub2 <- sort(phylobase::descendants(phy2, sub2, type = "all"))
    descenSub2 <- 
      as.numeric(descenSub2[(length(geiger::tips(phy, sub2)) + 1):length(descenSub2)])
    
    if(length(geiger::tips(phy, sub2)) == length(descenSub2)){ # the siblng node has no descendent nodes 
      descenSub2 <- sub2
    }
    
    ancesSub1 <- NULL
    for(i in 1:length(ances)){
      if(any(ances[i] == descenSub1))
      {
        ancesSub1 <- c(ancesSub1, ances[i])
      }
    }
    
    ancesSub2 <- NULL
    for(i in 1:length(ances)){
      if(any(ances[i] == descenSub2))
      {
        ancesSub2 <- c(ancesSub2, ances[i])
      }
    }
    
    ancesRest <- NULL
    for(i in 1:length(ances)){
      if(any(any(ances[i] == descenSub1) |
             any(ances[i] == descenSub2)) == FALSE)
      {
        ancesRest <- c(ancesRest, ances[i])
      }
    }
    
    return(
      list(
        states = states,
        ancesSub1 = ancesSub1,
        ancesSub2 = ancesSub2,
        ancesRest = ancesRest,
        forTime = forTime
      )
    )
  }

lemad_loglik <- function(parameter,
                         phy,
                         traits,
                         num_concealed_states,
                         use_fortran = TRUE,
                         methode = "ode45",
                         cond = "proper_cond",
                         root_state_weight = "proper_weights",
                         sampling_fraction,
                         run_parallel = FALSE,
                         setting_calculation = NULL,
                         setting_parallel= NULL,
                         see_ancestral_states = FALSE,
                         loglik_penalty = 0){
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  
  if(run_parallel == TRUE){ 
    if(is.null(setting_calculation)){
      check_input(traits,phy,sampling_fraction,root_state_weight)
      setting_calculation <- 
        build_initStates_time_bigtree(phy, traits, num_concealed_states, sampling_fraction)
    }
    
    states <- setting_calculation$states
    forTime <- setting_calculation$forTime
    ancesSub1 <- setting_calculation$ancesSub1
    ancesSub2 <- setting_calculation$ancesSub2
    ancesRest <- setting_calculation$ancesRest
    
    if(num_concealed_states != round(num_concealed_states)){ # for testing 
      d <- ncol(states) / 2 
      new_states <- states[,c(1:sqrt(d),(d + 1):((d + 1) + sqrt(d) - 1))]
      new_states <- states[,c(1,2,3,10,11,12)]
      states <- new_states
    }
    
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states) / 2
    take_ancesSub <- list(ancesSub1, ancesSub2)
    
    if(is.null(setting_parallel)){
      cl <- parallel::makeCluster(2)
      doParallel::registerDoParallel(cl)
    }
    
    statesNEW <- doParalThing(take_ancesSub,
                              states,
                              loglik,
                              forTime,
                              parameter,
                              use_fortran,
                              methode,
                              phy
    )
    
    comingfromSub1 <- statesNEW[[1]][[1]]
    comingfromSub2 <- statesNEW[[2]] [[1]]
    loglik <- statesNEW[[1]][[2]] + statesNEW[[2]][[2]]
    
    thoseCalculated <- 
      which(comingfromSub2[, ncol(comingfromSub2)] > 0 &
              comingfromSub2[, ncol(comingfromSub2)] < 1 &
              (is.na(comingfromSub2[, ncol(comingfromSub2)]) == FALSE))
    
    comingfromSub1[thoseCalculated, ] <- comingfromSub2[thoseCalculated, ]
    states <- comingfromSub1
    
    for(i in 1:length(ancesRest)){
      calcul <- 
        calThruNodes(ancesRest[i], states, loglik, forTime, parameter, use_fortran = use_fortran,methode=methode, phy= phy)
      states <- calcul$states
      loglik <- calcul$loglik
      
    }
  } else {
    
    if(is.null(setting_calculation)){
      check_input(traits,phy,sampling_fraction,root_state_weight)
      setting_calculation <- build_initStates_time(phy,traits,num_concealed_states,sampling_fraction)
    } 
    
    states <- setting_calculation$states
    forTime <- setting_calculation$forTime
    ances <- setting_calculation$ances
    
    if(num_concealed_states != round(num_concealed_states)){ # for testing 
      d <- ncol(states) / 2 
      new_states <- states[,c(1:sqrt(d),(d + 1):((d + 1) + sqrt(d) - 1))]
      new_states <- states[,c(1,2,3,10,11,12)]
      states <- new_states
    }
    
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states) / 2
    
    for(i in 1:length(ances)){
      calcul <- calThruNodes(ances[i],states,loglik,forTime,parameter,use_fortran = use_fortran,methode=methode, phy=phy)
      states <- calcul$states
      loglik <- calcul$loglik
      nodeN <- calcul$nodeN
    }
  }
  
  mergeBranch <- calcul$mergeBranch
  nodeM <- calcul$nodeM
  ## At the root
  mergeBranch2 <- (mergeBranch)
  if(class(root_state_weight) == "numeric"){
    giveWeights <- root_state_weight/num_concealed_states
    weightStates <- rep(giveWeights,num_concealed_states)
  } else {
    if(root_state_weight == "maddison_weights"){  
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }
    
    if(root_state_weight == "proper_weights"){
      weightStates <- (mergeBranch2/(lambdas * (1 - nodeM[1:d]) ^ 2))/sum((mergeBranch2/(lambdas * (1 - nodeM[1:d]) ^ 2)))
    }
    
    if(root_state_weight == "equal_weights"){  
      weightStates <- rep(1/length(mergeBranch2),length(mergeBranch2))
    }
  }  
  
  if(cond == "maddison_cond"){
    mergeBranch2 <- 
      mergeBranch2 / sum(weightStates * lambdas * (1 - nodeM[1:d]) ^ 2)
  }
  
  if(cond == "proper_cond"){
    mergeBranch2 <- mergeBranch2 / (lambdas * (1 - nodeM[1:d]) ^ 2)
  }
  
  atRoot <- ((mergeBranch2) * (weightStates))
  
  wholeLike <- sum(atRoot)
  LL <- log(wholeLike) + loglik - penalty(pars = parameter,loglik_penalty = loglik_penalty)
  
  if(see_ancestral_states == TRUE){
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):nrow(states),]
    ancestral_states <- ancestral_states[,-(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- sort(ances)
    return(list(ancestral_states = ancestral_states,LL = LL))
  } else {
    return(LL)
  }
}
