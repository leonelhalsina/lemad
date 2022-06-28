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


cla_calThruNodes <- function(
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
  hmax <- NULL
  ly <- ncol(states)
  d <- ncol(states)/2
  
  #ances <- ances[1] #########################  REMOVE!!!!!!!!!!!!
  #desIndex <- 1
  
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
      nodeMN <- deSolve::ode(y = y, func = cla_lemad_loglik_rhs,
                             times = c(0,timeInte), parms = parameter,rtol = reltol, atol = abstol,
                             hmax = NULL,method = methode)
    } else {
      nodeMN <- ode_FORTRAN(y = y, func = "cla_secsse_runmod",
                            times = c(0,timeInte), parms = parameter, rtol = reltol, atol = abstol,
                            method = methode)
    }
    if(desIndex == 1){
      nodeN <- nodeMN
    }
    if(desIndex == 2){
      nodeM <- nodeMN 
    }
    # print(nodeMN)
  }
  ## At the node
  nodeM <- as.numeric(nodeM[2,-1])
  nodeN <- as.numeric(nodeN[2,-1])
  ff <- lemad_normalize_loglik(nodeM[(1:d) + d],loglik); nodeM[(1:d) + d] <- ff$probs; loglik <- ff$loglik
  ff <- lemad_normalize_loglik(nodeN[(1:d) + d],loglik); nodeN[(1:d) + d] <- ff$probs; loglik <- ff$loglik
  # cat(rbind(nodeM[(d+1):length(nodeM)],nodeN[(d+1):length(nodeN)]),"\n")
  
  all_states <- cbind(nodeM[(d + 1):length(nodeM)],nodeN[(d + 1):length(nodeN)])
  a <- cbind(all_states[,2],all_states[,1])
  b <- t(all_states)
  cross_M_N <- a%*%b
  
  # mergeBranch <- NULL
  #for(iii in 1:d){
  #combProb <- 0.5*sum(lapply(lambdas,"*",cross_M_N)[[1]]) ## multiplication of probabilities of both branches
  #mergeBranch <- c(mergeBranch,combProb)
  mergeBranch <- 0.5 * (unlist(lapply(lapply(lambdas,"*",cross_M_N),sum)))
  #}
  ff <- lemad_normalize_loglik(mergeBranch,loglik); mergeBranch <- ff$probs; loglik <- ff$loglik
  #sumD <- sum(mergeBranch)
  #mergeBranch <- mergeBranch/sumD
  #loglik <- loglik + log(sumD)
  #cat(mergeBranch,"\n")
  newstate <- nodeM[1:d] ## extinction probabilities
  newstate <- c(newstate,mergeBranch)
  states[focal,] <- newstate
  #print(parameter); print(loglik)
  return(list(states = states,loglik = loglik,mergeBranch = mergeBranch,nodeM = nodeM))
}

#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%

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
      #cat("I will check your input \n")
      check_input(traits,phy,num_max_multiregion,sampling_fraction,root_state_weight)
      setting_calculation <- 
        build_initStates_time_bigtree(phy, traits,num_max_multiregion, sampling_fraction)
    }
    
    states <- setting_calculation$states
    forTime <- setting_calculation$forTime
    ancesSub1 <- setting_calculation$ancesSub1
    ancesSub2 <- setting_calculation$ancesSub2
    ancesRest <- setting_calculation$ancesRest
    
    
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states) / 2
    take_ancesSub <- list(ancesSub1, ancesSub2)
    
    if(is.null(setting_parallel)){
      cl <- parallel::makeCluster(2)
      doParallel::registerDoParallel(cl)
    }
    
    statesNEW <- cla_doParalThing(take_ancesSub,
                                  states,
                                  loglik,
                                  forTime,
                                  parameter,
                                  use_fortran,
                                  methode,
                                  phy
    )
    
    comingfromSub1 <- statesNEW[[1]][[1]]
    comingfromSub2 <- statesNEW[[2]][[1]]
    loglik <- statesNEW[[1]][[2]] + statesNEW[[2]][[2]]
    
    thoseCalculated <- 
      which(comingfromSub2[, ncol(comingfromSub2)] > 0 &
              comingfromSub2[, ncol(comingfromSub2)] < 1 &
              (is.na(comingfromSub2[, ncol(comingfromSub2)]) == FALSE))
    
    comingfromSub1[thoseCalculated, ] <- comingfromSub2[thoseCalculated, ]
    states <- comingfromSub1
    
    for(i in 1:length(ancesRest)){
      calcul <- 
        cla_calThruNodes(ancesRest[i], states, loglik, forTime, parameter, use_fortran = use_fortran,methode = methode, phy = phy)
      states <- calcul$states
      loglik <- calcul$loglik
      
    }
  } else {
    if(is.null(setting_calculation)){
     # cat("I will check your input \n")
      check_input(traits,phy,num_max_multiregion,sampling_fraction,root_state_weight)
      setting_calculation <- build_initStates_time(phy,traits,num_max_multiregion,sampling_fraction)
    } 
    
    states <- setting_calculation$states
    forTime <- setting_calculation$forTime
    ances <- setting_calculation$ances
    
    
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states)/2
    
    for(i in 1:length(ances)){
      # cat("partialloglik",loglik,"\n")########
      calcul <- cla_calThruNodes(ances[i],states,loglik,forTime,parameter,use_fortran = use_fortran,methode = methode,phy = phy)
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
    weightStates <- root_state_weight/length(root_state_weight)

    
  } else {
    if(root_state_weight == "maddison_weights"){  
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }
    
    if(root_state_weight == "proper_weights"){
      numerator <- NULL
      for(j in 1:length(mergeBranch2)){
        numerator <- c(numerator,
                       (mergeBranch2[j]/(sum(lambdas[[j]] * (1 - nodeM[1:d][j])^2))))
      }
      denomin <- NULL
      for(j in 1:length(mergeBranch2)){
        denomin <- c(denomin,(mergeBranch2[j]/(sum(lambdas[[j]] * (1 -nodeM[1:d][j]) ^2))))
      }
      
      weightStates <- numerator/sum(denomin)
    }
    #     if(root_state_weight == "proper_weights"){
    #   weightStates <- NULL
    #   for(j in 1:length(mergeBranch2)){
    #     weightStates <- c(weightStates,
    #                      (mergeBranch2[j]/(sum(lambdas[[j]] * (1 - nodeM[1:d][j]) ^2)))/
    #                        sum((mergeBranch2[j]/(sum(lambdas[[j]] * (1 -nodeM[1:d][j]) ^2)))))
    #   }
    # }
    
    if(root_state_weight == "equal_weights"){  
      weightStates <- rep(1/length(mergeBranch2),length(mergeBranch2))
    }
  }  
  
  if(cond == "maddison_cond"){
    preCond <- NULL
    for(j in 1:length(weightStates)){
      preCond <- c(preCond,
                   sum(weightStates[j] * lambdas[[j]] *  (1 - nodeM[1:d][j]) ^ 2)
      )
    }
    mergeBranch2 <- 
      mergeBranch2/(sum(preCond))
  }
  
  if(cond == "proper_cond"){
    preCond <- NULL
    for(j in 1:length(mergeBranch2)){
      # preCond <- c(preCond,
      #            sum((lambdas[[j]] *  (1 - nodeM[1:d][j]) ^ 2))) 
      preCond <- c(preCond,
                   sum(lambdas[[j]] * ((1 - nodeM[1:d]) %o% (1 - nodeM[1:d])))) # new conditioning based on probabilities of different daughter combinations
    }
    mergeBranch2 <- 
      mergeBranch2/preCond
  }
  
  atRoot <- ((mergeBranch2) * (weightStates))
  
  wholeLike <- sum(atRoot)
  LL <- log(wholeLike) + loglik - penalty(pars = parameter,loglik_penalty = loglik_penalty)
  #print(unique(unlist(parameter[[1]]))); print(LL);  
  if(see_ancestral_states == TRUE){
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips+1):nrow(states),]
    ancestral_states <- ancestral_states[,-(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- sort(ances)
    return(list(ancestral_states=ancestral_states,LL=LL))
  } else {
    return(LL)
  }
}
