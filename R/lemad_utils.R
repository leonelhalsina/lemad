
lemad_id_paramPos <- function(traits,num_concealed_states){
  idparslist <- list()
  if(is.matrix(traits)){
    traits <- traits[,1]
  }
  
  ly <- length(sort(unique(traits))) * 2 * num_concealed_states
  d <- ly/2
  idparslist[[1]] <- 1:d
  idparslist[[2]] <- (d + 1):ly
  toMatrix <- 1
  matPos <- (ly + 1): (((d^2)-d) +  d * 2)
  for (i in 1:d){
    toMatrix <- c(toMatrix,matPos[(i*d-(d-1)):((i*d-(d-1)) + d)])
    
  }
  toMatrix <- toMatrix[1:d^2]
  Q <- matrix(toMatrix,ncol = d,nrow = d,byrow = T)
  diag(Q) <- NA
  idparslist[[3]] <- Q
  
  lab_states <- rep(as.character(sort(unique(traits))),num_concealed_states)
  
  lab_conceal <- NULL
  for (i in 1:num_concealed_states){
    
    lab_conceal <- c(lab_conceal,rep(LETTERS[i],length(sort(unique(traits)))))
  }
  
  
  statesCombiNames <- character()
  for ( i in 1:length(lab_states)){
    statesCombiNames <- c(statesCombiNames,paste0(lab_states[i],lab_conceal[i]))
    
  }
  colnames(idparslist[[3]]) <- statesCombiNames
  rownames(idparslist[[3]]) <- statesCombiNames
  names(idparslist) <- c("lambdas","mus","Q")
  names(idparslist[[1]]) <- statesCombiNames
  names(idparslist[[2]]) <- statesCombiNames
  return(idparslist)
}



lemad_q_doubletrans <- function(traits,masterBlock,diff.conceal){
  
  if(diff.conceal == TRUE && all(floor(masterBlock) == masterBlock,na.rm =  TRUE) == FALSE){
    integersmasterBlock <- floor(masterBlock)
    factorBlock <- signif(masterBlock-integersmasterBlock,digits =  2)
    
    factorstoExpand <- unique(sort(c(factorBlock)))
    factorstoExpand <- factorstoExpand[factorstoExpand>0]
    newshareFac <- (max(factorstoExpand*10) + 1):(max(factorstoExpand*10) + length(factorstoExpand))
    newshareFac <- newshareFac/10
    
    for(iii in 1:length(newshareFac)){  
      factorBlock[which(factorBlock == factorstoExpand[iii])] <- newshareFac[iii]  
      
    }
    
    ntraits <- length(sort(unique(traits)))
    uniqParQ <- sort(unique(c(floor(masterBlock))))
    uniqParQ2 <- uniqParQ[which(uniqParQ>0)]
    concealnewQ <- (max(uniqParQ2) + 1):(max(uniqParQ2) + length(uniqParQ2))
    
    for(iii in 1:length(concealnewQ)){  
      integersmasterBlock[which(integersmasterBlock == uniqParQ2[iii])] <- concealnewQ[iii]  
      
    }
    concealnewQMatr <- integersmasterBlock + factorBlock
    
    Q <- NULL
    for(i in 1:ntraits){
      Qrow <- NULL
      for(ii in 1:ntraits){
        entry <- masterBlock[i,ii]
        if(is.na(entry)){
          Qrow <- cbind(Qrow,masterBlock)
        } else{
          entry <- concealnewQMatr[i,ii]
          
          outDiagBlock <- matrix(0,ncol = ntraits,nrow = ntraits,byrow = T)
          diag(outDiagBlock) <- entry
          Qrow <- cbind(Qrow,outDiagBlock)
        }
        
      }
      Q <- rbind(Q,Qrow)
    }
  } else {
    
    ntraits <- length(sort(unique(traits)))
    uniqParQ <- sort(unique(c(masterBlock)))
    uniqParQ2 <- uniqParQ[which(uniqParQ>0)]
    concealnewQ <- (max(uniqParQ2) + 1):(max(uniqParQ2) + length(uniqParQ2))
    concealnewQMatr <- masterBlock
    for( I in 1: length(uniqParQ2)){
      uniqParQ2
      concealnewQMatr[concealnewQMatr == uniqParQ2[I]] <- concealnewQ[I]
    }
    
    Q <- NULL
    for(i in 1:ntraits){
      Qrow <- NULL
      for(ii in 1:ntraits){
        entry <- masterBlock[i,ii]
        if(is.na(entry)){
          
          Qrow <- cbind(Qrow,masterBlock)
        } else{
          
          if(diff.conceal == TRUE){
            entry <- concealnewQMatr[i,ii]
          }
          outDiagBlock <- matrix(0,ncol = ntraits,nrow = ntraits,byrow = T)
          diag(outDiagBlock) <- entry
          Qrow <- cbind(Qrow,outDiagBlock)
        }
        
      }
      Q <- rbind(Q,Qrow)
    }
  }
  return(Q)
}



lemad_sortingtraits <- function(traitinfo,phy){
  traitinfo <- as.matrix(traitinfo)
  if(length(phy$tip.label)!= nrow(traitinfo)){
    stop("Number of species in the tree must be the same as in the trait file")
  }
  
  if(identical(as.character(sort(phy$tip.label)),
               as.character(sort(traitinfo[,1]))) == FALSE){
    mismatch <- match( as.character(sort(traitinfo[,1])),as.character(sort(phy$tip.label)))
    mismatched <- (sort(traitinfo[,1]))[which(is.na(mismatch))]
    stop(cat("Mismatch on tip labels and taxa names, check the species:",
             mismatched))
  }
  
  traitinfo <- traitinfo[match(phy$tip.label,traitinfo[,1]),]
  traitinfo[,1] == phy$tip.label
  
  if(ncol(traitinfo) == 2){
    traits <- as.numeric(traitinfo[,2])
  }
  
  if(ncol(traitinfo)>2){
    traits <- NULL
    for (i in 1:(ncol(traitinfo)-1)){
      traits <- cbind(traits,as.numeric(traitinfo[,1 + i]))
    }
  }
  return(traits)
}



lemad_cla_id_paramPos <- function(traits,num_concealed_states){
  idparslist <- list()
  if(is.matrix(traits)){
    traits <- traits[,1]
  }
  
  ly <- length(sort(unique(traits))) * 2 * num_concealed_states
  d <- ly/2
  #idparslist[[1]] <- 1:d
  toMatrix <- 1
  matPos <- (ly + 1): (((d^2)-d) +  d * 2)
  for (i in 1:d){
    toMatrix <- c(toMatrix,matPos[(i*d-(d-1)):((i*d-(d-1)) + d)])
    
  }
  toMatrix <- toMatrix[1:d^2]
  Q <- matrix(toMatrix,ncol = d,nrow = d,byrow = T)
  diag(Q) <- NA
  lab_states <- rep(as.character(sort(unique(traits))),num_concealed_states)
  
  lab_conceal <- NULL
  for (i in 1:num_concealed_states){
    
    lab_conceal <- c(lab_conceal,rep(LETTERS[i],length(sort(unique(traits)))))
  }
  
  
  statesCombiNames <- character()
  for ( i in 1:length(lab_states)){
    statesCombiNames <- c(statesCombiNames,paste0(lab_states[i],lab_conceal[i]))
  }

  idparslist[[1]] <- matrix(0,ncol = d,nrow = 4)
  idparslist[[2]] <- (d + 1):ly
  idparslist[[3]] <- Q

  rownames(idparslist[[1]]) <- c("dual_inheritance",
                      "single_inheritance",
                      "dual_symmetric_transition",
                      "dual_asymmetric_transition")
  
  colnames(idparslist[[1]]) <- statesCombiNames
  colnames(idparslist[[3]]) <- statesCombiNames
  rownames(idparslist[[3]]) <- statesCombiNames
  names(idparslist) <- c("lambdas","mus","Q")
  names(idparslist[[2]]) <- statesCombiNames
  return(idparslist)
}


lemad_prepare_full_lambdas <- function(traits,
                                 num_concealed_states,
                                 lambd_and_modeSpe){
num_exami <- length(sort(unique(traits)))
  mat_size <- num_exami*num_concealed_states
posib_trans <- matrix(1,ncol = num_exami,nrow = num_exami,byrow = TRUE)
diag(posib_trans) <- NA
posib_trans <- lemad_q_doubletrans(traits,masterBlock = posib_trans,diff.conceal = FALSE)

full_lambdas <- list()

for(jj in 1:mat_size){
  #dual_state_inhe
  m1 <- matrix(0,ncol = mat_size,nrow = mat_size)
  m1[jj,jj] <- as.numeric(lambd_and_modeSpe[,jj][1])
  
  #single_state_inhe
  m2 <- matrix(0,ncol = mat_size,nrow = mat_size)
  m2[,jj] <- posib_trans[jj,]
  m2[jj,jj] <- 0
  m2[m2 == 1] <- as.numeric(lambd_and_modeSpe[,jj][2])
  #symet_state_emerge
  
  m3 <- matrix(0,ncol = mat_size,nrow = mat_size)
  
  diag(m3) <- posib_trans[jj,]
  m3[jj,jj] <- 0
  m3[m3 == 1] <- as.numeric(lambd_and_modeSpe[,jj][3])
  #symet_state_emerge
  
  m4 <- matrix(0,ncol = mat_size,nrow = mat_size)
  for(i in 1:length(which(posib_trans[jj,] == 1))){
    m4[which(posib_trans[jj,] == 1)[i],] <- posib_trans[jj,] 
  }
  m4[,jj] <- 0
  m4[upper.tri(m4)] <- 0
  diag(m4) <- 0
  m4[is.na(m4)] <- 0
  m4[m4 == 1] <- as.numeric(lambd_and_modeSpe[,jj][4])
  
  full_lambdas[[jj]] <- m1 + m2 + m3 + m4
}
return(full_lambdas)
}

lemad_normalize_loglik <- function(probs,loglik)
{
  sumabsprobs <- sum(abs(probs))
  probs <- probs/sumabsprobs
  loglik <- loglik + log(sumabsprobs)
  return(list(probs = probs, loglik = loglik))
}

penalty <- function(pars,loglik_penalty = 0)
{
  pars <- unlist(unlist(pars))
  return(loglik_penalty * sum(pars^2)/(2 * length(pars)))
}

### Lemad functions


#' It sorts the traits coming from simulation to be in the same order than tree tips. 
#' @title Sorting traits from simulation.
#' @param phylotree simulated phylogeny; species coded start with "t".
#' @param traits_to_sort vector with trait states.
#' @return sorted traits. 
#' @export

lemad_sort_the_traits <- function (phylotree,traits_to_sort){
  
  # there is a problem with function sort() when working with species id of the shape "t33"  
  id_species_withT <- phylotree$tip.label
  species_id_numbersOnly <- NULL
  for(i in 1:length(phylotree$tip.label)){
    
    all_species_id <- strsplit(id_species_withT,split="t")
    species_id_numbersOnly <- c(species_id_numbersOnly,as.numeric(all_species_id[[i]][2]))
  }
  species_id_new <- NULL
  for(i in 1:length(phylotree$tip.label)){
    species_id_new <- c(species_id_new, paste0("t",sort(species_id_numbersOnly)[i]))
  }
  
  table_traits_sp <- as.data.frame(cbind(species_id_new,traits_to_sort))
  
  
  sorted_traits <- table_traits_sp[order(match(table_traits_sp$species_id_new,id_species_withT)),]
  
  
  return(as.character(sorted_traits$traits_to_sort))
}


#' Prepares a q matrix for range expansion and contraction. 
#' @title Transition matrix preparation.
#' @param all_area_combination Matrix with all regions combinations coming from function prepare_full_lambdas_vicariance().
#' @param matrices_names Vector with all regions names coming from function prepare_full_lambdas_vicariance().
#' @param q_expansion Parameter ID for colonization rate.
#' @param q_contraction Parameter ID for contraction (local extinction/extirpation) rate.
#' @return the transtion matrix. 
#' @export
lemad_prepare_q_matrix <- function(all_area_combination,matrices_names,
                                   id_q_expansion,id_q_contraction){
  qs <- matrix(0,ncol=ncol(all_area_combination),nrow = ncol(all_area_combination))
  for(iji in 1:ncol(all_area_combination)){
    
    taken_sub_area <- all_area_combination[,iji]
    
   # if(iji != ncol(all_area_combination)){
      taken_sub_area <- taken_sub_area[-which(is.na(taken_sub_area))]
   # }
    
    for(ij in 1:ncol(all_area_combination)){
      vector_all_area_combination <- all_area_combination[,ij]
      
     # if(ij != ncol(all_area_combination)){
        vector_all_area_combination <-
          vector_all_area_combination[-which(is.na(vector_all_area_combination))]
        
     # }
      
      if(length(taken_sub_area) == (length(vector_all_area_combination)-1)){
        matches <- 0
        for(jj in 1:length(taken_sub_area)){
          
          if(any(taken_sub_area[jj]==vector_all_area_combination)){
            matches <- matches + 1
          }
        }
        
        if(matches == length(taken_sub_area)){
          
          state_id_complement1 <- find_stateID(vector_all_area_combination,all_area_combination)
          state_id_complement2 <- find_stateID(taken_sub_area,all_area_combination)
          qs[state_id_complement1,state_id_complement2] <- id_q_expansion
          qs[state_id_complement2,state_id_complement1] <- id_q_contraction
          
        } 
        
      }
    }
  }
  
  colnames(qs) <- matrices_names
  rownames(qs) <- matrices_names
  diag(qs) <- NA
  return(qs)
}





#' Prepares the lambda (speciation) matrices . 
#' @title Making lambda matrices.
#' @param areas Vector with single regions. It is NOT the species presence.
#' @param id_rate_insitu Parameter ID for the rate of insitu speciation.
#' @param id_rate_vicariance Parameter ID for the rate of vicariance
#' @return A list of matrices but also the names of the states (the combinations of all locations). 
#' @examples 
#' all_locations <- c("A","B","C") 
#' id_rate_insitu <- 1
#' id_rate_vicariance <- 2
#' num_max_multiregion <- 3
#' DEC_events <- FALSE
#' prepare_full_lambdas_vicariance(all_locations,num_max_multiregion,id_rate_insitu,id_rate_vicariance,DEC_events)
#' @export

prepare_full_lambdas_vicariance <- function (areas,num_max_multiregion,id_rate_insitu,id_rate_vicariance,DEC_events){
  all_matrices <- list()
  if (num_max_multiregion[1] > length(areas)) {
    stop("num_max_multiregion cannot be longer than the vector of areas")
  }
  all_area_combination <- NULL
  for (i in 1:num_max_multiregion[1]) {
    combination_from_function <- combn(areas, i)
    all_area_combination <- cbind(all_area_combination, combination_from_function)
    all_area_combination <- rbind(all_area_combination, NA)
  }
  all_area_combination
  matrices_names <- NULL
  for (i in 1:(ncol(all_area_combination))) {
    matrices_names <- c(matrices_names, paste0(all_area_combination[,i][-which(is.na(all_area_combination[, i]))], collapse = ""))
  }
  for (i in 1:ncol(all_area_combination)) {
    one_lambda <- matrix(0, ncol = ncol(all_area_combination), nrow = ncol(all_area_combination))
    state <- all_area_combination[, i]
    state_No_NA <- state[-which(is.na(state))]
    number_regions_state <- length(state_No_NA)
    if (number_regions_state == 1) {
      state_id_taken_sub_area <- find_stateID(state_No_NA,all_area_combination)
      one_lambda[state_id_taken_sub_area, state_id_taken_sub_area] <- id_rate_insitu
    }
    if (number_regions_state > 1) {
      sub_area_combination <- NULL
      for(ii in 1:(length(state_No_NA) - 1)){
        subcombination_from_function <- combn(state_No_NA,ii)
        sub_area_combination <- cbind(sub_area_combination,subcombination_from_function)
        sub_area_combination <- rbind(sub_area_combination,NA)
      }
      
      if(DEC_events == TRUE){
        # cat(state_No_NA, "\n")
        for(iii in 1:ncol(sub_area_combination)){ # splits always in single-area distributions. Vicariance
          taken_sub_area <- sub_area_combination[,iii]
          taken_sub_area <- taken_sub_area[-which(is.na(taken_sub_area))]
          if(length(taken_sub_area) == 1){#splits always in single-area distributions.
            all_complement_areas <- unique(c(taken_sub_area,state_No_NA))
            complementary_area <-  all_complement_areas[(length(taken_sub_area) + 1):length(all_complement_areas)]
            state_id_complement <- find_stateID(complementary_area,all_area_combination)
            state_id_taken_sub_area <- find_stateID(taken_sub_area,all_area_combination)
            one_lambda[state_id_taken_sub_area,state_id_complement] <- id_rate_vicariance 
          }
        }
        
        for(iii in 1:ncol(sub_area_combination)){ #  Sympatric subset
          taken_sub_area <- sub_area_combination[,iii]
          taken_sub_area <- taken_sub_area[-which(is.na(taken_sub_area))]
          if(length(taken_sub_area) == 1){#splits always in single-area distributions.
            state_id_state_No_NA <- find_stateID(state_No_NA,all_area_combination)
            state_id_taken_sub_area <- find_stateID(taken_sub_area,all_area_combination)
            one_lambda[state_id_taken_sub_area, state_id_state_No_NA] <- id_rate_insitu 
          }
        }
      } else {
        #sub_area_combination <- sub_area_combination[-nrow(sub_area_combination),]
        for(iii in 1:ncol(sub_area_combination)){
          taken_sub_area <- sub_area_combination[,iii]
          taken_sub_area <- taken_sub_area[-which(is.na(taken_sub_area))]
          all_complement_areas <- unique(c(taken_sub_area,state_No_NA))
          complementary_area <-  all_complement_areas[(length(taken_sub_area) + 1):length(all_complement_areas)]
          state_id_complement <- find_stateID(complementary_area,all_area_combination)
          state_id_taken_sub_area <- find_stateID(taken_sub_area,all_area_combination)
          one_lambda[state_id_complement,state_id_taken_sub_area] <- id_rate_vicariance 
        }
      }
    }
    colnames(one_lambda) <- matrices_names
    rownames(one_lambda) <- matrices_names
    one_lambda[lower.tri(one_lambda)] <- 0
    all_matrices[[i]] <- one_lambda
  }
  names(all_matrices) <- matrices_names
  return(list(all_matrices = all_matrices, matrices_names = matrices_names,
              all_area_combination = all_area_combination))
  
}





find_stateID <- function(region_to_compare,all_area_combination){
  
  # if(length(region_to_compare)==nrow(all_area_combination)){# it is the cosmopolitan state
  #   state_id <- ncol(all_area_combination)
  # } else {
  #   
    for( i in 1:ncol(all_area_combination)){
      the_state_No_NA <- all_area_combination[,i]
      the_state_No_NA <- the_state_No_NA[-which(is.na(the_state_No_NA))]
      if(length(the_state_No_NA)==length(region_to_compare)){
        
        if(all(sort(region_to_compare) == sort(the_state_No_NA))){
          state_id <- i
        }
        
      }
    }
  #}
  return(state_id)
  
}




#' Combines the regions into all the possible multi-region distribution. 
#' @title All possible region combination.
#' @param areas Vector with single region It is NOT the species presence.
#' @param num_max_multiregion integer indicating the maximum number of regions a lineage can possible take at a given point in time. It can go from 2 to length(areas).
#' @return A vector with all the combinations of regions 
#' @examples
#' areas <- c("A", "B", "C")
#' give_me_states_combination(areas,num_max_multiregion = 3)
#' @export

give_me_states_combination <- function(areas,num_max_multiregion){
  
  if(num_max_multiregion[1] > length(areas)){
    stop("num_max_multiregion cannot be longer than the vector of areas")
  }
  all_area_combination <- NULL
  for(i in 1:num_max_multiregion[1]){
    combination_from_function <- combn(areas,i)
    all_area_combination <- cbind(all_area_combination,combination_from_function)
    all_area_combination <- rbind(all_area_combination,NA)
    
  }
  all_area_combination
  matrices_names <- NULL
  for( i in 1: (ncol(all_area_combination))){
    matrices_names <- c(matrices_names,paste0(all_area_combination[,i][-which(is.na(all_area_combination[,i]))],
                                              collapse=""))
    
  }
  return(matrices_names)

}


#' @rawNamespace useDynLib(lemad, .registration = TRUE)
#' @rawNamespace import(Rcpp)
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
#' @keywords internal
normalize_loglik <- function(probs, loglik) {
    sumabsprobs <- sum(abs(probs))
    probs <- probs / sumabsprobs
    loglik <- loglik + log(sumabsprobs)
    cat(probs, loglik, "\n")
    return(list(probs = probs, loglik = loglik))
}

penalty <- function(pars, loglik_penalty = 0) {
    pars <- unlist(unlist(pars))
    return(loglik_penalty * sum(pars^2)/(2 * length(pars)))
}

calc_mus <- function(is_complete_tree,
                     idparslist,
                     idparsfix,
                     parsfix,
                     idparsopt,
                     initparsopt) {
    mus <- NULL
    if (is_complete_tree) {
        mus <- rep(NA, length(idparslist[[2]]))
        for (i in seq_along(idparslist[[2]])) {
            mus[i] <- c(parsfix[which(idparsfix == idparslist[[2]][i])], initparsopt[which(idparsopt == idparslist[[2]][i])])
        }
    }
    return(mus)
}