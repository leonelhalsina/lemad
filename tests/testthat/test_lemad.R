context("lemad_test")

test_that("lemad gives the same result as Geosse", {
library(DDD)
  Sys.unsetenv("R_TESTS")
  #test geosse 
  pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
  names(pars) <- c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
  #set.seed(5)
  #phy <- diversitree::tree.geosse(pars, max.t=4, x0=0)
  phy <- NULL; rm(phy);
  utils::data('example_phy_GeoSSE', package = 'lemad');
  traits <- as.numeric(phy$tip.state)
  testit::assert(!is.null(phy))
  lik.g <- diversitree::make.geosse(phy, phy$tip.state)
  pars.g <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
  names(pars.g) <- diversitree::argnames(lik.g)
  lik.c <- diversitree::make.classe(phy, phy$tip.state+1, 3)
  pars.c <- 0 * diversitree::starting.point.classe(phy, 3)
  pars.c['lambda222'] <- pars.c['lambda112'] <- pars.g['sA']
  pars.c['lambda333'] <- pars.c['lambda113'] <- pars.g['sB']
  pars.c['lambda123'] <- pars.g['sAB']
  pars.c['mu2'] <- pars.c['q13'] <- pars.g['xA']
  pars.c['mu3'] <- pars.c['q12'] <- pars.g['xB']
  pars.c['q21'] <- pars.g['dA']
  pars.c['q31'] <- pars.g['dB']
  lik.g(pars.g) # -175.7685
  classe_diversitree_LL <- lik.c(pars.c) # -175.7685
  
  ## lemad part 
  
  traits_coded_letters <- NULL
  for(ii in 1:length(traits)){
    if(traits[ii] == 0){
      traits_coded_letters <- c(traits_coded_letters,"A")
    }
    if(traits[ii] == 1){
      traits_coded_letters <- c(traits_coded_letters,"B")
    }
    if(traits[ii] == 2){
      traits_coded_letters <- c(traits_coded_letters,"AB")
    }
  }
  
  
  lambdas<-list()
  lambdas[[1]] <- matrix(0,ncol = 3,nrow = 3,byrow = TRUE)
  #lambdas[[1]][1,1] <- 1.5
  lambdas[[1]][2,1] <- 1.5
  lambdas[[1]][3,1] <- 0.5
  lambdas[[1]][3,2] <- 1
  lambdas[[2]]<-matrix(0,ncol = 3,nrow = 3,byrow = TRUE)
  lambdas[[2]][2,2] <- 1.5
  #lambdas[[2]][2,2] <- 1.1
  lambdas[[3]] <- matrix(0,ncol=3,nrow=3,byrow=TRUE)
  lambdas[[3]][3,3] <- 0.5
  #lambdas[[3]][3,3] <- 1
  
  mus<-c(0,0.7,0.7)
  
  q<-matrix(0,ncol = 3,nrow = 3,byrow = TRUE)
  q[2,1] <- 1.4
  q[3,1] <- 1.3
  q[1,2] <- 0.7
  q[1,3] <- 0.7
  
  parameter <- list()
  parameter[[1]] <- lambdas
  parameter[[2]] <- mus
  parameter[[3]] <- q
  sampling_fraction <- rep(1,length(lambdas))
  
  num_max_multiregion <- 2
  cond <- "maddison_cond"
  root_state_weight <- "maddison_weights"
  
  setting_calculation  <- NULL
  see_ancestral_states <- FALSE
  loglik_penalty <- 0
  is_complete_tree <- FALSE
  num_threads <- 1
  
  lemad_LL <- cla_lemad_loglik(parameter, phy,  traits_coded_letters,
                               num_max_multiregion = num_max_multiregion,
                             cond = cond,
                             root_state_weight = root_state_weight,
                             sampling_fraction = sampling_fraction,
                             setting_calculation = setting_calculation,
                             see_ancestral_states = see_ancestral_states,
                             loglik_penalty = loglik_penalty,
                             is_complete_tree = is_complete_tree,
                             num_threads = num_threads,
                             method = ifelse(num_threads == 1,
                                             "odeint::bulirsch_stoer",
                                             "odeint::runge_kutta_fehlberg78"),
                             atol = 1e-16,
                             rtol = 1e-16)
  
 
  testthat::expect_equal(classe_diversitree_LL,lemad_LL)

})

test_that("lemad does a consistent ML value", {

  Sys.unsetenv("R_TESTS")
  library(DDD)
   parenthesis <- "(3:2.331354041,((2:0.02860010167,4:0.02860010167):1.434648311,1:1.463248413):0.8681056285);"
  phylotree_recons <- ape::read.tree(file="",parenthesis)
  species_presence <- c("B","C","A","AC")
  areas <- c("A","B","C")
  missing_spp_areas <- list()
  missing_spp_areas[[1]] <- c("A","AC")
  missing_spp_areas[[2]] <- c(0.8,0.95)
  
  output1 <- lemad_analysis(
    phylotree_recons,
    species_presence,
    areas,
    num_max_multiregion = 3,
	condition_on_origin = NULL,
  DEC_events = TRUE,
    missing_spp_areas,
    lineage_extinction = 0.005,
    initial_lambda = NULL,
    initial_disperextirpation = NULL)
  
  cat("ML",round(output1$model_ml,digits = 4),"\n")
 testthat::expect_equal(round(output1$model_ml,digits = 4),-10.0905)
  #testthat::expect_equal(-2,-2)
})