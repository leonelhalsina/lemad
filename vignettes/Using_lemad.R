## -----------------------------------------------------------------------------
rm(list = ls())
library(DDD)
library(lemad)

## -----------------------------------------------------------------------------
# amazilia_distribution <- read.csv(file="humm_distribution.csv")
data(distributioninfo)
head(amazilia_distribution)

## -----------------------------------------------------------------------------
# amazilia_tree <- read.nexus("amazilia_phylo.nex") # example of how to load a tree
data("phylo_Vign2")

## -----------------------------------------------------------------------------
rownames(amazilia_distribution) <- amazilia_distribution$species
geiger:::name.check(amazilia_tree, amazilia_distribution, data.names=NULL)
amazilia_distribution$region <- toupper(amazilia_distribution$region)
amazilia_distribution <- amazilia_distribution[order(match(amazilia_distribution$species,amazilia_tree$tip.label)),]

(amazilia_distribution$species) == (amazilia_tree$tip.label)


## -----------------------------------------------------------------------------
phylotree_recons <- amazilia_tree
species_presence <- amazilia_distribution$region

## -----------------------------------------------------------------------------
#for(ik in 1:length(species_presence)){
#  string_to_sort <- species_presence[ik]
#  species_presence[ik] <- paste(sort(unlist(strsplit(string_to_sort, ""))), collapse = "")
#}

## -----------------------------------------------------------------------------
all_areas <- c("A","B","C")

## -----------------------------------------------------------------------------
missing_spp_areas <- list()
missing_spp_areas[[1]] <- c("A","B","AB")
missing_spp_areas[[2]] <- c(0.8,0.57,0.83)

## -----------------------------------------------------------------------------
num_max_multiregion <- 3

## -----------------------------------------------------------------------------
#lineage_extinction <- "free"
# or
lineage_extinction <- 0.001

## -----------------------------------------------------------------------------
initial_lambda <- c(0.01,0.01)
initial_disperextirpation <- 0.002

## -----------------------------------------------------------------------------
DEC_events <- FALSE # to choose DIVAevents
# or when choosing DECevents:
# DEC_events <- TRUE 

## -----------------------------------------------------------------------------
object_a <- prepare_full_lambdas_vicariance(all_areas,num_max_multiregion,1,2,DEC_events)

my_dispersal_extirpationMatrix <- lemad_prepare_q_matrix(
  object_a$all_area_combination,object_a$matrices_names,4,4)
my_dispersal_extirpationMatrix

## -----------------------------------------------------------------------------
prepare_full_lambdas_vicariance(c(all_areas,"D"),num_max_multiregion,1,2,DEC_events)

## -----------------------------------------------------------------------------
# my_dispersal_extirpationMatrix[1,4] <- 0

## -----------------------------------------------------------------------------
initial_disperextirpation <- my_dispersal_extirpationMatrix

## -----------------------------------------------------------------------------
# condition_on_origin <- A

## -----------------------------------------------------------------------------
condition_on_origin <- NULL

## -----------------------------------------------------------------------------
optimizer <- "simplex"

## -----------------------------------------------------------------------------
output_vig <- lemad_analysis(
 phylotree_recons,
 species_presence,
 areas = all_areas,
 num_max_multiregion,
 condition_on_origin,
 DEC_events,
 missing_spp_areas = missing_spp_areas,
 lineage_extinction = lineage_extinction,
 initial_lambda = initial_lambda,
 initial_disperextirpation = initial_disperextirpation,
 optimizer = optimizer)


## -----------------------------------------------------------------------------
output_vig$model_ml
output_vig$number_free_pars

## -----------------------------------------------------------------------------
output_vig$estimated_rates

## -----------------------------------------------------------------------------
head(output_vig$ancestral_states)

