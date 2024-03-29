% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lemad_utils.R
\name{lemad_sort_the_traits}
\alias{lemad_sort_the_traits}
\title{Sorting traits from simulation.}
\usage{
lemad_sort_the_traits(phylotree, traits_to_sort)
}
\arguments{
\item{phylotree}{simulated phylogeny; species coded start with "t".}

\item{traits_to_sort}{vector with trait states.}
}
\value{
sorted traits.
}
\description{
It sorts the traits coming from simulation to be in the same order than tree tips.
}
