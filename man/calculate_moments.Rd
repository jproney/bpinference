% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_moments.R
\name{calculate_moments}
\alias{calculate_moments}
\title{Computes the first and second moments of the the joint population-size distribution of a Markov branching proccess}
\usage{
calculate_moments(e_mat, p_vec, r_vec, z0_vec, tf)
}
\arguments{
\item{e_mat}{The matrix of birth events that can occur in the brnaching process. Dimensions \code{nevents} x \code{ntypes}}

\item{p_vec}{A vector containing the parent type for each of the birth events in \code{e_mat}. Dimensions \code{nevents} x 1}

\item{r_vec}{A vector containing the rate at which each of the birth events in \code{e_mat} occurs. Dimensions \code{nevents} x 1}

\item{z0_vec}{The initial population vector at time 0. Dimensions \code{ntypes} x 1}
}
\value{
A list of two elements: \code{mu_mat} contains the mean number of each type at time \code{tf}, \code{sigma_mat} is the covariance matrix.
}
\description{
Computes the first and second moments of the the joint population-size distribution of a Markov branching proccess
}
