% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vash.R
\name{est_prior}
\alias{est_prior}
\title{Fit the mixture inverse-gamma prior of variance}
\usage{
est_prior(sehat, df, betahat, randomstart, singlecomp, unimodal, prior, g,
  maxiter, estpriormode, priormode, completeobs)
}
\arguments{
\item{sehat}{a p vector of observed standard errors}

\item{df}{appropriate degrees of freedom for (chi-square) distribution of sehat}

\item{betahat}{a p vector of estimates (optional)}

\item{randomstart}{logical, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)}

\item{singlecomp}{logical, indicating whether to use a single inverse-gamma distribution as the prior distribution for the variances}

\item{unimodal}{put unimodal constraint on the prior distribution of variances ("variance") or precisions ("precision")}

\item{prior}{string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or 1,1...,1; also can be "nullbiased" 1,1/k-1,...,1/k-1 to put more weight on first component)}

\item{g}{the prior distribution for variances (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)}

\item{maxiter}{maximum number of iterations of the EM algorithm}

\item{estpriormode}{logical, indicating whether to estimate the mode of the unimodal prior}

\item{priormode}{specified prior mode (only works when estpriormode=FALSE).}

\item{completeobs}{a p vector of non-missing flags}
}
\value{
The fitted mixture prior (g) and convergence info
}
\description{
Fit the mixture inverse-gamma prior of variance, given the variance estimates (sehat^2).
}

