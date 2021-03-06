% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priorT.R
\name{priorT}
\alias{priorT}
\title{MCMC algorithm for the estimation of local model}
\usage{
priorT(y, X, c_init, beta_init, zeta_init, nusq_init, w_init, v0, a1, a2, q,
  Nmcmc, ind)
}
\arguments{
\item{y}{binary response correspoding to \eqn{n} subjects}

\item{X}{data at time \eqn{t} - dimension \eqn{n} x \eqn{L}}

\item{c_init}{initial value for \eqn{c}}

\item{beta_init}{initial values for \eqn{\beta}}

\item{zeta_init}{initial values for \eqn{\zeta}}

\item{nusq_init}{initial values for \eqn{\nu^2}}

\item{w_init}{initial value for \eqn{w} - the complexity parameter}

\item{v0}{hyperparameter \eqn{v0}}

\item{a1}{hyperparameter \eqn{a1}}

\item{a2}{hyperparameter \eqn{a2}}

\item{q}{hyperparameter \eqn{q}}

\item{Nmcmc}{number of MCMC samples to generate}

\item{ind}{indices of MCMC samples to use after burnin ad thinning}
}
\description{
Gibbs sampling algorithm for the estiamtion of a local model
}
\examples{
priorT()
}
\keyword{priorT()}
