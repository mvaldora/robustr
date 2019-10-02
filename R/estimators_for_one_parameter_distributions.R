#' computes the psi function of the poisson MI-estimator
#' @param mu A positive real number
#' @param y A non negative integer
#' @return A real number in [0,1]
psi_pois <- function(mu, y) {
  stats::ppois(y, mu) - 0.5 * stats::dpois(y, mu) - 0.5
}

#' computes the score function of the poisson MI-estimator
#' @param mu A positive real number
#' @param mstr A vector of non negative integers
#' @return A real number

psi_pois_s <- function(mu, mstr) {
  sum(psi_pois(mu, mstr))
}


#' Computes a robust estimate of the parameter of a Poisson random variable.
#' @param x A vector of univariate observations
#' @examples
#'   x <- rpois(100,0.2)
#'   rob_estimate_poisson(x)
#'   MASS::fitdistr(x, dens="Poisson")
#' @return The estimated parameter of the Poisson distribution
#' @export
rob_estimate_poisson <- function(x) {
    if (max(x) == 0)
    est_mi <- 0
  else
    est_mi <- stats::uniroot(psi_pois_s, c(min(x), max(x)), mstr = x)$root
  est_mi
}

#' computes the psi function of the exponential MI-estimator
#' @param lambda A positive real number
#' @param y A non negative real number
#' @return A real number in [0,1]
psi_exp <- function(lambda, y) {
  stats::pexp(y, lambda) - 0.5
}


#' computes the score function of the exponential MI-estimator
#' @param lambda A positive real number
#' @param mstr A vector of non negative real numbers
#' @return A real number

psi_exp_s <- function(lambda, mstr) {
  sum(psi_exp(lambda, mstr))
}

#' Computes a robust estimate of the parameter of an Exponential random variable.
#' @param x A vector of univariate observations
#' @examples
#'   x <- rexp(50,2)
#'   rob_estimate_exp(x)
#'   MASS::fitdistr(x,dens="exponential")
#' @return The estimated parameter of the exponential distribution
#' @export
rob_estimate_exp <- function(x) {
  if (max(x) == 0)
    est_mi <- 0
  else
    est_mi <-
    stats::uniroot(psi_exp_s, c(1 / max(x), 1 / min(x)), mstr = x)$root
  est_mi
}

#' computes the psi function of the geometric MI-estimator
#' @param mu A positive real number
#' @param y A non negative integer
#' @return A real number in [0,1]
psi_geom <- function(mu, y) {
  stats::pgeom(y, mu) - 0.5 * stats::dgeom(y, mu) - 0.5
}

#' computes the score function of the geometric MI-estimator
#' @param mu A positive real number
#' @param mstr A vector of non negative integers
#' @return A real number

psi_geom_s <- function(mu, mstr) {
  sum(psi_geom(mu, mstr))
}


#' Computes a robust estimate of the parameter of a geometric random variable Y
#' @param x A vector of univariate observations
#' @examples
#'   x <- rgeom(50,1/2)
#'   rob_estimate_geom(x)
#'   MASS::fitdistr(x,dens="geometric")
#' @return The estimated parameter of the Poisson distribution
#' @export
rob_estimate_geom <- function(x) {
  if (max(x) == 0)
    est_mi <-
    0
  else
    est_mi <- stats::uniroot(psi_geom_s, c(1 / max(x), 1), mstr = x)$root
  est_mi
}

#' computes the psi function of the binomial MI-estimator
#' @param y An integer in [0,n]
#' @param n A positive integer
#' @param p A real number in (0,1)
#' @return A real number in [0,1]

psi_binom <- function(y, n, p) {
  stats::pbinom(y, n, p) - 0.5 * stats::dbinom(y, n, p) - 0.5
}

#' computes the score function of the binomial MI-estimator
#' @param mstr A  vector of integers in [0,n]
#' @param n A positive integer
#' @param p A real number in (0,1)
#' @return A real number

psi_binom_s <- function(mstr, n, p) {
  sum(psi_binom(mstr, n, p))
}


#' Computes a robust estimate of the probability of success of a binomial random variable Y with known n
#' @param x A vector of univariate observations
#' @param n An positive integer
#' @examples
#'   n<-3
#'   x <- rbinom(20,n,0.1)
#'   rob_estimate_binom(x,n)
#'   mean(x)/n
#' @return The estimated parameter of the Poisson distribution
#' @export

rob_estimate_binom <- function(x, n) {
  est_mi <- stats::uniroot(psi_binom_s, c(0, 1), mstr = x, n = n)$root
  est_mi
}
