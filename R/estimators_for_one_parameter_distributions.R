#' Computes a robust estimate of the parameter of a Poisson random variable.
#' @param x A vector of univariate observations
#' @examples
#'   x <- rpois(100,0.2)
#'   rob_estimate_poisson(x)
#'   fitdistr(x, dens="Poisson")
#' @return The estimated parameter of the Poisson distribution
#' @importFrom("stats", "uniroot")
#' @export
rob_estimate_poisson <- function(x) {
  phi_pois <- function(mu, y) {
    ppois(y, mu) - 0.5 * dpois(y, mu) - 0.5
  }
  phi_pois_s <- function(mu, mstr) {
    sum(phi_pois(mu, mstr))
  }
  if (max(x) == 0)
    est_mi <- 0
  else
    est_mi <- stats::uniroot(phi_pois_s, c(min(x), max(x)), mstr = x)$root
  est_mi
}

#' Computes a robust estimate of the parameter of an Exponential random variable.
#' @param y A vector of univariate observations
#' @examples
#'   x <- rexp(50,2)
#'   rob_estimate_exp(x)
#'   fitdistr(x,dens="exponential")
#' @return The estimated parameter of the exponential distribution
#' @export
rob_estimate_exp <- function(x) {
  phi_exp <- function(lambda, y) {
    pexp(y, lambda) - 0.5
  }
  phi_exp_s <- function(lambda, mstr) {
    sum(phi_exp(lambda, mstr))
  }
  if (max(x) == 0)
    est_mi <- 0
  else
    est_mi <-
    stats::uniroot(phi_exp_s, c(1 / max(x), 1 / min(x)), mstr = x)$root
  est_mi
}
#' Computes a robust estimate of the parameter of a geometric random variable Y
#' @param x A vector of univariate observations
#' @examples
#'   x <- rgeom(50,1/2)
#'   rob_estimate_geom(x)
#'   fitdistr(x,dens="geometric")
#' @return The estimated parameter of the Poisson distribution
#' @importFrom("stats", "uniroot")
#' @export
rob_estimate_geom <- function(x) {
  phi_geom <- function(mu, y) {
    pgeom(y, mu) - 0.5 * dgeom(y, mu) - 0.5
  }
  phi_geom_s <- function(mu, mstr) {
    sum(phi_geom(mu, mstr))
  }
  if (max(x) == 0)
    est_mi <-
    0
  else
    est_mi <- stats::uniroot(phi_geom_s, c(1 / max(x), 1), mstr = x)$root
  est_mi
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
  phi_binom <- function(y, n, p) {
    pbinom(y, n, p) - 0.5 * dbinom(y, n, p) - 0.5
  }
  phi_binom_s <- function(mstr, n, p) {
    sum(phi_binom(mstr, n, p))
  }
  est_mi <- stats::uniroot(phi_binom_s, c(0, 1), mstr = x, n = n)$root
  est_mi
}
#' Computes a robust estimate of the number of trials of a binomial random variable Y with known p
#' @param x A vector of univariate observations
#' @param p A real number in (0,1)
#' @examples
#'   p<-0.2
#'   x <- rbinom(100,3,p)
#'   rob_estimate_binom_for_n(x,p)
#' @return The estimated number of trials
#' @export
rob_estimate_binom_for_n <- function(x, p) {
  phi_binom <- function(y, n, p) {
    pbinom(y, n, p) - 0.5 * dbinom(y, n, p) - 0.5
  }
  phi_binom_s <- function(n) {
    sum(phi_binom(x, n, p))
  }
  phi_binom_s <- Vectorize(phi_binom_s)
  if (phi_binom_s(max(x)) < 0)
    est_mi <- 1
  else{
    i <- 0
    while ( (phi_binom_s(max(x) + i) * phi_binom_s(max(x) + i + 1)) > 0)
      i <- i + 1
    minimo_abs <-
      min(abs(phi_binom_s(max(x) + i)), abs(phi_binom_s(max(x) + i + 1)))
    indice_minimo <-
      which(abs(c(
        phi_binom_s(max(x) + i), phi_binom_s(max(x) + i + 1)
      )) == minimo_abs)
    est_mi <- max(x) + i + indice_minimo - 1
  }
  est_mi
}
