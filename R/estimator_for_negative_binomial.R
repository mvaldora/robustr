#' Computes the first entry of the score function of a MI-estimator for negative binomial distribution.

#' @param x A vector of univariate observations
#' @param mu A possitive real number
#' @param alpha A possitive real number
#' @examples
#'   mu <- 2
#'   alpha <- 1
#'   x <- rnbinom(100,mu=mu,size=1/alpha)
#'   mean(psibinneg1(mu,alpha,x))
#' @return A real number

psibinneg1 <- function( mu, alpha, x) {
  facu <- stats::pnbinom( x, mu = mu, size = 1 / alpha)
  fpun <- stats::dnbinom( x, mu = mu, size = 1 / alpha)
  mean( facu - 1 / 2 * fpun) - 1 / 2
}

#' Computes the  second entry of the score function of a MI-estimator for negative binomial distribution.
#' @param x A vector of univariate observations
#' @param mu A possitive real number
#' @param alpha A possitive real number
#' @examples
#'   mu <- 2
#'   alpha <- 1
#'   x <- rnbinom(100,mu=mu,size=1/alpha)
#'   mean(psibinneg2(mu,alpha,x))
#' @return A real number


psibinneg2 <- function(mu, alpha, x) {
  facu <- stats::pnbinom(x, mu = mu, size = 1 / alpha)
  fpun <- stats::dnbinom(x, mu = mu, size = 1 / alpha)
  mean(facu ^ 2 + 1 / 3 * fpun ^ 2 - facu * fpun) - 1 / 3
}


#' Computes the score function of a MI-estimator for negative binomial distribution.
#' @param par A two dimensional vector of parameters for negative binomial distribution. The first entry is the mean and the second is the dispertion parameter, that equals 1/size.
#' @param x A vector of univariate observations
#' @examples
#'   mu <- 2
#'   alpha <- 1
#'   x <- rnbinom(100,mu=mu,size=1/alpha)
#'   mean(psibinneg2(mu,alpha,x))
#' @return A real number

score_mi_negbin <- function(par, x) {
  mu <- par[1]
  alpha <- par[2]
  (psibinneg1(mu, alpha, x)) ^ 2 + (psibinneg2(mu, alpha, x)) ^ 2
}

#' Computes a robust estimate of the parameters of a negative binomial random variable.
#' @param x A vector of univariate observations
#' @param start Initial vector of parameters to start the optimization
#' @examples
#'   library(MASS)
#'   x <- rnbinom(100,mu=2,size=1)
#'   MASS::fitdistr(x,dens="negative binomial")
#'   rob_estimate_negbin(x)
#' @return The estimated vector of parameters. The first entry is the mean and the second is the dispersion parameter, that equals 1/size.

rob_estimate_negbin <- function(x, start = NULL) {
  if (is.null(start)) {
    m_estim <- MASS::rlm(x ~ 1)
    startmu <- max(m_estim$coefficients, 0.1)
    startalpha <-
      max( (summary(m_estim)$sigma ^ 2 - startmu) / startmu ^ 2, 0.01)
  } else{
    startmu <- start[1]
    startalpha <- start[2]
  }
  pars_mi_2 <- stats::optim(par = c(startmu, startalpha),
                   score_mi_negbin,
                   x = x)
  pars_mi_2$par
}
