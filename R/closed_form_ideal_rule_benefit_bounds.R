#' General Closed Form Ideal Rule Benefit Bound
#'
#' @param bounded_outcome logical; `TRUE` if the outcome is bounded and `FALSE` otherwise,
#'     default is `FALSE`
#' @param conf.upper logical; whether a upper confidence bound should be provided
#' @inheritParams teh_closed_form_bounded
#'
#' @returns
#' A list with two elements:
#'
#' * `estimate` the estimated benefit of the ideal treatment rule over the best uniform treatment
#' * `conf.upper` The `level` upper confidence bound if `conf.upper` is TRUE, `NULL` otherwise
#'
#' @examples #tbd
#'
#' @export
closed_form_ideal_rule_benefit <- function(s2_1,
                                           s2_0,
                                           n_1 = NULL,
                                           n_0 = NULL,
                                           mean_1 = NULL,
                                           mean_0 = NULL,
                                           m = NULL,
                                           M = NULL,
                                           level = 0.95,
                                           bounded_outcome = F,
                                           conf.upper = F)
{
  if (bounded_outcome) {
    if (is.null(m) | is.null(M)){
      stop("For bounded outcomes, bounds must be given by the m and M arguments")
    }

    if (m >= M) stop("For bounded outcomes, m must be less than M")

  } else{
    if (conf.int)
      warning("No confidence interval can be provided for unbounded outcomes")
  }

  teh_bound <- treatment_effect_heterogeneity_bound(s2_1 = s2_1,
                                                    s2_0 = s2_0,
                                                    n_1 = n_1,
                                                    n_0 = n_0,
                                                    mean_1 = mean_1,
                                                    mean_0 = mean_0,
                                                    m = m,
                                                    M = M,
                                                    bounded_outcome = bounded_outcome,
                                                    conf.int = FALSE)$estimates[2]

  estimate <- 0.5 * sqrt( teh_bound + (mean_1 - mean_0)^2 )

  if(conf.upper){

    n <- n_1 + n_0

    r_1 <- n_1/n

    r_0 <- n_0/n

    nu_hat <- sqrt(s2_1/s2_0)

    mean_delta <- mean_1 - mean_0

    max_sqdev_1 <- max((M-mean_1)^2,(m-mean_1)^2)

    max_sqdev_0 <- max((M-mean_0)^2,(m-mean_0)^2)

    point <- 1/2*sqrt(s2_0*(nu_hat+1)^2 +mean_delta^2)

    var <- (16*(s2_0*(nu_hat+1)^2+mean_delta^2))^(-1)*
      (
        (max_sqdev_1*s2_1-s2_1^2)*(nu_hat^(-1)+ 1)^2/r_1 +
          (max_sqdev_0*s2_0-s2_0^2)*(nu_hat+ 1)^2/r_0 +
          4*(s2_1/r_1+s2_0/r_0)*mean_delta^2+
          (max(abs(M),abs(m))*(s2_1+mean_1^2)-3*mean_1*s2_1-mean_1^3)*2*mean_delta*(1+nu_hat^(-1))/r_1+
          (max(abs(M),abs(m))*(s2_0+mean_0^2)+3*mean_0*s2_0+mean_0^3)*2*mean_delta*(1+nu_hat)/r_0
      )

    ci.upper <- stats::qnorm(level, mean = point, sd = sqrt(var/n))

  } else {ci.upper <- NULL}

  return(list(estimate = estimate, conf.upper = ci.upper))

}
