#' General Closed Form Ideal Rule Benefit Bound
#'
#' @param bounded_outcome logical; `TRUE` if the outcome is bounded and `FALSE` otherwise,
#'     default is `FALSE`
#' @param conf.upper logical; whether a upper confidence bound should be provided
#' @param binary_outcome logical; whether the outcome is binary
#' @inheritParams teh_closed_form_bounded
#'
#' @returns
#' A list with two elements:
#'
#' * `estimate` the estimated benefit of the ideal treatment rule over the best uniform treatment
#' * `conf.upper` The `level` upper confidence bound if `conf.upper` is TRUE, `NULL` otherwise
#'
#' @examples
#' #outcome is 1-5, variance 1 in each arm. treatment effect is 1
#' #   100 people in each arm
#' closed_form_ideal_rule_benefit(s2_1 = 1,
#'                                s2_0 = 1,
#'                                n_1 = 100,
#'                                n_0 = 100,
#'                                mean_1 = 3,
#'                                mean_0 = 2,
#'                                m = 1,
#'                                M = 5,
#'                                bounded_outcome = TRUE,
#'                                conf.upper = TRUE)
#'
#' #outcome is positive, variance 1 in each arm. treatment effect is 1
#' #   100 people in each arm
#' closed_form_ideal_rule_benefit(s2_1 = 1,
#'                                s2_0 = 1,
#'                                n_1 = 100,
#'                                n_0 = 100,
#'                                mean_1 = 3,
#'                                mean_0 = 2)
#'
#' #outcome is binary, treatment effect is 0.3
#' #   100 people in each arm
#' closed_form_ideal_rule_benefit(n_1 = 100,
#'                                n_0 = 100,
#'                                mean_1 = 0.4,
#'                                mean_0 = 0.7,
#'                                conf.upper = TRUE,
#'                                binary_outcome = TRUE)
#'
#' @export
closed_form_ideal_rule_benefit <- function(s2_1 = NULL,
                                           s2_0 = NULL,
                                           n_1 = NULL,
                                           n_0 = NULL,
                                           mean_1,
                                           mean_0,
                                           m = NULL,
                                           M = NULL,
                                           level = 0.95,
                                           bounded_outcome = FALSE,
                                           binary_outcome = FALSE,
                                           conf.upper = FALSE)
{
  #if outcome is binary return binary bound
  if (binary_outcome) {
    if (mean_1 > mean_0) {

      estimate <- min(mean_0, 1 - mean_1)

    } else{
      estimate <- min(mean_1, 1 - mean_0)

    }

  } else { #otherwise return bounded or general bound

    if(is.null(s2_1) | is.null(s2_0)){
      stop("For non-binary outcomes, must provide variance in each arm.")
    }

    if (bounded_outcome) {
      if (is.null(m) | is.null(M)) {
        stop("For bounded outcomes, bounds must be given by the m and M arguments")
      }

      if (m >= M)
        stop("For bounded outcomes, m must be less than M")

    } else{
      if (conf.upper)
        warning("No confidence interval can be provided for unbounded outcomes")
    }

    #get upper bound on treatment effect heterogeneity
    teh_bound <- treatment_effect_heterogeneity_bound(
      s2_1 = s2_1,
      s2_0 = s2_0,
      n_1 = n_1,
      n_0 = n_0,
      mean_1 = mean_1,
      mean_0 = mean_0,
      m = m,
      M = M,
      bounded_outcome = bounded_outcome,
      conf.int = FALSE
    )$estimates[2]

    #get point estimate of upper bound on ideal treatment rule benefit
    estimate <- 0.5 * sqrt(teh_bound + (mean_1 - mean_0) ^ 2)
  }
  if (conf.upper & (bounded_outcome | binary_outcome)) {
    #if the outcome is binary, set bounds at 0 and 1
    if (binary_outcome) {
      m <- 0
      M <- 1
    }

    #total sample size
    n <- n_1 + n_0

    #proportion of sample size in treatment arm
    r_1 <- n_1 / n

    #proportion of sample size in control arm
    r_0 <- n_0 / n

    #sample standard deviation ratio
    nu_hat <- sqrt(s2_1 / s2_0)

    #sample treatment effect
    mean_delta <- mean_1 - mean_0

    #bounds on variance of squared sd estimates in each arm
    max_sqdev_1 <- max((M - mean_1) ^ 2, (m - mean_1) ^ 2)
    max_sqdev_0 <- max((M - mean_0) ^ 2, (m - mean_0) ^ 2)

    #point estimate based on general bound
    point <- 1 / 2 * sqrt(s2_0 * (nu_hat + 1) ^ 2 + mean_delta ^ 2)

    #bound on the variance of the point estimate
    var <-
      (16 * (s2_0 * (nu_hat + 1) ^ 2 + mean_delta ^ 2)) ^ (-1) *
      (
        (max_sqdev_1 * s2_1 - s2_1 ^ 2) * (nu_hat ^ (-1) + 1) ^ 2 / r_1 +
          (max_sqdev_0 * s2_0 - s2_0 ^ 2) * (nu_hat + 1) ^ 2 / r_0 +
          4 * (s2_1 / r_1 + s2_0 / r_0) * mean_delta ^ 2 +
          (max(abs(M), abs(m)) * (s2_1 + mean_1 ^ 2) - 3 * mean_1 * s2_1 -
             mean_1 ^ 3) * 2 * mean_delta * (1 + nu_hat ^ (-1)) / r_1 +
          (max(abs(M), abs(m)) * (s2_0 + mean_0 ^ 2) + 3 * mean_0 * s2_0 +
             mean_0 ^ 3) * 2 * mean_delta * (1 + nu_hat) / r_0
      )

    q_high = 0.5 + level / 2

    #upper confidence bound at the specified level
    ci.upper <-
      stats::qnorm(q_high, mean = point, sd = sqrt(var / n))

  } else {
    ci.upper <- NULL
  }

  return(list(estimate = estimate, conf.upper = ci.upper))

}
