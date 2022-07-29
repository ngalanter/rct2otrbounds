


#' Closed Form Treatment Effect Heterogeneity Bound for Bounded Outcomes
#'
#' @param m,M two real numbers defining the outcome range with `m`<`M`
#' @param s2_1 estimated outcome variance of the treatment group
#' @param s2_0 estimated outcome variance of the control group
#' @param n_1 sample size of the treatment group
#' @param n_0 sample size of the control group
#' @param mean_1 mean outcome of the treatment group
#' @param mean_0 mean outcome of the control group
#' @param conf.int whether a confidence interval is generated
#' @param level confidence level between 0 and 1 for confidence interval;
#'     default is 0.95
#'
#' @returns
#' A list with two elements:
#'
#' * `estimates` a two element vector of the estimated lower and upper bounds on
#'     treatment effect heterogeneity
#' * `conf.int` a two element vector of the `level` confidence interval around the bounds,
#'     for bounded outcomes when the `conf.int` argument is set to TRUE
#' @examples #internal function
#'
teh_closed_form_bounded <-
  function(m,
           M,
           s2_1,
           s2_0,
           n_1 = NULL,
           n_0 = NULL,
           mean_1 = NULL,
           mean_0 = NULL,
           conf.int = FALSE,
           level = 0.95)
  {
    #calculating the upper and lower bounds on treatment effect heterogeneity for bounded outcomes
    point_low <- s2_1 + s2_0 - 2 * min(sqrt(s2_1) * sqrt(s2_0),
                                       (mean_1 - m) * (M - mean_0),
                                       (M - mean_1) * (mean_0 - m))

    point_high <- s2_1 + s2_0 + 2 * min(sqrt(s2_1) * sqrt(s2_0),
                                        (M - mean_1) * (M - mean_0),
                                        (mean_1 - m) * (mean_0 - m))


    if (conf.int) {
      if (any(sapply(list(n_1, n_0, mean_1, mean_0), is.null))) {
        stop("Values for n_1,n_0,mean_1,and mean_0 are required to form confidence intervals.")
      }

      q_low = 0.5 - level / 2
      q_high = 0.5 + level / 2

      #total sample size
      n <- n_1 + n_0

      #proportion of sample size in each trial arm
      r_1 <- n_1 / n
      r_0 <- n_0 / n

      #sample variance ratio of treatment to control arm
      nu_hat <- sqrt(s2_1 / s2_0)

      #upper bound on the variance of general lower bound on teh
      var_low <-
        (max((m - mean_1) ^ 2, (M - mean_1) ^ 2) * s2_1 - s2_1 ^ 2) * (1 / nu_hat - 1) ^
        2 / r_1 +
        (max((m - mean_0) ^ 2, (M - mean_0) ^ 2) * s2_0 - s2_0 ^ 2) * (nu_hat - 1) ^
        2 / r_0

      #upper bound on the variance of general upper bound on teh
      var_high <-
        (max((m - mean_1) ^ 2, (M - mean_1) ^ 2) * s2_1 - s2_1 ^ 2) * (1 / nu_hat + 1) ^
        2 / r_1 +
        (max((m - mean_0) ^ 2, (M - mean_0) ^ 2) * s2_0 - s2_0 ^ 2) * (nu_hat + 1) ^
        2 / r_0

      #lower confidence interval limit for the bounds
      ci_low <-
        stats::qnorm(q_low, mean = point_low, sd = sqrt(var_low / n))

      #upper confidence interval limit for the bounds
      ci_high <-
        stats::qnorm(q_high, mean = point_high, sd = sqrt(var_high / n))
    } else {
      ci_low = NULL
      ci_high = NULL
    }
    return(list(
      estimates = c(point_low, point_high),
      conf.int = c(ci_low, ci_high)
    ))



  }

#' Closed Form Treatment Effect Heterogeneity Bound for Unbounded Outcomes
#'
#' @param s2_1 estimated outcome variance of the treatment group
#' @param s2_0 estimated outcome variance of the control group
#'
#' @returns
#' A a vector with the lower bound on treatment effect variance followed by the upper bound
#'
#' @examples #internal function
#'
teh_closed_form_unbounded <- function(s2_1, s2_0) {
  #sample variance ratio of treatment to control arm
  nu_hat <- sqrt(s2_1 / s2_0)

  #calculating the upper and lower bounds on treatment effect heterogeneity for unbounded outcomes
  low <- s2_0 * (nu_hat - 1) ^ 2

  high <- s2_0 * (nu_hat + 1) ^ 2


  return(c(low, high))

}

#' General Closed Form Treatment Effect Heterogeneity Bound
#'
#' @param bounded_outcome logical; `TRUE` if the outcome is bounded and `FALSE` otherwise,
#'     default is `FALSE`
#' @inheritParams teh_closed_form_bounded
#'
#' @returns
#' A list with two elements:
#'
#' * `estimates` a two element vector of the estimated lower and upper bounds on
#'     treatment effect heterogeneity
#' * `conf.int` a two element vector of the `level` confidence interval around the bounds,
#'     for bounded outcomes when the `conf.int` argument is set to TRUE
#'
#' @examples
#'
#' #outcome is 1-5, variance 1 in each arm. treatment effect is 1
#' #   100 people in each arm
#' treatment_effect_heterogeneity_bound(s2_1 = 1,
#'                                      s2_0 = 1,
#'                                      n_1 = 100,
#'                                      n_0 = 100,
#'                                      mean_1 = 3,
#'                                      mean_0 = 2,
#'                                      m = 1,
#'                                      M = 5,
#'                                      conf.int = TRUE,
#'                                      bounded_outcome = TRUE)
#'
#' #outcome is positive, variance 1 in each arm
#' treatment_effect_heterogeneity_bound(s2_1 = 1,
#'                                      s2_0 = 1)
#'
#' @export
treatment_effect_heterogeneity_bound <- function(s2_1,
                                                 s2_0,
                                                 n_1 = NULL,
                                                 n_0 = NULL,
                                                 mean_1 = NULL,
                                                 mean_0 = NULL,
                                                 m = NULL,
                                                 M = NULL,
                                                 level = 0.95,
                                                 bounded_outcome = FALSE,
                                                 conf.int = FALSE)
{
  if (bounded_outcome) {
    if (bounded_outcome & (is.null(m) | is.null(M))) {
      stop("For bounded outcomes, bounds must be given by the m and M arguments")
    }

    if (bounded_outcome & (m >= M)) {
      stop("For bounded outcomes, m must be less than M")
    }

    return(
      teh_closed_form_bounded(
        m = m,
        M = M,
        s2_1 = s2_1,
        s2_0 = s2_0,
        n_1 = n_1,
        n_0 = n_0,
        mean_1 = mean_1,
        mean_0 = mean_0,
        conf.int = conf.int,
        level = level
      )
    )
  } else{
    if (conf.int)
      warning("No confidence interval can be provided for unbounded outcomes")

    return(list(
      estimates = teh_closed_form_unbounded(s2_1, s2_0),
      conf.int = c(NULL, NULL)
    ))
  }

}
