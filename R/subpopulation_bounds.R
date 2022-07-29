#' Plot rule benefit bounds for subpopulation containing all treatment effect heterogeneity.
#'
#' @param bounded_outcome logical; `TRUE` if the outcome is bounded and `FALSE` otherwise,
#'     default is `FALSE`
#' @param conf.upper logical; whether a upper confidence bound should be provided
#' @param mean_1_s estimated treatment arm mean in the subpopulation
#' @param mean_0_s estimated control arm mean in the subpopulation
#' @param s_prop_low lowest relative size of the subpopulation to plot, default is 0.1
#' @param s_prop_high highest relative size of the subpopulation to plot, default is 0.9
#' @inheritParams teh_closed_form_bounded
#'
#' @returns
#' A list with two elements:
#'
#' * `estimate` the estimated benefit of the ideal treatment rule over the best uniform treatment
#' * `conf.upper` The `level` upper confidence bound if `conf.upper` is TRUE, `NULL` otherwise
#'
#' @examples
#'
#' #outcome is 1-5, variance 1 in each arm. treatment effect is 1
#' #   100 people in each arm
#' #   subpopulation has the same effect as the overall population
#' subpop_ideal_rule_benefits(s2_1 = 1,
#'                            s2_0 = 1,
#'                            n_1 = 100,
#'                            n_0 = 100,
#'                            mean_1 = 3,
#'                            mean_0 = 2,
#'                            mean_1_s = 3,
#'                            mean_0_s = 2,
#'                            m = 1,
#'                            M = 5,
#'                            bounded_outcome = TRUE)
#'
#'
#' #outcome is 1-5, variance 1 in each arm. treatment effect is 1, 100 people in each arm
#' #   100 people in each arm
#' #   subpopulation has no effect
#' subpop_ideal_rule_benefits(s2_1 = 1,
#'                            s2_0 = 1,
#'                            n_1 = 100,
#'                            n_0 = 100,
#'                            mean_1 = 3,
#'                            mean_0 = 2,
#'                            mean_1_s = 2,
#'                            mean_0_s = 2,
#'                            m = 1,
#'                            M = 5,
#'                            bounded_outcome = TRUE)
#'
#'
#'
#' @export
subpop_ideal_rule_benefits <- function(s2_1,
                                       s2_0,
                                       n_1 = NULL,
                                       n_0 = NULL,
                                       mean_1_s,
                                       mean_0_s,
                                       mean_1,
                                       mean_0,
                                       s_prop_low = 0.1,
                                       s_prop_high = 0.9,
                                       m = NULL,
                                       M = NULL,
                                       level = 0.95,
                                       bounded_outcome = FALSE,
                                       conf.upper = FALSE)
{
  #get x values
  s_prop <- seq(s_prop_low, s_prop_high, by = 0.05)

  mean_1_not_s <- (mean_1 - mean_1_s * s_prop) / (1 - s_prop)

  mean_0_not_s <- (mean_0 - mean_0_s * s_prop) / (1 - s_prop)

  delta_not_s <- mean_1_not_s - mean_0_not_s

  delta_s <- mean_1_s - mean_0_s

  delta <- mean_1 - mean_0

  n <- n_1 + n_0

  r_1 <- n_1 / n

  r_0 <- n_0 / n

  nu_hat <- sqrt(s2_1 / s2_0)

  max_sqdev_1 <- max((M - mean_1) ^ 2, (m - mean_1) ^ 2)

  max_sqdev_0 <- max((M - mean_0) ^ 2, (m - mean_0) ^ 2)


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

  estimates <-
    0.5 * sqrt(-1 * (1 - s_prop) * (delta_not_s - delta_s) ^ 2 + teh_bound /
                 s_prop + delta_s ^ 2)

  var <- (16 * (s2_0 * (nu_hat + 1) ^ 2 / s_prop + delta ^ 2)) ^ (-1) *
    (
      (max_sqdev_1 * s2_1 - s2_1 ^ 2) * (nu_hat ^ (-1) + 1) ^ 2 / (r_1 * s_prop ^
                                                                     2) +
        (max_sqdev_0 * s2_0 - s2_0 ^ 2) * (nu_hat + 1) ^ 2 / (r_0 * s_prop ^
                                                                2) +
        4 * (s2_1 / r_1 + s2_0 / r_0) * delta ^ 2 +
        (max(abs(M), abs(m)) * (s2_1 + mean_1 ^ 2) - 3 * mean_1 * s2_1 -
           mean_1 ^ 3) * 2 *
        delta * (1 + nu_hat ^ (-1)) / (r_1 * s_prop) +
        (max(abs(M), abs(m)) * (s2_0 + mean_0 ^ 2) + 3 * mean_0 * s2_0 +
           mean_0 ^ 3) * 2 *
        delta * (1 + nu_hat) / (r_0 * s_prop)
    )

  q_high = 0.5 + level / 2

  ci.upper <-
    stats::qnorm(q_high, mean = estimates, sd = sqrt(var / n))

  return(list(estimates = estimates, conf.upper = ci.upper))

}

#' Plot rule benefit bounds for subpopulation containing all treatment effect heterogeneity.
#'
#' @inheritParams subpop_ideal_rule_benefits
#'
#' @returns
#' A list with two elements:
#'
#' * `estimate` the estimated benefit of the ideal treatment rule over the best uniform treatment
#' * `conf.upper` The `level` upper confidence bound if `conf.upper` is TRUE, `NULL` otherwise
#' @importFrom rlang .data
#'
#' @examples
#'
#' #outcome is 1-5, variance 1 in each arm. treatment effect is 1,
#' #   100 people in each arm
#' #   subpopulation has the same effect as the overall population
#' subpop_ideal_rule_benefits(s2_1 = 1,
#'                            s2_0 = 1,
#'                            n_1 = 100,
#'                            n_0 = 100,
#'                            mean_1 = 3,
#'                            mean_0 = 2,
#'                            mean_1_s = 3,
#'                            mean_0_s = 2,
#'                            m = 1,
#'                            M = 5,
#'                            bounded_outcome = TRUE)
#'
#'
#' #outcome is 1-5, variance 1 in each arm. treatment effect is 1,
#' #   100 people in each arm
#' #   subpopulation has no effect
#' subpop_ideal_rule_benefits(s2_1 = 1,
#'                            s2_0 = 1,
#'                            n_1 = 100,
#'                            n_0 = 100,
#'                            mean_1 = 3,
#'                            mean_0 = 2,
#'                            mean_1_s = 2,
#'                            mean_0_s = 2,
#'                            m = 1,
#'                            M = 5,
#'                            bounded_outcome = TRUE)
#' @export
plot_ideal_rule_benefits <- function(s2_1,
                                     s2_0,
                                     n_1 = NULL,
                                     n_0 = NULL,
                                     mean_1_s,
                                     mean_0_s,
                                     mean_1,
                                     mean_0,
                                     s_prop_low = 0.1,
                                     s_prop_high = 0.9,
                                     m = NULL,
                                     M = NULL,
                                     level = 0.95,
                                     bounded_outcome = F,
                                     conf.upper = F) {
  temp <- subpop_ideal_rule_benefits(
    s2_1 = s2_1,
    s2_0 = s2_0,
    n_1 = n_1,
    n_0 = n_0,
    mean_1_s = mean_1_s,
    mean_0_s = mean_0_s,
    mean_1 = mean_1,
    mean_0 = mean_0,
    s_prop_low = s_prop_low,
    s_prop_high = s_prop_high,
    m = m,
    M = M,
    level = level,
    bounded_outcome = bounded_outcome,
    conf.upper = conf.upper
  )

  #get x values
  s_prop <- seq(s_prop_low, s_prop_high, by = 0.05)

  dat <- data.frame(
    Benefit = c(temp$estimate, temp$conf.upper),
    Type = c(rep("Estimate", length(s_prop)), rep("Confidence Bound", length(s_prop))),
    p = rep(s_prop, 2)
  )

  print(
    ggplot2::ggplot(
      data = dat,
      ggplot2::aes(
        x = .data$p,
        y = .data$Benefit,
        linetype = .data$Type
      )
    ) +
      ggplot2::scale_linetype_manual(
        breaks = c("Estimate", "Confidence Bound"),
        values = c(1, 2)
      ) +
      ggplot2::geom_line(size = 1) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(x = "Size of Sub-population",
                    y = "Maximum Benefit",
                    linetype = "")
  )


}
