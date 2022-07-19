#' Plot rule benefit bounds for subpopulation containing all treatment effect heterogeneity.
#'
#' @param bounded_outcome logical; `TRUE` if the outcome is bounded and `FALSE` otherwise,
#'     default is `FALSE`
#' @param conf.upper logical; whether a upper confidence bound should be provided
#' @param mean_1_s estimated treatment arm mean in the subpopulation
#' @param mean_0_s estimated control arm mean in the subpopulation
#' @param s_prop_low lowest relative size of the subpopulation to plot, default is 0.1
#' @param s_prop_low highest relative size of the subpopulation to plot, default is 0.9
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
plot_subpop_ideal_rule_benefit <- function(s2_1,
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
                                           conf.upper = F)
{

  #tbc

  2 + 2

}
