#' Helper Function to Transform Outcome
#'
#'Transforms outcome range and means. First, if original scale means lower outcomes are beneficial,
#'     outcome scale is multiplied by -1. Then, once the scale is such that higher outcomes are beneficial,
#'     the range is shifted to start at 0.
#'     The intervals between possible values in the outcome range are not affected.
#'
#' @param y_range the range of outcomes
#' @param mean_worse mean of the worse treatment arm
#' @param mean_better mean of the better treatment arm
#' @param scale whether "higher" or "lower" outcomes are better
#'
#' @return a list of the new range, new mean of the better arm, and new mean of the worse arm
#'
#' @examples #this is an internal function
transform_range <- function(y_range,
                            mean_worse,
                            mean_better,
                            scale) {
  #flip outcome scale if lower outcomes are better, transform range to start at 0
  if (scale == "lower") {
    new_y_range <- -1 * y_range + max(y_range)
    new_mean_worse <- max(y_range) - mean_worse
    new_mean_better <- max(y_range) - mean_better

  } else{
    new_y_range <- y_range - min(y_range)
    new_mean_worse <- mean_worse - min(y_range)
    new_mean_better <- mean_better - min(y_range)

  }

  return(
    list(
      new_y_range = new_y_range,
      new_mean_worse = new_mean_worse,
      new_mean_better = new_mean_better
    )
  )
}

#' Linear Programming Bound on the Ideal Rule Benefit
#'
#' Provides tight bounds on the benefit of the ideal treatment rule over the best uniform
#'    treatment through linear programming. Also (optionally), provides an asymptotically
#'    conservative upper confidence bound at the specified confidence level (95% is the default).
#'
#' The estimate is the value at the solution
#'    to a linear programming problem over all possible joint distributions of Y^1,Y^0.
#'    The objective function is the treatment rule benefit, and the constraints are
#'    that all probabilities must add to 1 and be nonnegative, and that the joint distribution produces
#'    the mean and variance in each arm observed from the data.
#'
#' The upper confidence bound solves a linear program with the same objective, except that the contraints
#'     on the solution mean and variance are such that the solution must produce means and variances in each
#'     arm within the 95% confidence interval calculated from the data.
#'
#' @param scale Scale should be "higher" (default) if higher outcomes are beneficial,
#'     and "lower" otherwise
#' @param conf.int logical; whether a upper confidence interval should be provided
#' @param y_range the range of possible outcome values. either this or both m and M should be specified
#'    if m,M are specified this should be left as `NULL`
#' @param m,M the upper and lower bounds on the outcome range.
#'    Specifying these assumes a range of all integer values between and including m and M.
#'    Either m,M or y_range should be specified.
#'    If y_range is specified m,M should be left as `NULL`.
#' @inheritParams teh_closed_form_bounded
#'
#'
#' @returns
#' #' A list with two elements:
#'
#' * `estimate` the estimated benefit of the ideal treatment rule over the best uniform treatment
#' * `conf.int` The `level` confidence interval if `conf.int` is TRUE, `NULL` otherwise
#'
#' @examples
#'
#' #outcome is 1-5, variance 1 in each arm. treatment effect is 1
#' #   100 people in each arm
#' #   higher outcomes better
#' lp_ideal_rule_benefit(s2_1 = 1,
#'                       s2_0 = 1,
#'                       n_1 = 100,
#'                       n_0 = 100,
#'                       mean_1 = 3,
#'                       mean_0 = 2,
#'                       m = 1,
#'                       M = 5,
#'                       conf.int = TRUE)
#'
#' #outcome is 1-5, variance 1 in each arm. treatment effect is 1
#' #   100 people in each arm
#' #   lower outcomes better
#' lp_ideal_rule_benefit(s2_1 = 1,
#'                       s2_0 = 1,
#'                       n_1 = 100,
#'                       n_0 = 100,
#'                       mean_1 = 2,
#'                       mean_0 = 3,
#'                       m = 1,
#'                       M = 5,
#'                       scale = "lower",
#'                       conf.int = TRUE)
#'
#' #outcome range is (1,3,5,6,7,8,9,10)
#' lp_ideal_rule_benefit(s2_1 = 1,
#'                       s2_0 = 1,
#'                       n_1 = 100,
#'                       n_0 = 100,
#'                       mean_1 = 3,
#'                       mean_0 = 2,
#'                       y_range = c(1,3,5,6,7,8,9,10),
#'                       conf.int = TRUE)
#'
#' @export
lp_ideal_rule_benefit <- function(m = NULL,
                                  M = NULL,
                                  y_range = NULL,
                                  s2_1,
                                  s2_0,
                                  mean_1,
                                  mean_0,
                                  n_1 = NULL,
                                  n_0 = NULL,
                                  scale = "higher",
                                  conf.int = FALSE,
                                  level = 0.95)
{
  if ((is.null(m) | is.null(M)) & is.null(y_range)) {
    stop("Either the y_range or both m and M must be specified")
  }

  if ((!is.null(m) | !is.null(M)) & !is.null(y_range)) {
    stop("Only y_range or m,M can be specified")
  }

  if (!(scale %in% c("higher", "lower")))
    stop ("scale must be 'higher' or 'lower'")

  if (!(is.null(m) | is.null(M))) {
    if (m >= M)
      stop("m must be less than M")

    #if lower and upper bounds are specified, create range based on those
    y_range <- m:M

  }

  #assign which arms are worse or better based on scale
  if ((scale == "lower" &
       mean_1 <= mean_0) | (scale == "higher" & mean_1 >= mean_0)) {
    mean_better <- mean_1
    var_better <- s2_1
    n_better <- n_1

    mean_worse <- mean_0
    var_worse <- s2_0
    n_worse <- n_0

  } else {
    mean_better <- mean_0
    var_better <- s2_0
    n_better <- n_0

    mean_worse <- mean_1
    var_worse <- s2_1
    n_worse <- n_1

  }

  #transform the range to be nonnegative with higher outcomes beneficial
  temp <- transform_range (y_range,
                           mean_worse,
                           mean_better,
                           scale)

  new_y_range <- temp$new_y_range
  new_mean_worse <- temp$new_mean_worse
  new_mean_better <- temp$new_mean_better


  if (conf.int) {
    if (is.null(n_1) | is.null(n_0)) {
      stop("Values for n_1 and n_0 are required to form confidence intervals.")

    }

    #get confidence interval
    ci.int <- helper_lp_ci(
      y_range = new_y_range,
      new_mean_worse,
      new_mean_better,
      var_worse,
      var_better,
      n_worse,
      n_better,
      level = level
    )


  } else {
    ci.int <- NULL
  }

  #get bound estimate
  bound <- helper_lp_benefit(y_range = new_y_range,
                             new_mean_worse,
                             new_mean_better,
                             var_worse,
                             var_better)

  return(list(estimate = bound, conf.int = ci.int))


}

#' Helper Function to Calculate Linear Program Bound on Ideal Rule Benefit
#'
#'Calculates the estimated bound on the ideal treatment rule benefit as the value at the solution
#'    to a linear programming problem over all possible joint distributions of Y^1,Y^0.
#'    The objective function is the treatment rule benefit, and the constraints are
#'    that all probabilities must add to 1 and be nonnegative, and that the joint distribution produces
#'    the mean and variance in each arm observed from the data.
#'
#'
#' @param y_range The range of possible outcomes, should be nonnegative
#' @param mean_worse The estimated mean of the worse treatment
#' @param mean_better The estimated mean of the better treatment
#' @param var_worse The estimated variance of the worse treatment
#' @param var_better The estimated variance of the better treatment
#'
#' @return the estimated benefit of the ideal treatment rule over the best uniform treatment
#'
#' @examples #this is an internal function
helper_lp_benefit <- function(y_range,
                              mean_worse,
                              mean_better,
                              var_worse,
                              var_better) {
  #get joint distribution of outcome values in both arms
  support <- expand.grid("y_w" = y_range, "y_b" = y_range)

  #objective function - absolute value of treatment effect when less than 0
  objective.in <-
    (support$y_w - support$y_b) * (support$y_b < support$y_w)

  #constraint that all probs add to 1
  c1 <- rep(1, length(y_range) ^ 2)

  #left hand side of constraint on worse mean
  c2 <- support$y_w

  #left hand side of constraint on better mean
  c3 <- support$y_b

  #left hand side of constraint on worse second moment
  c4 <- support$y_w ^ 2

  #left hand side of constraint on better second moment
  c5 <- support$y_b ^ 2

  #matrix of left hand sides of constraints
  const.mat <- rbind(c1, c2, c3, c4, c5)

  #all constraints are equalities
  const.dir <- rep("==", 5)

  #right hand sides of constraints
  const.rhs <- c(1,
                 mean_worse,
                 mean_better,
                 var_worse + mean_worse ^ 2,
                 var_better + mean_better ^ 2)

  #solvea the linear programming problem
  #   returns the max value of the objective
  return(c(
    lpSolve::lp(
      direction = "min",
      objective.in = objective.in,
      const.mat = const.mat,
      const.dir = const.dir,
      const.rhs = const.rhs
    )$objval
    ,
    lpSolve::lp(
      direction = "max",
      objective.in = objective.in,
      const.mat = const.mat,
      const.dir = const.dir,
      const.rhs = const.rhs
    )$objval
  ))

}

#' Helper Function to Calculate Linear Program Upper Confidence Bound on Ideal Rule Benefit
#'
#'Calculates a (likely conservative) level% confidence bound on the ideal treatment rule benefit
#'    as the value of the solution to a linear programming problem over
#'    all possible joint distributions of Y^1,Y^0.
#'    The objective function is the treatment rule benefit, and the constraints are
#'    that all probabilities must add to 1 and be nonnegative, and that the joint distribution produces
#'    a mean and second absolute moment within the 95% confidence intervals calculated
#'    in each arm based on the data.
#'
#'
#' @param n_worse Sample size of worse treatment
#' @param n_better Sample size of better treatment
#' @param level confidence level between 0 and 1 for confidence interval;
#'     default is 0.95
#' @inheritParams helper_lp_benefit
#'
#' @return The `level` upper confidence bound
#' @export
#'
#' @examples #this is an internal function
helper_lp_ci <- function(y_range,
                         mean_worse,
                         mean_better,
                         var_worse,
                         var_better,
                         n_worse,
                         n_better,
                         level) {
  #get low and high quantiles for given confidence level
  q_low = 0.5 - level / 2
  q_high = 0.5 + level / 2

  #get joint distribution of outcome values in both arms
  support <- expand.grid("y_w" = y_range, "y_b" = y_range)

  #get bounds on range for use in se bound calculations
  range_low = min(y_range)
  range_high = max(y_range)

  #se for sample means of worse and better
  m_w_se <- sqrt(var_worse / n_worse)
  m_b_se <- sqrt(var_better / n_better)

  #se for sample var + sample mean squared of worse arm
  second_moment_w_se <- max((range_high - mean_worse) ^ 2,
                            (range_low - mean_worse) ^ 2) * var_worse -
    var_worse ^ 2 +
    4 * mean_worse * max(abs(range_high), abs(range_low)) * (var_worse +
                                                               mean_worse ^ 2) -
    8 * var_worse * mean_worse ^ 2 - 4 * mean_worse ^ 4

  second_moment_w_se <- sqrt(second_moment_w_se / n_worse)

  #se for sample var + sample mean squared of better arm
  second_moment_B_se <- max((range_high - mean_better) ^ 2,
                            (range_low - mean_better) ^ 2) * var_better -
    var_better ^ 2 +
    4 * mean_better * max(abs(range_high), abs(range_low)) * (var_better +
                                                                mean_better ^ 2) -
    8 * var_better * mean_better ^ 2 - 4 * mean_better ^ 4

  second_moment_B_se <- sqrt(second_moment_B_se / n_better)

  #objective function - absolute value of treatment effect when less than 0
  objective.in <-
    (support$y_w - support$y_b) * (support$y_b < support$y_w)

  #left hand side of constraint that all probs add to 1
  c1 <- rep(1, length(y_range) ^ 2)

  #left hand side of constraint on worse mean (used twice)
  c2 <- support$y_w

  #left hand side of constraint on better mean (used twice)
  c3 <- support$y_b

  #left hand side of constraint on worse second moment (used twice)
  c4 <- support$y_w ^ 2

  #left hand side of constraint on better second moment (used twice)
  c5 <- support$y_b ^ 2


  #left hand side of constraints
  const.mat <- rbind(c1, c2, c2, c3, c3, c4, c4, c5, c5)

  #c1 is equality, rest are within lower and upper bounds
  const.dir <- c("==", rep(c(">=", "<="), 4))

  #right hand side of constraints
  const.rhs <- c(
    1,
    stats::qnorm(c(q_low, q_high),
                 mean = mean_worse, sd = m_w_se),
    stats::qnorm(c(q_low, q_high),
                 mean = mean_better, sd = m_b_se),
    stats::qnorm(
      c(q_low, q_high),
      mean = var_worse + mean_worse ^ 2,
      sd = second_moment_w_se
    ),
    stats::qnorm(
      c(q_low, q_high),
      mean = var_better + mean_better ^ 2,
      sd = second_moment_B_se
    )
  )

  #solvea the linear programming problem
  #   returns the max value of the objective

  return(c(
    lpSolve::lp(
      direction = "min",
      objective.in = objective.in,
      const.mat = const.mat,
      const.dir = const.dir,
      const.rhs = const.rhs
    )$objval,
    lpSolve::lp(
      direction = "max",
      objective.in = objective.in,
      const.mat = const.mat,
      const.dir = const.dir,
      const.rhs = const.rhs
    )$objval
  ))

}

#' Linear programming rule benefit bounds with stratification variable
#'
#' @param m The lowest possible outcome value
#' @param M The highest possible outcome value
#' @param y_range The range of outcomes
#' @param s2_1 The vector of variances in arm 1 from each stratum
#' @param s2_0 The vector of variances in arm 0 from each stratum
#' @param mean_1 The vector of means in arm 1 from each stratum
#' @param mean_0 The vector of means in arm 0 from each stratum
#' @param strata_props The vector of proportions of each stratum in the population
#' @param scale Scale should be "higher" (default) if higher outcomes are beneficial,
#'     and "lower" otherwise
#'
#' @return A length 2 vector of the lower and upper bounds on treatment rule benefit
#' @export
#'
#' @examples #TBC
strata_lp_ideal_rule_benefit <- function(m = NULL,
                                         M = NULL,
                                         y_range = NULL,
                                         s2_1,
                                         s2_0,
                                         mean_1,
                                         mean_0,
                                         strata_props,
                                         scale = "higher"){

  if ((is.null(m) | is.null(M)) & is.null(y_range)) {
    stop("Either the y_range or both m and M must be specified")
  }

  if ((!is.null(m) | !is.null(M)) & !is.null(y_range)) {
    stop("Only y_range or m,M can be specified")
  }

  if (!(scale %in% c("higher", "lower")))
    stop ("scale must be 'higher' or 'lower'")

  if(!(length(s2_1)==length(s2_0) &
       length(s2_1)==length(mean_1) &
       length(s2_1)==length(mean_0))){
    stop("Length of all mean, variance, and strata proportions vector must be the same.")
  }

  if(sum(strata_props) != 1 | sum(strata_props < 0) != 0 | sum(strata_props > 1) != 0){
    stop ("Strata proportions must be between 0 and 1 and sum to 1.")
  }


  if (!(is.null(m) | is.null(M))) {
    if (m >= M)
      stop("m must be less than M")

    #if lower and upper bounds are specified, create range based on those
    y_range <- m:M

  }

  #get the number of strata
  n_strata <- length(mean_1)

  mean_worse <- rep(NA,n_strata)
  mean_better <- mean_worse
  var_worse <- mean_worse
  var_better <- mean_better


  if ((scale == "lower" &
       strata_props %*% mean_1 <= strata_props %*% mean_0) |
      (scale == "higher" & strata_props %*% mean_1 >= strata_props %*% mean_0)) {
    overall_mean_better <- strata_props %*% mean_1
    overall_mean_worse <- strata_props %*% mean_0

  } else {
    overall_mean_better <- strata_props %*% mean_0
    overall_mean_worse <- strata_props %*% mean_1

  }

  for(i in 1:n_strata){

    #assign which arms are worse or better based on scale
    if ((scale == "lower" &
         mean_1[i] <= mean_0[i]) | (scale == "higher" & mean_1[i] >= mean_0[i])) {
      mean_better[i] <- mean_1[i]
      var_better[i] <- s2_1[i]

      mean_worse[i] <- mean_0[i]
      var_worse[i] <- s2_0[i]

    } else {
      mean_better[i] <- mean_0[i]
      var_better[i] <- s2_0[i]

      mean_worse[i] <- mean_1[i]
      var_worse[i] <- s2_1[i]

    }
  }


  temp <- transform_range (y_range,
                           mean_worse,
                           mean_better,
                           scale)

  new_y_range <- temp$new_y_range
  new_mean_worse <- temp$new_mean_worse
  new_mean_better <- temp$new_mean_better

  #get bounds on further benefit beyond strata based rule

  further_benefit <- matrix(NA,nrow = n_strata,ncol = 2)

  for(i in 1:n_strata){

    further_benefit[i,] <- helper_lp_benefit(y_range = new_y_range,
                                             new_mean_worse[i],
                                             new_mean_better[i],
                                             var_worse[i],
                                             var_better[i])

  }


  if (scale == "lower") {
    uniform_value <- max(y_range) - overall_mean_better

  } else{
    uniform_value  <- overall_mean_better - min(y_range)

  }

  browser()


  strata_rule_benefit <- round(strata_props %*% new_mean_better -
                                             uniform_value,10)

  further_benefit_total <- as.vector(strata_props %*% further_benefit)

  return(c(strata_rule_benefit) + further_benefit_total)

}
