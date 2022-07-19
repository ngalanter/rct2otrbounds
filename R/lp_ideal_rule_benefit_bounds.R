#' Helper Function to Transform Outcome
#'
#' @param y_range the range of outcomes
#' @param mean_worse mean of the worse treatment arm
#' @param mean_better mean of the better treatment arm
#' @param scale whether "higher" or "lower" outcomes are better
#'
#' @return a list of the new range, new mean of the better arm, and new mean of the worse arm
#'
#' @examples
#' transform_range(1:10,4,6,scale = "higher")
#' transform_range(-1*(1:10),-4,-6,scale = "lower")
transform_range <- function(y_range,
                            mean_worse,
                            mean_better,
                            scale){

  #flip outcome scale if lower outcomes are better, transform range to start at 0
  if(scale == "lower"){

    new_y_range <- -1*y_range + max(y_range)
    new_mean_worse <- max(y_range) - mean_worse
    new_mean_better <- max(y_range) - mean_better

  } else{

    new_y_range <- y_range - min(y_range)
    new_mean_worse <- mean_worse - min(y_range)
    new_mean_better <- mean_better - min(y_range)

  }

  return(list(new_y_range = new_y_range,
              new_mean_worse = new_mean_worse,
              new_mean_better = new_mean_better))
}

#' Linear Programming Bound on the Ideal Rule Benefit
#'
#' @param scale Scale should be "higher" (default) if higher outcomes are beneficial,
#'     and "lower" otherwise
#' @param conf.upper logical; whether a upper confidence bound should be provided
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
#' * `conf.upper` The `level` upper confidence bound if `conf.upper` is TRUE, `NULL` otherwise
#'
#' @examples #tbc
#'
#' @export
lp_ideal_rule_benefit <- function(m = NULL,
                                          M = NULL,
                                          y_range = NULL,
                                              s2_1,
                                              s2_0,
                                              mean_1,
                                              mean_0,
                                              n_1=NULL,
                                              n_0=NULL,
                                              scale = "higher",
                                              conf.upper = FALSE,
                                              level = 0.95)
{

  if((is.null(m) | is.null(M)) & is.null(y_range)){
    stop("Either the y_range or both m and M must be specified")
  }

  if((!is.null(m) | !is.null(M)) & !is.null(y_range)){
    stop("Only y_range or m,M can be specified")
  }

  if(!(scale %in% c("higher","lower"))) stop ("scale must be 'higher' or 'lower'")

  if(!(is.null(m) | is.null(M))){

    if(m >= M) stop("m must be less than M")

    #if lower and upper bounds are specified, create range based on those
    y_range <- m:M

  }

  #assign which arms are worse or better based on scale
  if((scale == "lower" & mean_1 <= mean_0) | (scale == "higher" & mean_1 >= mean_0)){

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


  if(conf.upper){

    if(is.null(n_1) | is.null(n_0)){

      stop("Values for n_1 and n_0 are required to form confidence intervals.")

    }

    #get upper confidence bound
    ci.upper <- helper_lp_ci(y_range = new_y_range,
                               new_mean_worse,
                               new_mean_better,
                               var_worse,
                               var_better,
                               n_worse,
                               n_better,
                               level = level)


  } else {ci.upper <- NULL}

  #get bound estimate
  bound <- helper_lp_benefit(y_range = new_y_range,
                             new_mean_worse,
                             new_mean_better,
                             var_worse,
                             var_better)

  return(list(estimate = bound, conf.upper = ci.upper))


}

#' Helper Function to Calculate Linear Program Bound on Ideal Rule Benefit
#'
#' @param y_range The range of possible outcomes, should be nonnegative
#' @param mean_worse The estimated mean of the worse treatment
#' @param mean_better The estimated mean of the better treatment
#' @param var_worse The estimated variance of the worse treatment
#' @param var_better The estimated variance of the better treatment
#'
#' @return the estimated benefit of the ideal treatment rule over the best uniform treatment
#'
#' @examples #tbc
helper_lp_benefit <- function(y_range,
                               mean_worse,
                               mean_better,
                               var_worse,
                               var_better){


  support <- expand.grid("y_w" = y_range,"y_b" = y_range)

  #absolute value of treatment effect when less than 0
  objective.in <- (support$y_w - support$y_b) * (support$y_b < support$y_w)

  #constraint that all probs add to 1
  c1 <- rep(1,length(y_range)^2)

  #constraint on worse mean
  c2 <- support$y_w

  #contsraint on better mean
  c3 <- support$y_b

  #constraint on worse second moment
  c4 <- support$y_w^2

  #constraint on better second moment
  c5 <- support$y_b^2

  const.mat <- rbind(c1,c2,c3,c4,c5)


  const.dir <- rep("==",5)

  const.rhs <- c(1,mean_worse,mean_better,
                 var_worse+mean_worse^2,var_better+mean_better^2)

  return(
    lpSolve::lp(direction = "max",objective.in = objective.in,
       const.mat = const.mat, const.dir = const.dir,
       const.rhs = const.rhs)$objval
  )

}

#' Helper Function to Calculate Linear Program Upper Confidence Bound on Ideal Rule Benefit
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
#' @examples #tbc
helper_lp_ci <- function(y_range,
                                        mean_worse,
                                        mean_better,
                                        var_worse,
                                        var_better,
                                        n_worse,
                                        n_better,
                         level){

  q_low = 0.5 - level/2
  q_high = 0.5 + level/2

  support <- expand.grid("y_w" = y_range,"y_b" = y_range)

  range_low = min(y_range)

  range_high = max(y_range)

  #se for sample means of worse and better
  m_w_se <- sqrt(var_worse/n_worse)
  m_b_se <- sqrt(var_better/n_better)

  #se for sample var + sample mean squared
  second_moment_w_se <- max((range_high-mean_worse)^2,
                            (range_low-mean_worse)^2)*var_worse-var_worse^2 +
    4*mean_worse*max(abs(range_high),abs(range_low))*(var_worse+mean_worse^2)-
    8*var_worse*mean_worse^2 - 4*mean_worse^4

  second_moment_w_se <- sqrt(second_moment_w_se/n_worse)

  second_moment_B_se <- max((range_high-mean_better)^2,
                            (range_low-mean_better)^2)*var_better-var_better^2 +
    4*mean_better*max(abs(range_high),abs(range_low))*(var_better+mean_better^2)-
    8*var_better*mean_better^2 - 4*mean_better^4

  second_moment_B_se <- sqrt(second_moment_B_se/n_better)

  #absolute value of treatment effect when less than 0
  objective.in <- (support$y_w - support$y_b) * (support$y_b < support$y_w)

  #constraint that all probs add to 1
  c1 <- rep(1,length(y_range)^2)

  #constraint on worse mean (used twice)
  c2 <- support$y_w

  #contsraint on better mean (used twice)
  c3 <- support$y_b

  #constraint on worse second moment (used twice)
  c4 <- support$y_w^2

  #constraint on better second moment (used twice)
  c5 <- support$y_b^2


  const.mat <- rbind(c1,c2,c2,c3,c3,c4,c4,c5,c5)

  #c1 is equality, rest are within lower and upper bounds
  const.dir <- c("==",rep(c(">=","<="),4))

  const.rhs <- c(1,
                 stats::qnorm(c(q_low,q_high),
                       mean = mean_worse, sd = m_w_se),
                 stats::qnorm(c(q_low,q_high),
                       mean = mean_better, sd = m_b_se),
                 stats::qnorm(c(q_low,q_high),
                       mean =var_worse+mean_worse^2, sd = second_moment_w_se),
                 stats::qnorm(c(q_low,q_high),
                       mean = var_better+mean_better^2, sd = second_moment_B_se))


  lpSolve::lp(direction = "max",objective.in = objective.in,
     const.mat = const.mat, const.dir = const.dir,
     const.rhs = const.rhs)

}





