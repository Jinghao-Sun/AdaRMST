#' Continuous-Time Adaptive Restricted Mean Survival Time Method: Point Estimate
#'
#' Calculates the treatment effect point estimate using an adaptive Restricted Mean Survival Time (RMST) approach
#' over a continuous-time restriction candidate set. By default, it searches over the range
#' between \code{min_time_fix} and \code{max_time_fix} to find an optimal restriction time that maximizes
#' a penalized objective function.
#'
#' @param survival_data A data frame containing survival information. It should have three columns:
#'   \itemize{
#'     \item \strong{X}: Observed event or censoring time. Default time unit: year.
#'     \item \strong{A}: A binary value indicating treatment arm (0 = control, 1 = treatment).
#'     \item \strong{delta}: An event indicator where 1 = event (not censored), 0 = censored.
#'   }
#'
#' @param min_time_fix The lower bound of the time range over which to search for the RMST (default is \code{0.2}).
#' @param max_time_fix The upper bound of the time range over which to search for the RMST (default is \code{1}).
#' @param L.tilde The initial guess (starting value) for the optimal restriction time. By default, it is
#'   set to \code{(min_time_fix + max_time_fix)/2}.
#' @param c A numeric penalty coefficient. A larger value imposes a bigger penalty on deviating from
#'   \code{L.tilde}. By default, it is set to \code{0.032 / (max_time_fix - min_time_fix)^2}
#'   (intended for time units in years). If you set \code{c = 0}, there is no penalty.
#' @param tol The tolerance level for the optimization routine (default is \code{1e-7}).
#' @param maxeval The maximum number of evaluations allowed during optimization (default is \code{70}).
#' @param verbose A logical value indicating whether to print optimization progress messages (default is \code{FALSE}).
#'
#' @importFrom nloptr nloptr
#' @importFrom survRM2 rmst2
#'
#' @return A list with the following elements:
#'   \itemize{
#'     \item \strong{K_star}: The estimated treatment effect (difference in RMST between arms) at \code{L_star}.
#'     \item \strong{L_star}: The estimated optimal restriction time.
#'   }
#'
#' @references
#' This function implements the method described in:
#'
#' Jinghao Sun, Douglas E. Schaubel, Eric J. Tchetgen Tchetgen (2025).
#' \emph{Beyond Fixed Restriction Time: Adaptive Restricted Mean Survival Time Methods in Clinical Trials}.
#'
#' @export
AdaRMST.ct <- function(
    survival_data,
    min_time_fix = 0.2,
    max_time_fix = 1,
    L.tilde = (min_time_fix + max_time_fix) / 2,
    c = 0.032 / (max_time_fix - min_time_fix)^2,
    tol = 1e-7,
    maxeval = 70,
    verbose = FALSE
) {
  # Function body
  n = nrow(survival_data)
  time = survival_data$X
  status = survival_data$delta
  arm = survival_data$A
  max_followup_time = min(
    max(survival_data$X[survival_data$A == 1]),
    max(survival_data$X[survival_data$A == 0])
  )
  if (max_followup_time < min_time_fix) {
    return(list(L_star = NA, K_star = NA))
  }

  min_time = min_time_fix
  max_time = max_time_fix
  lower_bound = min_time
  upper_bound = max_time

  # Vectorized function to compute RMST difference
  RMST.diff.sample = Vectorize(function(t){
    rmst.tmp = rmst2(time, status, arm, t)
    rmstdiff = rmst.tmp$unadjusted.result[1,1]
    return(rmstdiff)
  })

  # Objective function for optimization
  Mn_get_safe <- function(tau) {
    tryCatch({
      rmst.tmp <- rmst2(time, status, arm, tau)
      Mn <- (rmst.tmp$unadjusted.result[1,1])^2 /
        (rmst.tmp$RMST.arm1$rmst.var + rmst.tmp$RMST.arm0$rmst.var) / n -
        c * (tau - L.tilde)^2
      return(Mn)
    }, error = function(e) {
      return(-999)
    }, warning = function(w) {
      return(-999)
    })
  }

  # Use NLOPT for global+local optimization
  algorithm = "NLOPT_GN_DIRECT"
  result = nloptr(
    x0 = L.tilde,
    eval_f = function(x) -Mn_get_safe(x),  # we minimize the negative
    lb = lower_bound,
    ub = upper_bound,
    opts = list(
      algorithm = algorithm,
      xtol_rel = tol,
      xtol_abs = tol,
      ftol_abs = tol,
      ftol_rel = tol,
      maxeval = maxeval,
      maxtime = 1,
      local_optimizer = "NLOPT_LD_LBFGS",
      print_level = ifelse(verbose, 1, 0)
    )
  )

  L_star = result$solution
  Mn_star = -result$objective
  K_star = RMST.diff.sample(L_star)
  res = list(K_star = K_star, L_star = L_star)
  return(res)
}


#' Bootstrap-Based Confidence Intervals for Adaptive RMST in Continuous Time
#'
#' Applies a bootstrap procedure to estimate the confidence intervals of the adaptive Restricted Mean Survival Time (RMST) point estimates
#' obtained from \code{\link{AdaRMST.ct}}. The function resamples the input \code{survival_data} a specified
#' number of times (\code{B.boot}), recalculates the estimated optimal restriction time (\code{L_star}) and
#' the corresponding treatment effect (\code{K_star}), then uses percentile-based intervals from
#' \code{\link[boot]{boot.ci}} to produce confidence intervals.
#'
#' @param survival_data A data frame containing survival information. Must have the same structure
#'   required by \code{\link{AdaRMST.ct}}, typically including:
#'   \itemize{
#'     \item \strong{X}: Observed event or censoring time (default time unit: years).
#'     \item \strong{A}: A binary indicator for treatment arm (0 = control, 1 = treatment).
#'     \item \strong{delta}: An event indicator where 1 = event (not censored) and 0 = censored.
#'   }
#' @param alpha A numeric value for the nominal significance level (default is \code{0.05}). \code{boot.ci} calculates
#'   (1 - alpha) confidence intervals.
#' @param B.boot An integer specifying the number of bootstrap replicates (default is \code{1000}).
#' @param min_time_fix The lower bound of the time range used by \code{\link{AdaRMST.ct}}
#'   (default is \code{0.2}).
#' @param max_time_fix The upper bound of the time range used by \code{\link{AdaRMST.ct}}
#'   (default is \code{1}).
#' @param L.tilde A numeric value representing the initial guess of the optimal restriction time
#'   (default is \code{(min_time_fix + max_time_fix)/2}).
#' @param c A numeric penalty coefficient (default is \code{0.032 / (max_time_fix - min_time_fix)^2} when
#'   the time unit is years). A larger value imposes a bigger penalty for deviating from \code{L.tilde}.
#' @param tol A numeric value specifying the tolerance for the optimization method
#'   (default is \code{1e-7}).
#' @param maxeval The maximum number of evaluations for the optimization routine (default is \code{70}).
#' @param verbose A logical indicating whether to print progress messages
#'   (default is \code{FALSE}).
#'
#' @return A named list with the following components:
#'   \itemize{
#'     \item \strong{K_star}: The treatment effect (difference in RMST) estimated by \code{\link{AdaRMST.ct}}
#'       on the original data.
#'     \item \strong{L_star}: The estimated optimal restriction time (on the original data).
#'     \item{\strong{ci.K}}: {A numeric vector of length 2 giving the lower and upper bounds of the (1-alpha) bootstrap confidence interval for \code{K_star}.}
#'     \item{\strong{ci.L}}: {A numeric vector of length 2 giving the lower and upper bounds of the (1-alpha) bootstrap confidence interval for \code{L_star}.}
#'   }
#'
#' @references
#' This function implements the method described in:
#'
#' Jinghao Sun, Douglas E. Schaubel, Eric J. Tchetgen Tchetgen (2025).
#' \emph{Beyond Fixed Restriction Time: Adaptive Restricted Mean Survival Time Methods in Clinical Trials}.
#'
#' @seealso
#'   \code{\link{AdaRMST.ct}},
#'   \code{\link[boot]{boot}},
#'   \code{\link[boot]{boot.ci}}
#'
#' @importFrom boot boot boot.ci
#' @export
AdaRMST.ct.ci <- function(
    survival_data,
    alpha = 0.05,
    B.boot = 1000,
    min_time_fix = 0.2,
    max_time_fix = 1,
    L.tilde = (min_time_fix + max_time_fix)/2,
    c = 0.032 / (max_time_fix - min_time_fix)^2,
    tol = 1e-7,
    maxeval = 70,
    verbose = FALSE
) {
  # Internal statistic function for boot
  boot.est <- function(data, indices) {
    # Resample data using the indices given by the bootstrap
    survival_data_boot <- data[indices, ]
    # Call AdaRMST.ct to estimate L_star and K_star on the bootstrap sample
    res.tmp <- AdaRMST.ct(
      survival_data_boot,
      min_time_fix = min_time_fix,
      max_time_fix = max_time_fix,
      L.tilde = L.tilde,
      c = c,
      tol = tol,
      maxeval = maxeval,
      verbose = verbose
    )
    c(L.hat = res.tmp$L_star, K.hat = res.tmp$K_star)
  }

  # Perform bootstrap resampling
  bootstrap_results <- boot(
    data = survival_data,
    statistic = boot.est,
    R = B.boot,
    stype = "i"
  )

  # Compute percentile-based confidence intervals
  ci.L <- boot.ci(bootstrap_results, conf = 1- alpha, type = "perc", index = 1)
  ci.K <- boot.ci(bootstrap_results, conf = 1- alpha, type = "perc", index = 2)

  # Organize results
  res.ls <- list(
    K_star = bootstrap_results$t0[2],
    L_star = bootstrap_results$t0[1],
    ci.K = ci.K$percent[4:5],   # lower and upper percentile CI for K_star
    ci.L = ci.L$percent[4:5]  # lower and upper percentile CI for L_star
  )
  return(res.ls)
}



