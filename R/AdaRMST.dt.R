#' Discrete (Grid) Search for Adaptive Restricted Mean Survival Time
#'
#' Performs a discrete (grid-based) search for the optimal restriction time (\code{L_star})
#' in calculating the Restricted Mean Survival Time (RMST) difference between two arms.
#' The function evaluates a penalized objective function at each time in \code{L.grid}
#' (filtered by \code{max_followup_time}), selecting the time that maximizes this criterion
#' as \code{L_star}. It returns the corresponding treatment effect (\code{K_star}),
#'confidence interval, and a \emph{p}-value.
#'
#' @param survival_data A data frame containing survival information with columns:
#'   \itemize{
#'     \item \strong{X}: Observed event/censoring time. Default time unit: year.
#'     \item \strong{delta}: Event indicator (1 = event, 0 = censored).
#'     \item \strong{A}: A binary indicator of treatment arm (0 = control, 1 = treatment).
#'   }
#' @param alpha A numeric value for the nominal significance level (default is \code{0.05}). A (1 - alpha) asymptotic confidence interval will be returned.
#' @param min_time_fix The lower bound of the time range used in the grid search
#'   (default is \code{0.2}).
#' @param max_time_fix The upper bound of the time range used in the grid search
#'   (default is \code{1}).
#' @param c A numeric penalty coefficient applied to deviations from \code{L.tilde}
#'   (default is \code{0.08 / (max_time_fix - min_time_fix)^2} with year as unit).
#' @param L.grid A numeric vector of candidate restriction times (default is
#'   \code{seq(min_time_fix, max_time_fix, length.out = floor(1.5* nrow(survival_data)^0.25))} according to infill aymptotic regime results).
#' @param L.tilde.index An integer index indicating which element of \code{L.grid} is
#'   considered \code{L.tilde} (default is the midpoint:
#'   \code{floor((length(L.grid)+1)/2)}).
#'
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{\strong{K_star}}{The estimated treatment effect (difference in RMST) at \code{L_star}.}
#'   \item{\strong{L_star}}{The selected optimal restriction time from the grid.}
#'   \item{\strong{ci.K}}{Lower and upper bounds of the (1-alpha) confidence interval for \code{K_star}.}
#'   \item{\strong{p}}{A \emph{p}-value associated with the treatment effect measured by RMST difference.}
#' }
#'
#' @references
#' This function implements the method described in:
#'
#' Jinghao Sun, Douglas E. Schaubel, Eric J. Tchetgen Tchetgen (2025).
#' \emph{Beyond Fixed Restriction Time: Adaptive Restricted Mean Survival Time Methods in Clinical Trials}.
#'
#' @importFrom survRM2 rmst2
#' @export
AdaRMST.dt <- function(
    survival_data,
    alpha = 0.05,
    min_time_fix = 0.2,
    max_time_fix = 1,
    c = 0.08 / (max_time_fix - min_time_fix)^2,
    L.grid = seq(min_time_fix, max_time_fix, length.out = floor(1.5* nrow(survival_data)^0.25)),
    L.tilde.index = floor((length(L.grid) + 1) / 2)
) {
  n <- nrow(survival_data)
  time <- survival_data$X
  status <- survival_data$delta
  arm <- survival_data$A

  # Determine the maximum follow-up time common to both arms
  max_followup_time <- min(
    max(survival_data$X[survival_data$A == 1]),
    max(survival_data$X[survival_data$A == 0])
  )
  if (max_followup_time < min_time_fix) {
    return(list(
      K_star = NA,
      L_star = NA,
      ci.K = NA,
      p = NA
    ))
  }

  # Restrict the search space if the max follow-up time is < max_time_fix
  max_time_rmst <- min(max_followup_time, max_time_fix)
  L.tilde <- L.grid[L.tilde.index]

  # Valid times in L.grid cannot exceed max_time_rmst
  L.grid.valid <- L.grid[L.grid <= max_time_rmst]

  # Helper function to compute RMST difference
  RMST.diff.sample <- Vectorize(function(t) {
    rmst.tmp <- rmst2(time, status, arm, t)
    rmstdiff <- rmst.tmp$unadjusted.result[1, 1]
    return(rmstdiff)
  })

  # Objective function (penalized) at each candidate L
  Mn_get_safe <- function(tau) {
    tryCatch({
      rmst.tmp <- rmst2(time, status, arm, tau)
      # Square of the difference in RMST, scaled by variance
      Mn <- (rmst.tmp$unadjusted.result[1,1])^2 /
        (rmst.tmp$RMST.arm1$rmst.var + rmst.tmp$RMST.arm0$rmst.var) / n -
        c * (tau - L.tilde)^2
      Mn
    }, error = function(e) {
      -999
    }, warning = function(w) {
      -999
    })
  }

  # Grid search over valid L values
  algorithm <- "grid search"
  Mn_vals <- vapply(L.grid.valid, Mn_get_safe, numeric(1))

  # Identify L_star that maximizes the objective
  L_star <- L.grid.valid[which.max(Mn_vals)]
  Mn_star <- max(Mn_vals)

  # Recompute RMST difference at L_star
  rmst.tmp <- rmst2(time, status, arm, L_star, alpha = alpha)
  K_star <- rmst.tmp$unadjusted.result[1,1]

  # Assemble results
  res <- list(
    K_star = K_star,
    L_star = L_star,
    ci.K = rmst.tmp$unadjusted.result[1,2:3],
    p = rmst.tmp$unadjusted.result[1,4]
  )
  return(res)
}

