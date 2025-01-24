#' HulC-Based Confidence Interval for Adaptive RMST
#'
#'
#' @param survival_data A data frame containing survival times and event indicators.
#'   Must have:
#'   \itemize{
#'     \item \strong{X}: Observed time (event or censoring). Default time unit: years.
#'     \item \strong{delta}: Event indicator (1 = event, 0 = censored).
#'     \item \strong{A}: Treatment arm indicator (0 = control, 1 = treatment).
#'   }
#' @param alpha Numeric value for the nominal significance level (default is \code{0.05}).
#'   Used to decide \code{B} and report (1-alpha) confidence intervals.
#' @param type A character string indicating the rule for choosing \code{B}.
#'   Either \code{"anti-conservative"} or \code{"conservative"}.
#' @param min_time_fix Lower bound of the time range for RMST estimation
#'   (default is \code{0.2}).
#' @param max_time_fix Upper bound of the time range for RMST estimation
#'   (default is \code{1}).
#' @param tol A numeric tolerance for optimization (default is \code{1e-7}).
#' @param maxeval The maximum number of function evaluations for optimization
#'   (default is \code{70}).
#' @param verbose Logical indicating whether to print verbose messages (default is \code{FALSE}).
#'
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\strong{K_star}}{The overall RMST difference estimate from the full dataset.}
#'   \item{\strong{L_star}}{The overall optimal restriction time estimate from the full dataset.}
#'   \item{\strong{ci.K}}{Lower and upper bounds of the (1-alpha) confidence interval for \code{K_star}.}
#' }
#'
#' @seealso
#'   \code{\link{AdaRMST.ct}} for point estimation of the adaptive RMST,
#'   \code{\link{choose_B_floor}}, \code{\link{choose_B_ceiling}} for determining \code{B}.
#' @references
#' This function implements the method described in:
#'
#' Jinghao Sun, Douglas E. Schaubel, Eric J. Tchetgen Tchetgen (2025).
#' \emph{Beyond Fixed Restriction Time: Adaptive Restricted Mean Survival Time Methods in Clinical Trials}.
#' @export
AdaRMST.hulc <- function(
    survival_data,
    alpha = 0.05,
    type = "anti-conservative",
    min_time_fix = 0.2,
    max_time_fix = 1,
    tol = 1e-7,
    maxeval = 70,
    verbose = FALSE
) {
  # 1. Obtain a point estimate (L_star, K_star) from the entire dataset
  #    using c = 0 in AdaRMST.ct
  res.point <- AdaRMST.ct(
    survival_data = survival_data,
    min_time_fix = min_time_fix,
    max_time_fix = max_time_fix,
    L.tilde = (min_time_fix + max_time_fix) / 2,
    c = 0,
    tol = tol,
    maxeval = maxeval,
    verbose = verbose
  )

  # 2. Determine how many subsets B using external functions
  if (type == "anti-conservative") {
    B <- choose_B_floor(alpha)
  } else if (type == "conservative") {
    B <- choose_B_ceiling(alpha)
  } else {
    stop("Invalid 'type' argument. Must be 'anti-conservative' or 'conservative'.")
  }

  # 3. Randomly permute the row indices, then split into B groups
  n <- nrow(survival_data)
  group_assignments <- cut(
    sample(seq_len(n)),
    breaks = B,
    labels = FALSE
  )
  split_data <- split(survival_data, group_assignments)

  # 4. For each subset, re-estimate (L_star, K_star)
  L_vec <- numeric(B)
  K_vec <- numeric(B)
  for (i in seq_len(B)) {
    sub_data <- split_data[[i]]
    res.tmp <- AdaRMST.ct(
      survival_data = sub_data,
      min_time_fix = min_time_fix,
      max_time_fix = max_time_fix,
      L.tilde = (min_time_fix + max_time_fix) / 2,
      c = 0,
      tol = tol,
      maxeval = maxeval,
      verbose = verbose
    )
    L_vec[i] <- res.tmp$L_star
    K_vec[i] <- res.tmp$K_star
  }

  # 5. Construct HulC intervals from min and max of each subset estimate
  ci.L <- c(min(L_vec), max(L_vec))
  ci.K <- c(min(K_vec), max(K_vec))

  # 6. Return final results
  res <- list(
    K_star = res.point$K_star,
    L_star = res.point$L_star,
    ci.K   = ci.K
    # ci.L   = ci.L  # Uncomment if you want to return L-interval as well
  )
  return(res)
}

