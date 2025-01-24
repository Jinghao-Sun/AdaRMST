#' Compute the Number of Folds for HulC: Conservative Approach
#'
#' Calculates the number of split folds (\code{B}) required for HulC
#' method using a conservative rule based on the specified significance level (\code{alpha}).
#'
#' The conservative version ensures that the number of folds is rounded up to the nearest integer,
#' providing a higher number of splits which may lead to more robust confidence interval estimates.
#'
#' @param alpha Numeric. Significance level for determining the number of folds.
#'   Common choices are \code{0.05} for a 95\% confidence level or \code{0.01} for a 99\% confidence level.
#'   The default is \code{0.05}.
#'
#' @return Integer. The calculated number of split folds (\code{B}) for the HULC method.
#'
#' @export
choose_B_ceiling <- function(alpha = 0.05) {
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value between 0 and 1.")
  }
  return(ceiling(1 - log2(alpha)))
}


#' Compute the Number of Folds for HulC: Anti-Conservative Approach
#'
#' Calculates the number of split folds (\code{B}) required for HulC
#' method using an anti-conservative rule based on the specified significance level (\code{alpha}).
#'
#' The anti-conservative version ensures that the number of folds is rounded down to the nearest integer,
#' providing a lower number of splits which may increase power.
#'
#' @param alpha Numeric. Significance level for determining the number of folds.
#'   Common choices are \code{0.05} for a 95\% confidence level or \code{0.01} for a 99\% confidence level.
#'   The default is \code{0.05}.
#'
#' @return Integer. The calculated number of split folds (\code{B}) for the HULC method.
#'
#'
#' @export
choose_B_floor <- function(alpha = 0.05) {
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value between 0 and 1.")
  }
  return(floor(1 - log2(alpha)))
}
