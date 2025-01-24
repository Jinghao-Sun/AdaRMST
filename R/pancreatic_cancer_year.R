#' Reconstructed Pancreatic Cancer Trial Data in the Transformed Time Unit (years)
#'
#'
#' @format A data frame with 866 rows and 3 columns:
#' \describe{
#'   \item{\code{X}}{Observed event or censoring time (time unit: years)}
#'   \item{\code{A}}{A binary indicator for treatment arm (0 = control, 1 = treatment)}
#'   \item{\code{delta}}{An event indicator where 1 = event (not censored) and 0 = censored.}
#' }
#'
#' @source
#' Jinghao Sun, Douglas E. Schaubel, Eric J. Tchetgen Tchetgen (2025).
#' \emph{Beyond Fixed Restriction Time: Adaptive Restricted Mean Survival Time Methods in Clinical Trials}.
#'
#'
"pancreatic_cancer_year"
