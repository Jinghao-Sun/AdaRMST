% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AdaRMST.ct.R
\name{AdaRMST.ct}
\alias{AdaRMST.ct}
\title{Continuous-Time Adaptive Restricted Mean Survival Time Method: Point Estimate}
\usage{
AdaRMST.ct(
  survival_data,
  min_time_fix = 0.2,
  max_time_fix = 1,
  L.tilde = (min_time_fix + max_time_fix)/2,
  c = 0.032/(max_time_fix - min_time_fix)^2,
  tol = 1e-07,
  maxeval = 70,
  verbose = FALSE
)
}
\arguments{
\item{survival_data}{A data frame containing survival information. It should have three columns:
\itemize{
\item \strong{X}: Observed event or censoring time. Default time unit: year.
\item \strong{A}: A binary value indicating treatment arm (0 = control, 1 = treatment).
\item \strong{delta}: An event indicator where 1 = event (not censored), 0 = censored.
}}

\item{min_time_fix}{The lower bound of the time range over which to search for the RMST (default is \code{0.2}).}

\item{max_time_fix}{The upper bound of the time range over which to search for the RMST (default is \code{1}).}

\item{L.tilde}{The initial guess (starting value) for the optimal restriction time. By default, it is
set to \code{(min_time_fix + max_time_fix)/2}.}

\item{c}{A numeric penalty coefficient. A larger value imposes a bigger penalty on deviating from
\code{L.tilde}. By default, it is set to \code{0.032 / (max_time_fix - min_time_fix)^2}
(intended for time units in years). If you set \code{c = 0}, there is no penalty.}

\item{tol}{The tolerance level for the optimization routine (default is \code{1e-7}).}

\item{maxeval}{The maximum number of evaluations allowed during optimization (default is \code{70}).}

\item{verbose}{A logical value indicating whether to print optimization progress messages (default is \code{FALSE}).}
}
\value{
A list with the following elements:
\itemize{
\item \strong{K_star}: The estimated treatment effect (difference in RMST between arms) at \code{L_star}.
\item \strong{L_star}: The estimated optimal restriction time.
}
}
\description{
Calculates the treatment effect point estimate using an adaptive Restricted Mean Survival Time (RMST) approach
over a continuous-time restriction candidate set. By default, it searches over the range
between \code{min_time_fix} and \code{max_time_fix} to find an optimal restriction time that maximizes
a penalized objective function.
}
\references{
This function implements the method described in:

Jinghao Sun, Douglas E. Schaubel, Eric J. Tchetgen Tchetgen (2025).
\emph{Beyond Fixed Restriction Time: Adaptive Restricted Mean Survival Time Methods in Clinical Trials}.
}
