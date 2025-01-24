
# AdaRMST

<!-- badges: start -->
<!-- badges: end -->

The goal of **AdaRMST** is to implement a set of novel adaptive restricted mean survival time methods for clinical trials, especially under non-proportional hazards. See the methodology paper: *Jinghao Sun, Douglas E. Schaubel, Eric J. Tchetgen Tchetgen (2025). Beyond Fixed Restriction Time: Adaptive Restricted Mean Survival Time Methods in Clinical Trials.*

Restricted mean survival time (RMST) offers a compelling nonparametric alternative to hazard ratios for right-censored time-to-event data, particularly when the proportional hazards assumption is violated. By capturing the total event-free time over a specified horizon, RMST provides an intuitive and clinically meaningful measure of absolute treatment benefit. Nonetheless, selecting the restriction time $L$ poses challenges: choosing a small $L$ can overlook late-emerging benefits, whereas a large $L$, often underestimated in its impact, may inflate variance and undermine power. We propose a novel data-driven, adaptive procedure that identifies the optimal restriction time $L^*$ from a continuous range by maximizing a criterion balancing effect size and estimation precision. Consequently, our procedure is particularly useful when the pattern of the treatment effect is unknown at the design stage. We provide a rigorous theoretical foundation that accounts for variability introduced by adaptively choosing $L^*$. To address nonregular estimation under the null, we develop two complementary strategies: a convex-hull-based estimator, and a penalized approach that further enhances power. Additionally, when restriction time candidates are defined on a discrete grid, we propose a procedure that surprisingly incurs no asymptotic penalty for selection, thus achieving oracle performance. 

## Installation

You can install the development version of **AdaRMST** from GitHub using the `devtools` package. First, ensure that you have `devtools` installed:

``` r
install.packages("devtools")
```
Then, install **AdaRMST** from GitHub:
``` r
devtools::install_github("Jinghao-Sun/AdaRMST")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(AdaRMST)

# Load the application dataset
data(pancreatic_cancer_year)

alpha = 0.05
B.boot = 1000
min_time_fix = 0.25
max_time_fix = 4.4

res.ct.pointest = AdaRMST.ct(
  pancreatic_cancer_year,
  min_time_fix,
  max_time_fix,
)

res.ct.ci = AdaRMST.ct.ci(
    pancreatic_cancer_year,
    alpha
    B.boot,
    min_time_fix,
    max_time_fix
)

res.dt = AdaRMST.dt(
    pancreatic_cancer_year,
    alpha,
    min_time_fix
    max_time_fix,
    L.grid = seq(min_time_fix, max_time_fix, length.out = 10)
) 

res.hulc = AdaRMST.hulc(
  pancreatic_cancer_year,
  alpha,
  type = "anti-conservative",
  min_time_fix,
  max_time_fix
)
```

## Included Dataset

We reconstructed a recent oncology clinical trial by (Tempero et al, 2023, Journal of Clinical Oncology). This randomized, open-label phase III trial evaluated the efficacy and safety of adjuvant nab-paclitaxel plus gemcitabine (nab-P + Gem) versus gemcitabine (Gem) alone in 866 patients with surgically resected pancreatic ductal adenocarcinoma. 

## License

This package is licensed under the GPL (>= 3) License.

## References

Jinghao Sun, Douglas E. Schaubel, Eric J. Tchetgen Tchetgen (2025). *Beyond Fixed Restriction Time: Adaptive Restricted Mean Survival Time Methods in Clinical Trials.*

# Contact

For questions, issues, or suggestions, please open an issue on GitHub or contact me (jinghao.sun.io\@gmail.com).
