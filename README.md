# MultiMAP

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

Simultaneous credible intervals for multi-arm studies that borrow historical
controls or other external data via robustified meta-analytic-predictive (MAP)
priors. `MultiMAP` also provides MAP-based prediction intervals for a quality
check of the concurrent control, as desired by many OECD test guidelines for
regulatory toxicology studies. Normal, binomial and Poisson endpoints are
supported.

> **Developmental version.** Especially the methodology for Poisson endpoints
> needs further adaptions and is not yet ready for practical usage.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("MaxMenssen/MultiMAP")
```

## Example

Borrow historical control data into the concurrent control of a new trial and
compute simultaneous credible intervals for the dose-vs-control risk ratios
(binomial endpoint):

``` r
library(MultiMAP)

# Historical control data (one row per study)
histdat <- data.frame(events     = c(2, 2, 1, 5, 3, 5, 2, 2, 4, 2, 3, 3, 3, 2, 3, 4),
                      non_events = c(48,47,49,44,51,44,54,53,46,46,46,86,71,45,45,46))

# Current trial: row 1 is the concurrent control, the remaining rows are dose groups
newdat  <- data.frame(group      = c("0", "10", "50", "1500", "10000"),
                      events     = c(3, 4, 5, 7, 8),
                      non_events = c(44, 36, 39, 42, 41))

res <- SCI_binom(histdat, newdat, prior_est = "MAP", contr = "RR", seed = 81771)
summary(res)
#> Single-arm borrowing for the concurrent control
#> Arm '0'
#>   95 % prediction interval for the proportion: [0, 0.1277]
#>   observed proportion: 0.06383 (covered)
#>   ESS  MAP: 316 | robust. MAP: 243 | posterior: 369
#> 95 % simultaneous credible intervals for ratio of proportions
#>      comp estimate lower.ci upper.ci
#> 1    10/0     2.07    0.504     6.08
#> 2    50/0     2.30    0.595     6.44
#> 3  1500/0     2.79    0.888     7.26
#> 4 10000/0     3.16    1.050     7.93
```

Analogous functions cover normal (`SCI_norm`) and Poisson (`SCI_pois`) endpoints,
prediction-interval-only variants (`PI_norm`, `PI_binom`, `PI_pois`), and
borrowing for more than one arm (pass `histdat` as a named list). See
`?SCI_binom` for the full argument list.

## Reference

Menssen, M., Kneuer, C., Akyianu, G., Roever, C., Friede, T., & Schaarschmidt, F.
(2026). Including historical control data in simultaneous inference for
pre-clinical multi-arm studies. *arXiv preprint* arXiv:2603.11730.
<https://doi.org/10.48550/arXiv.2603.11730>

## License

GPL-3
