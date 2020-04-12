
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

A few functions used for creating and manipulating dendritic networks.

### Installation

``` r
# devtools::install_github("flee598/DenMet")
library(DenMet)
```

### create a dendritic network with 10 nodes

``` r
nwk <- fun_crtNwk(10, "int")
plot(nwk[[1]]$ig)
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->
