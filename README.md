
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

A few functions used for creating and manipulating dendritic networks,
mainly wrappers for various igraph functions.

### Installation

``` r
# devtools::install_github("flee598/DenMet")
library(DenMet)
```

### Create a dendritic network with 10 nodes

``` r
nwk <- DenMet::fun_crtNwk2(10, 0.5)
nwk2 <- igraph::graph.adjacency(nwk)
fun_pltNwk(nwk2, direct = "out", edge.arrow.size = 2)
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

### Downstream network, hide arrows

``` r
nwk <- DenMet::fun_convert_nwk(nwk, "dwn_adj")
g <- igraph::graph.adjacency(nwk)
fun_pltNwk(g, "in", edge.arrow.size = 0)
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

### get strahler order of graph nodes

``` r
fun_strahler_order(g)
#>  [1] 2 2 2 2 1 2 1 2 1 1
```

### Get downstream and upstream nodes of a given node

``` r
fun_downstream_nodes(g, 8)
#> [1] 1 2 3 4 6
fun_upstream_nodes(g, 2)
#> [1]  3  4  5  6  7  8  9 10
```
