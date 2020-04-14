
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

A few functions used for creating and manipulating dendritic networks,
mainly wrappers for various igraph functions.

### Installation

``` r
devtools::install_github("flee598/DenMet")
#> 
#>          checking for file 'C:\Users\Finn\AppData\Local\Temp\RtmpO2vjA2\remotes376839f46d87\flee598-DenMet-283193d/DESCRIPTION' ...  v  checking for file 'C:\Users\Finn\AppData\Local\Temp\RtmpO2vjA2\remotes376839f46d87\flee598-DenMet-283193d/DESCRIPTION'
#>       -  preparing 'DenMet':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   v  checking DESCRIPTION meta-information
#>       -  checking for LF line-endings in source and make files and shell scripts
#>   -  checking for empty or unneeded directories
#>       -  building 'DenMet_0.1.0.tar.gz'
#>      
#> 
library(DenMet)
library(igraph)
```

### Create a dendritic network with 10 nodes

``` r
nwk <-  fun_crtNwk(10, "int")
fun_pltNwk(nwk$ig, direct = "out", edge.arrow.size = 2)
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

# downstream network, hide arrows

``` r
g <- nwk$adjDwn
g <- igraph::graph.adjacency(g)
fun_pltNwk(g, "in", edge.arrow.size = 0)
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

# add strahler order graph nodes

``` r
g <- fun_strahler_order(g)
igraph::get.vertex.attribute(g)
#> $strahler
#>  [1] 3 2 2 2 1 1 1 1 1 1
```
