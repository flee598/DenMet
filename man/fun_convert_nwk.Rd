% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networks.R
\name{fun_convert_nwk}
\alias{fun_convert_nwk}
\title{Convert SBN adjacency matrices to other useful formats.}
\usage{
fun_convert_nwk(
  adj_mtx,
  out = c("dwn_adj", "unwtd_adj", "wtd_adj", "ig", "dwn_ig", "unwtd_ig", "n2n_dist",
    "unwtd_n2n_dist")
)
}
\arguments{
\item{adj_mtx}{SBN mtx}

\item{out}{One of 8 possible output formats}
}
\value{
An adjacency matrix
}
\description{
Convert SBN adjacency matrices to other useful formats.
}
\examples{
\dontrun{
fun_convert_nwk(mtx, "ig")
}
}
