#' used internally by fun_shreve_mag
#'
#' @param dwnAdj directed adjacency mtx
#' @param jAdded Vector of just added nodes
#' @param g graph object
#' @return integer - cut-off distance
#' @importFrom purrr map
fun_get_cut_off <- function(dwnAdj, jAdded, g) {
  tt <- expand.grid(dwnAdj, jAdded)
  xx <- mapply(igraph::all_simple_paths,
    from = tt$Var1,
    to = tt$Var2,
    MoreArgs = list(
      graph = g,
      mode = "in"))
  xx <- Filter(length, xx)
  min(unlist(purrr::map(purrr::map(xx, unlist), length)))
}


#' used internally by fun_shreve_mag
#'
#' @param g graph object
#' @return vector of "true" headwater nodes
#' @importFrom purrr map
fun_headW <- function(g) {
  h <- fun_headwater_nodes(g)
  dwnAdj <- unique(unlist(igraph::adjacent_vertices(g, v = h, mode = "out")))
  cutoff <- fun_get_cut_off(dwnAdj, h, g)
  for (i in dwnAdj) {
    #i = 9
    a <- igraph::all_simple_paths(g, from = i,  to = h, mode = "in")
    b <- max(unlist(purrr::map(a, length)))
    if (b > cutoff) # remove i from dwnAdj
      dwnAdj <- setdiff(dwnAdj, i)
  }
  unique(unlist(igraph::adjacent_vertices(g, v = dwnAdj, mode = "in")))
}

#' used internally by fun_shreve_mag
#'
#' @param x vector of 1+ length
#' @return vector
sample.vec <- function(x, ...) x[sample(length(x), ...)]


#' Create networks
#'
#' @param nNodes Intiger numbr of nodes.
#' @param shape Network shape, one of "lin", "int" or "den".
#' @return A list containing adjacency matrices and an igraph object.
#' @examples
#' \dontrun{
#' fun_crtNwk(10, "den")
#' }
#' @export
fun_crtNwk <- function(nNodes, shape = c("lin", "int", "den")) {

  #tr <- c(p1, p2, p3)

  if (shape == "lin") {
    tr <- c(0,1,0)
  }
  if (shape == "int") {
    tr <- runif(3, 0, 1)
  }

  if (shape == "den") {
    tr <- c(0,0,1)
  }

  if (nNodes == 2) {
    tr <- c(0, 1, 0)
  }

  # function to sample number of nodes to add
  fun_nNodesAdd <- function(br = c(0, 1, 2), times = 1) {
    sample(br, times, prob = tr, replace = TRUE)
  }

  # the + 10 is for over-spill, it will be culled later
  nodeLs <- vector(mode = "list", length = nNodes + 10)
  nodeLs <- setNames(nodeLs, 1:nNodes)

  # set first node
  node1 <- 1
  nodeLs[[1]] <- node1
  newNode <- nodeLs[[1]]

  # place holder for first round
  maxNode <- node1
  maxNode2 <- 0
  nodeVec2 <- vector()

  while (is.null(nodeLs[[nNodes]]) == TRUE) {
    # used to make sure all newly added nodes are remembered not
    # just the latest [i]teration in the for loop below
    preNodeVec <- unlist(nodeLs)

    if (length(newNode) == 0) {
      preMaxNode <- "dead"
    } else {
      preMaxNode <- maxNode
    }
    # if no new nodes have been added sample the "dead" nodes to find a new starting point
    if (preMaxNode == maxNode2) {
      NodeVec <- unlist(nodeLs)
      dNodes <- as.data.frame(subset(NodeVec, NodeVec == "dead"))
      newNode <- as.numeric(sample(rownames(dNodes), 1))
    } else {
      newNode <- newNode
    }
    if (preMaxNode == maxNode2) {
      for (i in newNode) {
        nNodesAdd <- fun_nNodesAdd()
        if (nNodesAdd > 0) {
          maxNode <- Position(function(nodeLs) !is.null(nodeLs), nodeLs, right = TRUE)
          nodeLs[[i]] <-  maxNode + 1:nNodesAdd
          nodeVec2 <- unlist(nodeLs)
          newNode <- nodeLs[[i]]
          nodeID <- max(newNode)
        } else {
          nodeLs[[i]] <- "dead"
        }
      }
    } else {
      for (i in min(newNode):max(newNode)) {
        nNodesAdd <- fun_nNodesAdd()
        if (nNodesAdd > 0) {
          maxNode <- max(newNode)
          nodeLs[[i]] <-  maxNode + 1:nNodesAdd
          nodeVec2 <- unlist(nodeLs)
          newNode <- nodeLs[[i]]
          nodeID <- max(newNode)
        } else {
          nodeLs[[i]] <- "dead"
        }
      }
    }

    # determine which new nodes have been added
    findNewNodes <- setdiff(nodeVec2, preNodeVec)
    newNode <- as.numeric(findNewNodes[!findNewNodes %in% "dead"])
    if (length(newNode) == 0) {
      maxNode2 <- "dead"
    } else {
      maxNode2 <- max(newNode)
    }
  }
  # convert nodeLs to numeric and "dead" to NA's
  sNodeLs <- nodeLs[1:nNodes]
  sNodeVec <- unname(sNodeLs)
  suppressWarnings(sNodeLs2 <- lapply(sNodeVec, as.numeric))

  # convert nodeLs to adjacency matrix
  if (nodeID >= nNodes) {
    adj_mtx <- matrix(0, nrow = nNodes, ncol = nNodes)
    for (i in 1:length(sNodeLs2)) {
      filt <- sNodeLs2[[i]] <= nNodes
      adj_mtx[i, sNodeLs2[[i]][filt]] <- 1
    }
    # adjacency to igraph structure
    g <- igraph::graph_from_adjacency_matrix(t(adj_mtx), mode = "directed")


    # convert adjaceny matrix to useable form for metasteBySp model
    unweighted_mtx <- t(adj_mtx) + adj_mtx

    # weighted adjacency matrix (downstream = 1, upstream = 2)
    weighted_mtx <- unweighted_mtx + adj_mtx

    ig <- igraph::graph.adjacency(weighted_mtx, mode = "directed", weighted = TRUE)
    nodeToNodeDist <- igraph::shortest.paths(ig)

    # where disperses can go.
    dispersalDf  <- igraph::graph.adjacency(weighted_mtx, weighted = TRUE)
    dispersalDf <- igraph::get.data.frame(dispersalDf)
    g2 <- igraph::graph.data.frame(dispersalDf, directed = TRUE)
    wgtNodetoNodeDist <- igraph::distances(g2, v = igraph::V(g2), to = igraph::V(g2),
      mode = "out", weights = dispersalDf$direct,
      algorithm = "dijkstra")

    ## store all useful network info
    list(edges = sNodeLs2, adjUp = adj_mtx, adjDwn = t(adj_mtx), nwk = g, tr = tr,
      unweighted_mtx = unweighted_mtx, weighted_mtx = weighted_mtx,
      ig = ig, nodeToNodeDist = nodeToNodeDist, wgtNodetoNodeDist = wgtNodetoNodeDist)
  } else {
    list(edges = NULL, adjUp = NULL, adjDwn = NULL, nwk = NULL, tr = NULL,
      unweighted_mtx = NULL, weighted_mtx = NULL, ig = NULL,
      nodeToNodeDist = NULL, wgtNodetoNodeDist = NULL)
  }
}



#' Rescale vector to desired range.
#'
#' @param x vector
#' @param L Lower limit
#' @param H Upper limit
#' @return A vector
#' @examples
#' \dontrun{
#' rescale.to.range2(seq(1, 20), 100, 200)
#' }
#' @export
rescale.to.range2 <- function(x, L, H)
{
  x.range <- range(x)
  if (x.range[1] == x.range[2]) {return(x)}
  mfac <- (H - L)/(x.range[2] - x.range[1])
  return(L + (x - x.range[1]) * mfac)
}


#' Generate noise signal
#'
#' @param alpha autocorrelation coefficient
#' @param beta xx
#' @param C xx
#' @param L Lower limit for rescaling, default = 0
#' @param H Upper limit for rescaling, default = 1
#' @param num_gens sequence length
#' @return A vector
#' @examples
#' \dontrun{
#' generate.noise(alpha = 0, beta = 1, C = 1, L = 0, H = 1, num_gens = 100)
#' }
#' @export
generate.noise <- function(alpha = 0, beta = 1, C = 1, L = 0, H = 1, num_gens) {
  series <- vector(length = num_gens, mode = 'numeric')
  err.series <- vector(length = num_gens, mode = 'numeric')
  # Rescale to the actual variance (Eq 6-8 in)
  series[1] <- 0.0
  for (i in 2:num_gens)
  {
    err.series[i] <- rnorm(1)
    series[i] <- alpha * series[i - 1] + beta * err.series[i]
  }
  series <- (sqrt(C) / sd(series)) * (series - mean(series))
  series <- rescale.to.range2(series, L, H)
  series
}


#' reverse direction of a directed graph
#'
#' @param g a dendritic network as an igraph object
#' @return an igraph object
#' @examples
#' \dontrun{
#' fun_reverse_graph(g)
#' }
#' @export
fun_reverse_graph <- function(g) {
  if (!igraph::is.directed(g)) stop("Graph must be directed.")
  el <- igraph::get.edgelist(g, names = FALSE)
  graph(rbind(el[, 2], el[, 1]))
}

#' find all upstream nodes of a vertex
#' taken and edited from https://github.com/robertness/lucy/blob/master/R/lucy.R
#'
#' @param g a dendritic network as an igraph object
#' @param w target node
#' @return a vector of upstream node id's
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' fun_upstream_nodes(g, 5)
#' }
#' @export
fun_upstream_nodes <- function(g, w){
  if (!igraph::is.directed(g)) stop("Graph must be directed.")
  igraph::V(g)$name <- paste(igraph::V(g))
  sp.mat <- igraph::shortest.paths(g, v = igraph::V(g), to = w, mode = "out")
  if (is.null(dimnames(sp.mat))) {
    dimnames(sp.mat) <- list(paste(1:igraph::vcount(g)), paste(w))
  }
  sp.mat[, igraph::V(g)[w]$name] %>%
    Filter(f = is.finite) %>%
    names %>%
    {igraph::V(g)[.]} %>%
    as.numeric %>%
    setdiff(w)
}


#' find all downstream nodes of a vertex
#' taken and edited from https://github.com/robertness/lucy/blob/master/R/lucy.R
#'
#' @param g a dendritic network as an igraph object
#' @param w target node
#' @return a vector of downstream node id's
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' fun_downstream_nodes(g, 5)
#' }
#' @export
fun_downstream_nodes <- function(g, w){
  if (!igraph::is.directed(g)) stop("Graph must be directed.")
  igraph::V(g)$name <- paste(igraph::V(g))
  sp.mat <- igraph::shortest.paths(g, v = igraph::V(g), to = w, mode = "in")
  if (is.null(dimnames(sp.mat))) {
    dimnames(sp.mat) <- list(paste(1:igraph::vcount(g)), paste(w))
  }
  sp.mat[, igraph::V(g)[w]$name] %>%
    Filter(f = is.finite) %>%
    names %>%
    {igraph::V(g)[.]} %>%
    as.numeric %>%
    setdiff(w)
}

#' get all headwater nodes
#'
#' @param g a dendritic network as an igraph object
#' @return a vector of headwater node id's
#' @examples
#' \dontrun{
#' fun_headwater_nodes(g)
#' }
#' @export
fun_headwater_nodes <- function(g){
  if (!igraph::is.directed(g)) stop("Graph must be directed.")
  which(igraph::degree(g, v = igraph::V(g), mode = "in") == 0)
}


#' Get strahler order
#'
#' @param g a dendritic network as an igraph object
#' @return a vector of stream strahler orders
#' @examples
#' \dontrun{
#' fun_strahler_order(g)
#' }
#' @export
fun_strahler_order <- function(g) {
  # go from starting graph to directed towards the root
  g <- igraph::as_adjacency_matrix(igraph::as.undirected(g), type = "lower")
  g <- igraph::graph.adjacency(g)

  # set-up
  df <- data.frame(node = 1:igraph::gorder(g), d = igraph::degree(g, mode = "in"), strahler = 0)
  df$strahler <- ifelse(df$d == 0, 1, 0)
  df <- df[,c(1,3)]
  strt <- df[df$strahler == max(df$strahler),"node"]

  while(any(df$strahler == 0)) {

    sav <- unique(unlist(igraph::adjacent_vertices(g, strt, mode = "out")))

    for( i in seq_along(strt)) {

      # get downstream nodes (start with highest numbered node)
      xx <- unlist(igraph::adjacent_vertices(g, rev(strt)[i], mode = "out"))

      # upstream nodes
      nd2 <- unlist(igraph::adjacent_vertices(g, xx, mode = "in"))

      # order of upstream nodes
      strhl <- df[nd2, "strahler"]

      # figure out new strahler order
      res <- NA
      if (length(strhl) == 1) res <- strhl
      if (length(strhl) > 1) {
        res <- ifelse(strhl[1] == strhl[2], max(strhl) + 1, max(strhl))
      }
      #update dataframe
      df[xx, "strahler"] <- res

    }
    strt <- sav
  }
  df$strahler
}

#' Shreves stream magnitude
#'
#' @param g a dendritic network as an igraph object
#' @return an igraph object with stream Shreve's order added as an attribute
#' @export
fun_shreve_mag <- function(g) {
  y <- fun_headwater_nodes(g)
  z <- sapply(1:igraph::gorder(g), function(x) {
    x <- fun_upstream_nodes(g, x)
    length(intersect(y,x))
    }
    )
  z <- ifelse(z == 0, 1, z)
  g <- igraph::set_vertex_attr(g, "shreve", value = z)
  g
}


#' Quick and dirty plotting of dendritic networks, wrapper for
#' igraph::layout_as_tree and plot.igraph.
#'
#' @param g a dendritic network as an igraph object
#' @return plot
#' @examples
#' \dontrun{
#' fun_pltNwk(g)
#' }
#' @export
fun_pltNwk <- function(g, direct = c("in", "out"), ...) {
  l <- suppressWarnings(igraph::layout_as_tree(g, flip.y = FALSE, mode = direct))
  plot(g,
    vertex.label.color = "black",
    vertex.color = "darkgray",
    vertex.frame.color = " white",
    vertex.shape = "circle",
    edge.width = 5,
    # edge.arrow.size = 2,
    edge.color = "grey",
    layout = l,
    asp = 0,
    ...)
}

