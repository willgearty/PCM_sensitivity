# these are custom versions of FossilSim::sim.fossils.poisson() and
# FossilSim::sim.taxonomy() that use data.table::rbindlist for speed

library(TreeSim)
library(FossilSim)
sim.fossils.poisson.data.table <- function (rate, tree = NULL, taxonomy = NULL, fossils = NULL,
                                            ignore.taxonomy = FALSE, root.edge = TRUE)
{
  if (is.null(tree) && is.null(taxonomy))
    stop("Specify phylo or taxonomy object")
  if (!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if (!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
    stop("taxonomy must be an object of class \"taxonomy\"")
  if (!is.null(fossils) && !"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")
  if (!is.null(tree) && !is.null(taxonomy))
    warning("tree and taxonomy both defined, using taxonomy")
  if (!is.null(attr(rate, "from.taxonomy"))) {
    if (attr(rate, "from.taxonomy") && is.null(taxonomy)) {
      stop("rates simulated from taxonomy, matching \"taxonomy\" object also required")
    }
    if (!attr(rate, "from.taxonomy") && is.null(tree)) {
      stop("rates simulated from tree, matching \"tree\" object also required")
    }
  }
  if (is.null(taxonomy) && is.null(tree$edge.length))
    stop("tree must have edge lengths")
  if (is.null(taxonomy) && !ape::is.rooted(tree))
    stop("tree must be rooted")
  if (is.null(taxonomy)) {
    taxonomy = sim.taxonomy.data.table(tree, beta = 1, root.edge = root.edge) # WG: use data.table version
    if (length(rate) > 1) {
      if (is.null(tree$root.edge))
        rate = c(0, rate)
      rate = rate[order(c(root(tree), tree$edge[, 2]))]
      rate = rate[as.numeric(taxonomy$sp)]
    }
    from.taxonomy = FALSE
  }
  else from.taxonomy = TRUE
  if (!all(as.vector(na.omit(fossils$edge)) %in% taxonomy$edge))
    stop("Mismatch between fossils and taxonomy objects")
  if (length(rate) > 1 && length(rate) != length(unique(taxonomy$sp)))
    stop("The vector of rates provided doesn't correspond to the number of species")
  else if (length(rate) == 1)
    rate = rep(rate, length(unique(taxonomy$sp)))
  if (any(rate < 0))
    stop("Rates must be positive numbers")
  use.exact.times = TRUE
  if (is.null(fossils))
    fdf = list()
  else fdf = fossils
  lineages = unique(taxonomy$sp)
  for (i in 1:length(lineages)) {
    sp = lineages[i]
    start = max(taxonomy$start[which(taxonomy$sp == sp)])
    end = min(taxonomy$end[which(taxonomy$sp == sp)])
    edges = taxonomy[which(taxonomy$sp == sp), ] # WG: this subsetting would break if taxonomy were a data.table
    blength = start - end
    rand = rpois(1, blength * rate[i])
    if (rand > 0) {
      if (ignore.taxonomy)
        sp = NA
      h = runif(rand, min = end, max = start)
      edge = sapply(h, function(x) edges$edge[which(edges$start >
                                                      x & edges$end < x)])
      if (use.exact.times) {
        # WG: this was originally rbind()
        fdf[[i]] <- data.frame(sp = sp, edge = edge,
                               hmin = h, hmax = h, stringsAsFactors = F)
      }
      else {
        # WG: this was originally rbind()
        fdf[[i]] <- data.frame(sp = sp, edge = edge,
                               hmin = rep(end, rand), hmax = rep(start, rand),
                               stringsAsFactors = F)
      }
    }
  }
  fdf <- data.table::rbindlist(fdf) # WG: added this to combine the list of data.frames
  fdf <- as.fossils(fdf, from.taxonomy)
  return(fdf)
}

sim.taxonomy.data.table <- function (tree, beta = 0, lambda.a = 0, kappa = 0, root.edge = TRUE)
{
  if (!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if (!(beta >= 0 && beta <= 1))
    stop("beta must be a probability between 0 and 1")
  if (lambda.a < 0)
    stop("lambda.a must be zero or positive")
  if (!(kappa >= 0 && kappa <= 1))
    stop("kappa must be a probability between 0 and 1")
  if (is.null(tree$edge.length))
    stop("tree must have edge lengths")
  if (!ape::is.rooted(tree))
    stop("tree must be rooted")
  node.ages = FossilSim:::n.ages(tree)
  species <- data.frame(sp = integer(), edge = integer(), parent = integer(),
                        start = numeric(), end = numeric(), mode = character(),
                        cryptic = logical(), cryptic.id = integer())
  root = length(tree$tip.label) + 1
  if (root.edge && exists("root.edge", tree)) {
    start = node.ages[root] + tree$root.edge
    mode = "o"
  }
  else {
    start = node.ages[root]
    mode = "r"
  }
  species <- rbind(species, data.frame(sp = root, edge = root,
                                       parent = 0, start = start, end = node.ages[root], mode = mode,
                                       cryptic = 0, cryptic.id = root))
  aux = function(node, p) {
    descendants = tree$edge[which(tree$edge[, 1] == node),
                            2]
    if (length(descendants) == 0) {
      return(p)
    }
    d1 <- descendants[1]
    d2 <- descendants[2]
    if (beta == 1 || (beta > 0 && runif(1) > (1 - beta))) {
      a <- p$sp[which(p$edge == node)]
      p <- data.table::rbindlist(list( # WG: this was originally rbind()
        p,
        data.frame(sp = d1, edge = d1, parent = a,
                   start = node.ages[a], end = node.ages[d1], mode = "s",
                   cryptic = 0, cryptic.id = d1),
        data.frame(sp = d2, edge = d2, parent = a,
                   start = node.ages[a], end = node.ages[d2], mode = "s",
                   cryptic = 0, cryptic.id = d2)
        ))
    }
    else {
      p$sp[which(p$sp == node)] = d1
      p$parent[which(p$parent == node)] = d1
      p$cryptic.id[which(p$cryptic.id == node)] = d1
      a <- p$parent[which(p$edge == node)]
      m <- p$mode[which(p$sp == d1)][1]
      p <- data.table::rbindlist(list( # WG: this was originally rbind()
        p,
        data.frame(sp = d1, edge = d1, parent = a,
                   start = node.ages[FossilSim:::ancestor(d1, tree)], end = node.ages[d1],
                   mode = m, cryptic = 0, cryptic.id = d1),
        data.frame(sp = d2, edge = d2, parent = d1,
                   start = node.ages[FossilSim:::ancestor(d2, tree)], end = node.ages[d2], mode = "b",
                   cryptic = 0, cryptic.id = d2)
        ))
    }
    p = aux(d1, p)
    p = aux(d2, p)
    p
  }
  species = as.data.frame(aux(root, species)) # WG: convert back to data.frame
  species = species[order(species$sp), ]
  species = taxonomy(species)
  if (lambda.a > 0)
    species = sim.anagenetic.species(tree, species, lambda.a)
  if (kappa > 0)
    species = sim.cryptic.species(species, kappa)
  rownames(species) = NULL
  return(species)
}
