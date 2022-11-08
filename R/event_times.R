#' Times at which speciation or extinction occurs
#' @title Event times of a (possibly non-ultrametric) phylogenetic tree
#' @param phy phylogenetic tree of class phylo, without polytomies, rooted and 
#' with branch lengths. Need not be ultrametric.
#' @return times at which speciation or extinction happens.
#' @note This script has been modified from BAMMtools' internal function 
#' NU.branching.times
#' @export
event_times <- function(phy) {
    if (ape::is.ultrametric(phy)) {
        return(ape::branching.times(phy))
    } else {
        if (ape::is.binary(phy) == FALSE) {
            stop("error. Need fully bifurcating (resolved) tree\n")
        }
        phy$begin <- rep(0, nrow(phy$edge))
        phy$end <- rep(0, nrow(phy$edge))
        fx <- function(phy, node) {
            cur.time <- 0
            root <- length(phy$tip.label) + 1
            if (node > root) {
                cur.time <- phy$end[which(phy$edge[, 2] == node)]
            }
            dset <- phy$edge[, 2][phy$edge[, 1] == node]
            i1 <- which(phy$edge[, 2] == dset[1])
            i2 <- which(phy$edge[, 2] == dset[2])
            phy$end[i1] <- cur.time + phy$edge.length[i1]
            phy$end[i2] <- cur.time + phy$edge.length[i2]
            if (dset[1] > length(phy$tip.label)) {
                phy$begin[phy$edge[, 1] == dset[1]] <- phy$end[i1]
                phy <- fx(phy, node = dset[1])
            }
            if (dset[2] > length(phy$tip.label)) {
                phy$begin[phy$edge[, 1] == dset[2]] <- phy$end[i2]
                phy <- fx(phy, node = dset[2])
            }
            return(phy)
        }
        phy <- fx(phy, node = length(phy$tip.label) + 1)
        maxbt <- max(phy$end)
        nodes <- (length(phy$tip.label) + 1):(2 * length(phy$tip.label) - 1)
        bt <- numeric(length(nodes))
        names(bt) <- nodes
        for (i in seq_along(bt)) {
            tt <- phy$begin[phy$edge[, 1] == nodes[i]][1]
            bt[i] <- maxbt - tt
        }
        return(bt)
    }
}