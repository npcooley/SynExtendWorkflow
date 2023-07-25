###### -- BlockSize -----------------------------------------------------------
# Given a PairSummaries object return a vector of integers that gives
# the size of the syntenic block that each pair is in

BlockSize <- function(Pairs,
                      Verbose) {
  if (!is(object = Pairs,
          class2 = "PairSummaries")) {
    stop("Object must be of class PairSummaries")
  }
  if (nrow(Pairs) < 1L) {
    stop("Must include at least one pair.")
  }
  
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  POIDs <- paste(Pairs$p1,
                 Pairs$p2,
                 sep = "_")
  FeaturesMat <- do.call(rbind,
                         strsplit(x = POIDs,
                                  split = "_",
                                  fixed = TRUE))
  FeaturesMat <- data.frame("g1" = as.integer(FeaturesMat[, 1L]),
                            "i1" = as.integer(FeaturesMat[, 2L]),
                            "f1" = as.integer(FeaturesMat[, 3L]),
                            "g2" = as.integer(FeaturesMat[, 4L]),
                            "i2" = as.integer(FeaturesMat[, 5L]),
                            "f2" = as.integer(FeaturesMat[, 6L]))
  dr1 <- FeaturesMat[, 3L] + FeaturesMat[, 6L]
  dr2 <- FeaturesMat[, 3L] - FeaturesMat[, 6L]
  InitialBlocks1 <- unname(split(x = FeaturesMat,
                                 f = list(as.integer(FeaturesMat[, 1L]),
                                          as.integer(FeaturesMat[, 4L]),
                                          as.integer(FeaturesMat[, 2L]),
                                          as.integer(FeaturesMat[, 5L]),
                                          dr1),
                                 drop = TRUE))
  InitialBlocks2 <- unname(split(x = FeaturesMat,
                                 f = list(as.integer(FeaturesMat[, 1L]),
                                          as.integer(FeaturesMat[, 4L]),
                                          as.integer(FeaturesMat[, 2L]),
                                          as.integer(FeaturesMat[, 5L]),
                                          dr2),
                                 drop = TRUE))
  
  Blocks <- c(InitialBlocks1[sapply(InitialBlocks1,
                                    function(x) nrow(x),
                                    simplify = TRUE) > 1],
              InitialBlocks2[sapply(InitialBlocks2,
                                    function(x) nrow(x),
                                    simplify = TRUE) > 1])
  L01 <- length(Blocks)
  
  for (m1 in seq_along(Blocks)) {
    # blocks are guaranteed to contain more than 1 row
    
    sp1 <- vector(mode = "integer",
                  length = nrow(Blocks[[m1]]))
    sp2 <- Blocks[[m1]][, 3L]
    
    it1 <- 1L
    it2 <- sp2[1L]
    # create a map vector on which to split the groups, if necessary
    for (m2 in seq_along(sp1)) {
      it3 <- sp2[m2]
      if (it3 - it2 > 1L) {
        # if predicted pairs are not contiguous, update the iterator
        it1 <- it1 + 1L
      }
      sp1[m2] <- it1
      it2 <- it3
    }
    
    # if the splitting iterator was updated at all, a gap was detected
    if (it1 > 1L) {
      Blocks[[m1]] <- unname(split(x = Blocks[[m1]],
                                   f = sp1))
    } else {
      Blocks[[m1]] <- Blocks[m1]
    }
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / L01)
    }
  }
  if (Verbose) {
    close(pBar)
    cat("\n")
  }
  Blocks <- unlist(Blocks,
                   recursive = FALSE)
  
  Blocks <- Blocks[sapply(X = Blocks,
                          FUN = function(x) {
                            nrow(x)
                          },
                          simplify = TRUE) > 1L]
  L01 <- length(Blocks)
  
  AbsBlockSize <- rep(1L,
                      nrow(Pairs))
  
  for (m1 in seq_along(Blocks)) {
    pos <- as.integer(rownames(Blocks[[m1]]))
    val <- rep(nrow(Blocks[[m1]]),
               nrow(Blocks[[m1]]))
    keep <- AbsBlockSize[pos] < val
    
    if (any(keep)) {
      AbsBlockSize[pos[keep]] <- val[keep]
    }
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / L01)
    }
  }
  if (Verbose) {
    close(pBar)
    cat("\n")
  }
  
  if (Verbose) {
    TimeEnd <- Sys.time()
    # close(pBar)
    cat("\nBlock sizes identified.\n")
    print(TimeEnd - TimeStart)
  }
  
  return(AbsBlockSize)
}

