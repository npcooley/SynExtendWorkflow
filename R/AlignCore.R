###### -- do some alignments and other things ---------------------------------

AlignCore <- function(GeneCalls,
                      PotentialCoreGenes,
                      Verbose = FALSE) {
  
  CoreAli <- CoreAnnotations <- CoreDends <- vector(mode = "list",
                                                    length = length(PotentialCoreGenes))
  if (Verbose) {
    TSTART <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
    PBAR <- length(CoreAli)
  }
  for (m1 in seq_along(CoreAli)) {
    w1 <- do.call(rbind,
                  strsplit(x = names(PotentialCoreGenes[[m1]]),
                           split = "_",
                           fixed = TRUE))
    w1 <- matrix(data = as.integer(w1),
                 nrow = nrow(w1))
    w2 <- mapply(FUN = function(x, y) {
      GeneCalls[[x]]$Coding[y]
      # attr(x = p01,
      #      which = "GeneCalls")[[x]]$Coding[y]
    },
    x = w1[, 1L],
    y = w1[, 3L])
    CoreAnnotations[[m1]] <- mapply(FUN = function(x, y) {
      GeneCalls[[x]]$Product[y]
      # attr(x = p01,
      #      which = "GeneCalls")[[x]]$Product[y]
    },
    x = w1[, 1L],
    y = w1[, 3L])
    CoreAli[[m1]] <- PotentialCoreGenes[[m1]][order(w1[, 1L])]
    if (all(w2)) {
      CoreAli[[m1]] <- AlignTranslation(myXStringSet = CoreAli[[m1]],
                                        verbose = FALSE)
      x <- DistanceMatrix(CoreAli[[m1]],
                          includeTerminalGaps = TRUE,
                          verbose = FALSE)
      CoreDends[[m1]] <- TreeLine(myDistMatrix = x,
                                  method = "NJ",
                                  type = "dendrogram",
                                  verbose = FALSE)
    } else {
      CoreAli[[m1]] <- AlignSeqs(CoreAli[[m1]],
                                 verbose = FALSE)
      x <- DistanceMatrix(CoreAli[[m1]],
                          includeTerminalGaps = TRUE,
                          verbose = FALSE)
      CoreDends[[m1]] <- TreeLine(myDistMatrix = x,
                                  method = "NJ",
                                  type = "dendrogram",
                                  verbose = FALSE)
    }
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / PBAR)
    }
  }
  
  if (Verbose) {
    close(pBar)
    cat("\n")
    TEND <- Sys.time()
    print(TEND - TSTART)
  }
  
  return(list("Alignments" = CoreAli,
              "Dendrograms" = CoreDends,
              "Annotations" = CoreAnnotations))
}
