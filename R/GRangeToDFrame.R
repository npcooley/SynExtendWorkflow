###### -- an adhoc function to improve data capture from GFFs -----------------
# author: nicholas cooley
# contact: npc19@pitt.edu / npcooley@gmail.com

GRangeToDFrame <- function(GRangeObject,
                           FeaturesToCollect = c("gene",
                                                 "pseudogene"),
                           Verbose = FALSE) {
  # search until no more searches are necessary
  s1 <- as.character(GRangeObject$type)
  s2 <- GRangeObject$gene
  s3 <- GRangeObject@strand
  s4 <- GRangeObject@ranges
  s5 <- as.character(GRangeObject@seqnames)
  s6 <- GRangeObject$Note
  s7 <- GRangeObject$Parent
  s8 <- GRangeObject$ID
  s9 <- !is.na(s2)
  
  TOTAL <- sum(table(s1[s1 %in% FeaturesToCollect]))
  # print(TOTAL)
  CONTINUE <- TRUE
  KEEP <- vector(mode = "logical",
                 length = length(s1))
  START <- STOP <- vector(mode = "integer",
                          length = length(s1))
  NOTE <- CONTIG <- TYPE <- ID <- GENE <- vector(mode = "character",
                                                 length = length(s1))
  COUNT <- 1L
  FOUNDFEATURES <- 0L
  if (Verbose) {
    TSTART <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  while (CONTINUE) {
    # is the line a line to evaluate
    # check its children
    if (s1[COUNT] %in% FeaturesToCollect) {
      if (s1[COUNT] == "pseudogene") {
        w1 <- which(s7 == s8[COUNT])
        w1 <- which(lengths(w1) > 0L)
        # print(w1)
        # if the feature has any children
        if (length(w1) > 0L) {
          ph1 <- ""
          for (m2 in seq_along(w1)) {
            ph2 <- unlist(s6[w1[m2]])
            # print(ph2)
            if (length(ph2) > 0) {
              if (!is.na(ph2)) {
                # print(nchar(ph1))
                if (nchar(ph1) == 0) {
                  ph1 <- ph2
                } else {
                  ph1 <- paste(ph1, ph2, sep = "; ")
                }
              }
            } else {
              ph1 <- "pseudofeature with absent note"
            }
          }
          NOTE[COUNT] <- ph1
        } else {
          # feature has no children, what to do here?
          NOTE[COUNT] <- "child lines absent"
        }
      } else {
        # ph1 <- "normal feature"
        NOTE[COUNT] <- "normal feature"
      }
      START[COUNT] <- s4@start[COUNT]
      STOP[COUNT] <- s4@start[COUNT] + s4@width[COUNT] - 1L
      CONTIG[COUNT] <- s5[COUNT]
      TYPE[COUNT] <- s1[COUNT]
      ID[COUNT] <- s8[COUNT]
      
      if (s9[COUNT]) {
        GENE[COUNT] <- s2[COUNT]
      } else {
        GENE[COUNT] <- ""
      }
      KEEP[COUNT] <- TRUE
      FOUNDFEATURES <- FOUNDFEATURES + 1L
      
    } # end if s2 is a feature to collect or not
    
    if (FOUNDFEATURES >= TOTAL) {
      CONTINUE <- FALSE
    } else {
      COUNT <- COUNT + 1L
    }
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = FOUNDFEATURES / TOTAL)
    }
    
  }
  if (Verbose) {
    close(pBar)
    cat("\n")
    TEND <- Sys.time()
    print(TEND - TSTART)
  }
  
  # return(list(START,
  #             STOP,
  #             TYPE,
  #             CONTIG,
  #             ID,
  #             NOTE,
  #             KEEP))
  res <- DataFrame("Start" = START[KEEP],
                   "Stop" = STOP[KEEP],
                   "Type" = TYPE[KEEP],
                   "Contig" = CONTIG[KEEP],
                   "ID" = ID[KEEP],
                   "Gene" = GENE[KEEP],
                   "Note" = NOTE[KEEP])
  return(res)
}
