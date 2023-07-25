README
================
npc
2023-07-25

# Working With SynExtend

The SynExtend R package was originally built as a set of tools for
building out lists of orthologous gene pairs from biological data,
though it has since expanded. It is currently a work in progress. The
provided workflow here is a simple example built to showcase usage and
outputs of the orthology prediction tools housed within. The most
unfinished part of this workflow is in the `SelectByK` function that
uses a naive k-means based approach to drop assumed false positive
pairs.

## Initial workspace setup

This workflow is loading in two ad-hoc functions, one is an updated
version of a function within `SynExtend`, the other is an additional
function for providing context to orthologous pairs.

``` r
# start the timer
TIMESTART <- Sys.time()

# load in some libraries 
suppressMessages(library(SynExtend))
suppressMessages(library(igraph))
# SelectByK has had some minor adjustments, this imported version was added to v 1.13.5 (devel)
source(file = "R/SelectByK.R",
       echo = FALSE)
# The block size adjacency observer in PairSummaries() will eventually be superceded by this function
source(file = "R/BlockSize.R",
       echo = FALSE)
set.seed(1986)

# Sasha's distinct colors, hexcodes,
# from:
# https://sashamaps.net/docs/resources/20-colors/
ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")
ColVec3 <- paste0(ColVec1,
                  "50")
# specifically for the overlapping histogram
ColVec4 <- ColVec3[c(1:6, 13, 15, 18)]
```

## Data Collection Part 1

The NCBI Edirect tools can be used to collect small, or large sets of
data to evaluate. In this case a small example set of bacterial genomes
and their annotations can be selected.

``` r
# using NCBI edirect tools to grab a small digestible test set
Selection <- "deinococcus"
EntrezQuery <- paste("esearch -db assembly ",
                     "-query '",
                     Selection,
                     "[organism] ",
                     'AND "complete genome"[filter] ',
                     'AND "refseq has annotation"[properties] ',
                     "NOT anomalous[filter]' ",
                     '| ',
                     'esummary ',
                     '| ',
                     'xtract -pattern DocumentSummary -element FtpPath_RefSeq',
                     sep = "")
FtPPaths <- system(command = EntrezQuery,
                   intern = TRUE,
                   timeout = 500L)
if (length(FtPPaths) > 5L) {
  FtPPaths <- FtPPaths[sample(x = length(FtPPaths),
                              size = 5L,
                              replace = FALSE)]
}
```

## Data Collection Part 2

Pull the FNA and GFF files associated with the selected genomes, build a
DECIPHER DB, and manage the gene calls.

``` r
# grab associated FNAs and GFFs
FNAs <- unname(sapply(FtPPaths,
                      function(x) paste(x,
                                        "/",
                                        strsplit(x,
                                                 split = "/",
                                                 fixed = TRUE)[[1]][10],
                                        "_genomic.fna.gz",
                                        sep = "")))
GFFs <- unname(sapply(FtPPaths,
                      function(x) paste(x,
                                        "/",
                                        strsplit(x,
                                                 split = "/",
                                                 fixed = TRUE)[[1]][10],
                                        "_genomic.gff.gz",
                                        sep = "")))

GC01 <- vector(mode = "list",
               length = length(FNAs))

DBPATH <- paste0(getwd(),
                 "/data/",
                 Selection,
                 ".sqlite")

for (m1 in seq_along(FNAs)) {
  
  
  X <- readDNAStringSet(FNAs[m1])
  
  Seqs2DB(seqs = X,
          type = "XStringSet",
          dbFile = DBPATH,
          identifier = as.character(m1),
          verbose = TRUE)
  
  GC01[[m1]] <- gffToDataFrame(GFF = GFFs[m1],
                               Verbose = TRUE)
}
```

    ## Adding 7 sequences to the database.
    ## 
    ## 7 total sequences in table Seqs.
    ## Time difference of 0.14 secs
    ## 
    ## ================================================================================
    ## Time difference of 25.19963 secs
    ## Adding 1 sequences to the database.
    ## 
    ## Added 1 new sequence to table Seqs.
    ## 8 total sequences in table Seqs.
    ## Time difference of 0.11 secs
    ## 
    ## ================================================================================
    ## Time difference of 27.08691 secs
    ## Adding 1 sequences to the database.
    ## 
    ## Added 1 new sequence to table Seqs.
    ## 9 total sequences in table Seqs.
    ## Time difference of 0.09 secs
    ## 
    ## ================================================================================
    ## Time difference of 27.2416 secs
    ## Adding 7 sequences to the database.
    ## 
    ## Added 7 new sequences to table Seqs.
    ## 16 total sequences in table Seqs.
    ## Time difference of 0.12 secs
    ## 
    ## ================================================================================
    ## Time difference of 39.2006 secs
    ## Adding 4 sequences to the database.
    ## 
    ## Added 4 new sequences to table Seqs.
    ## 20 total sequences in table Seqs.
    ## Time difference of 0.1 secs
    ## 
    ## ================================================================================
    ## Time difference of 28.37736 secs

``` r
names(GC01) <- seq(length(GC01))
```

## Build a Synteny Map

See the DECIPHER documentation for further reading. This is a modestly
diverse set of genomes within the same genera, and what we see in the
below pairs plot, combined with a general summary of syntenic hits is
that though they share reasonable amounts of ordered information, those
individual blocks of order are still pretty fractious. This is a
modestly difficult test case for our tools in this workflow.

``` r
Syn01 <- FindSynteny(dbFile = DBPATH,
                     verbose = TRUE,
                     processors = NULL)
```

    ## ================================================================================
    ## 
    ## Time difference of 32 secs

``` r
print(Syn01)
```

    ##             1           2           3           4        5
    ## 1      7 seqs    13% hits    12% hits    13% hits 20% hits
    ## 2  935 blocks       1 seq    40% hits    36% hits 19% hits
    ## 3  936 blocks  587 blocks       1 seq    41% hits 18% hits
    ## 4 1092 blocks  817 blocks  897 blocks      7 seqs 18% hits
    ## 5 1297 blocks 1309 blocks 1258 blocks 1444 blocks   4 seqs

``` r
pairs(Syn01)
```

![](README_files/figure-gfm/synteny%20map-1.png)<!-- -->

## Build Predicted Pairs

Reconcile syntenic hits with gene calls to build out an initial list of
predicted pairs.

``` r
L01 <- NucleotideOverlap(SyntenyObject = Syn01,
                         GeneCalls = GC01,
                         LimitIndex = FALSE,
                         AcceptContigNames = TRUE,
                         Verbose = TRUE)
```

    ## 
    ## Reconciling genecalls.
    ## ================================================================================
    ## Finding connected features.
    ## ================================================================================
    ## Time difference of 5.673937 secs

``` r
P01 <- PairSummaries(SyntenyLinks = L01,
                     GeneCalls = GC01,
                     DBPATH = DBPATH,
                     PIDs = TRUE,
                     Score = TRUE,
                     Verbose = TRUE)
```

    ## 
    ## Preparing overhead data.
    ## Overhead complete.
    ## Aligning pairs.
    ## ================================================================================
    ## Time difference of 2.984553 mins

``` r
hist(P01$PID, breaks = 100)
```

![](README_files/figure-gfm/pairwise%20connections%20and%20evaluations-1.png)<!-- -->

## State of the Data

Some of these predictions are good, some are bad, and we have somewhat
limited data available to evaluate these pairings. In the `P01` object,
columns 3-8 represent generalized information about the blocks of
information that led to the prediction, while columns 10, 11, and 12
represent general evaluations of the pairings themselves.

``` r
head(P01)
```

    ##       p1       p2 ExactMatch TotalKmers MaxKmer Consensus p1FeatureLength
    ## 1 1_1_10  2_1_298        186          4      69 0.9628983             522
    ## 2 1_1_16 2_1_1575        158          6      54 0.9827541            1107
    ## 3 1_1_25 2_1_2019        306          8      60 0.9666740             912
    ## 4 1_1_28 2_1_1629        163          4      64 0.9702062            1197
    ## 5 1_1_29 2_1_2295        676         12     117 0.9849346            2181
    ## 6 1_1_30 2_1_2840        107          3      50 0.9721277             420
    ##   p2FeatureLength Adjacent    TetDist       PID     SCORE PIDType PredictedPID
    ## 1             474        0 0.06612293 0.5804598  479.1796      AA    0.6482591
    ## 2            1074        0 0.03944655 0.5336927  985.4697      AA    0.5932843
    ## 3             984        0 0.04529876 0.5903614 1028.5381      AA    0.6615063
    ## 4            1128        0 0.04862226 0.4901961  858.3170      AA    0.5415850
    ## 5            2127        0 0.03760356 0.6338798 2189.9736      AA    0.6845313
    ## 6             441        0 0.05721383 0.6394558  516.5079      AA    0.6131397

``` r
dim(P01)
```

    ## [1] 20655    14

``` r
pairs(P01[, c(3:8, 10:12)], pch = 46)
```

![](README_files/figure-gfm/view%20pairs-1.png)<!-- -->

## Evaluate Predicted Pairs

As seen in the above histogram, predicted pairs span the range of
available PIDs. Where users choose to threshold their PIDs is a
complicated series of value judgements based on their downstream needs,
and confidence preferences, however hard cutoffs are often inappropriate
for the data generated herein. To give users at least some tools to
avoid a hard threshold and effectively drop bad predictions, the
`SelectByK` function was constructed to use a K-means approach where
users provide a PID confidence that signifies the **lowest PID centroid
at which a cluster represents TRUE pairings**.

``` r
P02 <- suppressWarnings(SelectByK(Pairs = P01,
                                  ClusterScalar = 4L, # default is 1, though that isn't likely to be aligned with user preferences
                                  UserConfidence = 0.5, # this is the PID that the cluster centroid must be above to be retained
                                  Verbose = TRUE,
                                  ShowPlot = TRUE,
                                  ReturnAllCommunities = TRUE))
```

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10
    ## [1] 11
    ## [1] 12
    ## [1] 13
    ## [1] 14

![](README_files/figure-gfm/naively%20trim%20false%20positives-1.png)<!-- -->

## Evaluated Clusters

We can evaluate the identified clusters with a histogram in relation to
the user cutoff and show the outcome of the selection. The total number
of clusters selected by `SelectByK` is a function of the `ClusterScalar`
argument and the initial input data.

``` r
P03 <- P02[[2]]
P02 <- P02[[1]]

plts <- vector(mode = "list",
               length = length(P03))
maxcounts <- vector(mode = "integer",
                    length = length(P03))
brks <- seq(from = 0,
            to = 1,
            by = 0.01)

for (m1 in seq_along(P03)) {
  plts[[m1]] <- hist(P03[[m1]][, "PID"],
                     breaks = brks,
                     plot = FALSE)
  maxcounts[m1] <- max(plts[[m1]]$counts)
}

o1 <- order(maxcounts, decreasing = TRUE)

plts <- plts[o1]

for (m1 in seq_along(plts)) {
  if (m1 > 1L) {
    plot(plts[[m1]],
         col = ColVec4[m1],
         add = TRUE)
  } else {
    plot(plts[[m1]],
         col = ColVec4[m1],
         xlab = "PID",
         ylab = "Frequency",
         main = "Clustered Predictions")
  }
}

abline(v = 0.5, lwd = 2.5, lty = 2)
```

![](README_files/figure-gfm/clusters%20histogram-1.png)<!-- -->

``` r
rm(P03)
```

## Some Other Things

We can collect the single linkage sets of pairs for further evaluation.

``` r
Sets01 <- DisjointSet(Pairs = P02,
                      Verbose = TRUE)
```

    ## 
    ## Assigning initial root:
    ## ================================================================================
    ## Time difference of 0.103724 secs
    ## 
    ## Assigning final root:
    ##                                                                                 ================================================================================
    ## Time difference of 0.06499505 secs
    ## 
    ## Assigning single linkage clusters.
    ## Assignments complete.
    ## 
    ## Time difference of 0.1826658 secs

## End of workflow

Save off our data and print out some session information.

``` r
TIMEEND <- Sys.time()
# total time for the entire workflow
print(TIMEEND - TIMESTART)
```

    ## Time difference of 6.577236 mins

``` r
save(Syn01,
     L01,
     P01,
     P02,
     Sets01,
     file = paste0("data/",
                   Selection,
                   ".RData"),
     compress = "xz")

sessionInfo()
```

    ## R version 4.3.0 (2023-04-21)
    ## Platform: x86_64-apple-darwin20 (64-bit)
    ## Running under: macOS Monterey 12.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] igraph_1.4.2        SynExtend_1.12.0    DECIPHER_2.28.0    
    ##  [4] RSQLite_2.3.1       Biostrings_2.68.0   GenomeInfoDb_1.36.0
    ##  [7] XVector_0.40.0      IRanges_2.34.0      S4Vectors_0.38.1   
    ## [10] BiocGenerics_0.46.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.6.2             crayon_1.5.2            cli_3.6.1              
    ##  [4] knitr_1.42              rlang_1.1.1             xfun_0.39              
    ##  [7] highr_0.10              DBI_1.1.3               bit_4.0.5              
    ## [10] RCurl_1.98-1.12         htmltools_0.5.5         rmarkdown_2.21         
    ## [13] evaluate_0.21           bitops_1.0-7            fastmap_1.1.1          
    ## [16] yaml_2.3.7              memoise_2.0.1           compiler_4.3.0         
    ## [19] blob_1.2.4              pkgconfig_2.0.3         rstudioapi_0.14        
    ## [22] digest_0.6.31           GenomeInfoDbData_1.2.10 magrittr_2.0.3         
    ## [25] tools_4.3.0             bit64_4.0.5             zlibbioc_1.46.0        
    ## [28] cachem_1.0.8
