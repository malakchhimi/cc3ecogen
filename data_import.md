R Notebook
================

``` bash
# wget -i url_data
```

``` r
path <- "/home/rstudio/cc3ecogen/data"
list.files(path)
```

    ##  [1] "filtered"                "SRR19667519_R1.fastq.gz"
    ##  [3] "SRR19667519_R2.fastq.gz" "SRR19667520_R1.fastq.gz"
    ##  [5] "SRR19667520_R2.fastq.gz" "SRR19667521_R1.fastq.gz"
    ##  [7] "SRR19667521_R2.fastq.gz" "SRR19667522_R1.fastq.gz"
    ##  [9] "SRR19667522_R2.fastq.gz" "SRR19667523_R1.fastq.gz"
    ## [11] "SRR19667523_R2.fastq.gz" "SRR19667524_R1.fastq.gz"
    ## [13] "SRR19667524_R2.fastq.gz" "SRR19667525_R1.fastq.gz"
    ## [15] "SRR19667525_R2.fastq.gz" "SRR19667526_R1.fastq.gz"
    ## [17] "SRR19667526_R2.fastq.gz" "SRR19667527_R1.fastq.gz"
    ## [19] "SRR19667527_R2.fastq.gz" "SRR19667528_R1.fastq.gz"
    ## [21] "SRR19667528_R2.fastq.gz" "SRR19667529_R1.fastq.gz"
    ## [23] "SRR19667529_R2.fastq.gz" "SRR19667530_R1.fastq.gz"
    ## [25] "SRR19667530_R2.fastq.gz" "SRR19667531_R1.fastq.gz"
    ## [27] "SRR19667531_R2.fastq.gz" "SRR19667532_R1.fastq.gz"
    ## [29] "SRR19667532_R2.fastq.gz" "SRR19667533_R1.fastq.gz"
    ## [31] "SRR19667533_R2.fastq.gz" "SRR19667534_R1.fastq.gz"
    ## [33] "SRR19667534_R2.fastq.gz" "SRR19667535_R1.fastq.gz"
    ## [35] "SRR19667535_R2.fastq.gz" "SRR19667536_R1.fastq.gz"
    ## [37] "SRR19667536_R2.fastq.gz" "SRR19667537_R1.fastq.gz"
    ## [39] "SRR19667537_R2.fastq.gz" "SRR19667538_R1.fastq.gz"
    ## [41] "SRR19667538_R2.fastq.gz" "SRR19667539_R1.fastq.gz"
    ## [43] "SRR19667539_R2.fastq.gz" "SRR19667540_R1.fastq.gz"
    ## [45] "SRR19667540_R2.fastq.gz" "SRR19667541_R1.fastq.gz"
    ## [47] "SRR19667541_R2.fastq.gz" "SRR19667542_R1.fastq.gz"
    ## [49] "SRR19667542_R2.fastq.gz" "SRR19667543_R1.fastq.gz"
    ## [51] "SRR19667543_R2.fastq.gz" "SRR19667544_R1.fastq.gz"
    ## [53] "SRR19667544_R2.fastq.gz" "SRR19667545_R1.fastq.gz"
    ## [55] "SRR19667545_R2.fastq.gz" "SRR19667546_R1.fastq.gz"
    ## [57] "SRR19667546_R2.fastq.gz" "SRR19667547_R1.fastq.gz"
    ## [59] "SRR19667547_R2.fastq.gz" "SRR19667548_R1.fastq.gz"
    ## [61] "SRR19667548_R2.fastq.gz" "SRR19667549_R1.fastq.gz"
    ## [63] "SRR19667549_R2.fastq.gz" "SRR19667550_R1.fastq.gz"
    ## [65] "SRR19667550_R2.fastq.gz" "SRR19667551_R1.fastq.gz"
    ## [67] "SRR19667551_R2.fastq.gz" "SRR19667552_R1.fastq.gz"
    ## [69] "SRR19667552_R2.fastq.gz" "SRR19667553_R1.fastq.gz"
    ## [71] "SRR19667553_R2.fastq.gz" "SRR19667554_R1.fastq.gz"
    ## [73] "SRR19667554_R2.fastq.gz" "SRR19667555_R1.fastq.gz"
    ## [75] "SRR19667555_R2.fastq.gz"

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(phyloseq)
library(DECIPHER)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

``` r
library(phangorn)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

``` r
library(ggplot2)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
library(shiny)
library(miniUI)
library(caret)
```

    ## Loading required package: lattice

``` r
library(pls)
```

    ## 
    ## Attaching package: 'pls'

    ## The following object is masked from 'package:caret':
    ## 
    ##     R2

    ## The following object is masked from 'package:ape':
    ## 
    ##     mvr

    ## The following object is masked from 'package:stats':
    ## 
    ##     loadings

``` r
library(e1071)
library(ggplot2)
library(randomForest)
```

    ## randomForest 4.7-1.1

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:randomForest':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggrepel)
#library(nlme)
library(devtools)
```

    ## Loading required package: usethis

``` r
library(reshape2)
library(PMA)
#library(structSSI)
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
library(ggnetwork)
library(intergraph)
library(scales)
library(phyloseqGraphTest)
library(genefilter)
library(impute)
```

``` r
# nous lisons les noms des fichiers et effectuons quelques manipulations de chaînes pour triage selon R1 OU R2,La fonction sort permet d’obtenir le même ordre entre les fastq sens et antisens.
fnFs <- sort(list.files(path, pattern="_R1", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2", full.names = TRUE))
```

``` r
#Etant donné que les paires de fichiers fastq sens et antisens appartiennent au même échantillon, nous allons créer une variable qui extrait le nom de cet échantillon
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
```

``` r
#visualiser les profils de qualité des forward reads 
plotQualityProfile(fnFs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](data_import_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#Dans ces figures, la médianne est en vert, et les quartiles en orange pointille
```

``` r
plotQualityProfile(fnRs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](data_import_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# nous allons créer des objets (filtFs et filtRs) pour stoquer les séquences filtrées.
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

``` r
#450 pb
#filtrer amorces, procédons avec la fonction filterAndTrim, sa sortie sera stocké dans l’objet out.
#truncQ : définit un indice Q score minimale, a la première instance d’un score de qualité inférieur ou égal à truncQ, la séquence est tronquée.
#truncLen : définit à quelle longueur les séquences vont être tronquées,les séquences plus courtes que la longueur sont éliminées.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(18,18), truncLen=c(120,200), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE) 
head(out)
```

    ##                         reads.in reads.out
    ## SRR19667519_R1.fastq.gz   382439    283095
    ## SRR19667520_R1.fastq.gz   247892    188360
    ## SRR19667521_R1.fastq.gz   203290    173731
    ## SRR19667522_R1.fastq.gz   343027    297510
    ## SRR19667523_R1.fastq.gz   282482    247892
    ## SRR19667524_R1.fastq.gz   249004    219441

``` r
#estimer le taux d’erreur de séquençage, son but est de permettre de différencier les séquences mutantes et les séquences érronées
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 121439976 total bases in 1190588 reads from 5 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 117423852 total bases in 645186 reads from 3 samples will be used for learning the error rates.

``` r
#Visualisation du taux d’erreur
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](data_import_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
plotErrors(errR, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](data_import_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#Exemple d'inférence, appliquer l'algorithme d'inférence de l'échantillon de base aux données de séquence filtrées et découpées
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 283095 reads in 70626 unique sequences.
    ## Sample 2 - 188360 reads in 51485 unique sequences.
    ## Sample 3 - 173731 reads in 29374 unique sequences.
    ## Sample 4 - 297510 reads in 42953 unique sequences.
    ## Sample 5 - 247892 reads in 43105 unique sequences.
    ## Sample 6 - 219441 reads in 45284 unique sequences.
    ## Sample 7 - 158382 reads in 29337 unique sequences.
    ## Sample 8 - 139167 reads in 28404 unique sequences.
    ## Sample 9 - 176730 reads in 29693 unique sequences.
    ## Sample 10 - 516035 reads in 77732 unique sequences.
    ## Sample 11 - 74076 reads in 15632 unique sequences.
    ## Sample 12 - 59783 reads in 12732 unique sequences.
    ## Sample 13 - 799800 reads in 125929 unique sequences.
    ## Sample 14 - 196090 reads in 32951 unique sequences.
    ## Sample 15 - 394220 reads in 75895 unique sequences.
    ## Sample 16 - 143153 reads in 33162 unique sequences.
    ## Sample 17 - 177403 reads in 41161 unique sequences.
    ## Sample 18 - 110704 reads in 31563 unique sequences.
    ## Sample 19 - 205830 reads in 44708 unique sequences.
    ## Sample 20 - 150117 reads in 36002 unique sequences.
    ## Sample 21 - 125913 reads in 31842 unique sequences.
    ## Sample 22 - 228692 reads in 50559 unique sequences.
    ## Sample 23 - 180879 reads in 33749 unique sequences.
    ## Sample 24 - 98141 reads in 16944 unique sequences.
    ## Sample 25 - 118905 reads in 28353 unique sequences.
    ## Sample 26 - 1447823 reads in 232037 unique sequences.
    ## Sample 27 - 234576 reads in 49089 unique sequences.
    ## Sample 28 - 128496 reads in 31213 unique sequences.
    ## Sample 29 - 247015 reads in 61612 unique sequences.
    ## Sample 30 - 257913 reads in 54729 unique sequences.
    ## Sample 31 - 380831 reads in 67248 unique sequences.
    ## Sample 32 - 144681 reads in 43870 unique sequences.
    ## Sample 33 - 125299 reads in 31192 unique sequences.
    ## Sample 34 - 307071 reads in 76908 unique sequences.
    ## Sample 35 - 137242 reads in 39649 unique sequences.
    ## Sample 36 - 193888 reads in 32125 unique sequences.
    ## Sample 37 - 207470 reads in 44325 unique sequences.

  

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 283095 reads in 147121 unique sequences.
    ## Sample 2 - 188360 reads in 101407 unique sequences.
    ## Sample 3 - 173731 reads in 72337 unique sequences.
    ## Sample 4 - 297510 reads in 118886 unique sequences.
    ## Sample 5 - 247892 reads in 108742 unique sequences.
    ## Sample 6 - 219441 reads in 91593 unique sequences.
    ## Sample 7 - 158382 reads in 60800 unique sequences.
    ## Sample 8 - 139167 reads in 56033 unique sequences.
    ## Sample 9 - 176730 reads in 73547 unique sequences.
    ## Sample 10 - 516035 reads in 172383 unique sequences.
    ## Sample 11 - 74076 reads in 33965 unique sequences.
    ## Sample 12 - 59783 reads in 27497 unique sequences.
    ## Sample 13 - 799800 reads in 304817 unique sequences.
    ## Sample 14 - 196090 reads in 73345 unique sequences.
    ## Sample 15 - 394220 reads in 188212 unique sequences.
    ## Sample 16 - 143153 reads in 77314 unique sequences.
    ## Sample 17 - 177403 reads in 96075 unique sequences.
    ## Sample 18 - 110704 reads in 71691 unique sequences.
    ## Sample 19 - 205830 reads in 102849 unique sequences.
    ## Sample 20 - 150117 reads in 75675 unique sequences.
    ## Sample 21 - 125913 reads in 70812 unique sequences.
    ## Sample 22 - 228692 reads in 105140 unique sequences.
    ## Sample 23 - 180879 reads in 74736 unique sequences.
    ## Sample 24 - 98141 reads in 37923 unique sequences.
    ## Sample 25 - 118905 reads in 63585 unique sequences.
    ## Sample 26 - 1447823 reads in 754019 unique sequences.
    ## Sample 27 - 234576 reads in 151961 unique sequences.
    ## Sample 28 - 128496 reads in 102568 unique sequences.
    ## Sample 29 - 247015 reads in 158697 unique sequences.
    ## Sample 30 - 257913 reads in 163491 unique sequences.
    ## Sample 31 - 380831 reads in 263140 unique sequences.
    ## Sample 32 - 144681 reads in 85210 unique sequences.
    ## Sample 33 - 125299 reads in 66387 unique sequences.
    ## Sample 34 - 307071 reads in 163323 unique sequences.
    ## Sample 35 - 137242 reads in 76111 unique sequences.
    ## Sample 36 - 193888 reads in 82058 unique sequences.
    ## Sample 37 - 207470 reads in 93135 unique sequences.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1531 sequence variants were inferred from 70626 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
#L'algorithme DADA2 a déduit 1531 vrais variants de séquence à partir des 70626 séquences uniques du premier échantillon
```

``` r
dadaRs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1541 sequence variants were inferred from 147121 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
##L'algorithme DADA2 a déduit 1541 vrais variants de séquence à partir des 147121 séquences uniques du premier échantillon
```

``` r
#Tout l’intérêt du séquençage en paire réside dans le but de fusionner les deux brins afin d’accroitre notre confiance en leur fiabilité. La fusion permet également d’obtenir des amplicons plus long.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 4279 paired-reads (in 559 unique pairings) successfully merged out of 275673 (in 72457 pairings) input.

    ## 3066 paired-reads (in 413 unique pairings) successfully merged out of 181864 (in 52319 pairings) input.

    ## 1610 paired-reads (in 230 unique pairings) successfully merged out of 168551 (in 45155 pairings) input.

    ## 1304 paired-reads (in 275 unique pairings) successfully merged out of 290352 (in 61809 pairings) input.

    ## 1020 paired-reads (in 293 unique pairings) successfully merged out of 239949 (in 74257 pairings) input.

    ## 315 paired-reads (in 73 unique pairings) successfully merged out of 211213 (in 57047 pairings) input.

    ## 545 paired-reads (in 140 unique pairings) successfully merged out of 154456 (in 27609 pairings) input.

    ## 250 paired-reads (in 77 unique pairings) successfully merged out of 133935 (in 32691 pairings) input.

    ## 439 paired-reads (in 110 unique pairings) successfully merged out of 171952 (in 31450 pairings) input.

    ## 1406 paired-reads (in 296 unique pairings) successfully merged out of 505293 (in 85536 pairings) input.

    ## 71 paired-reads (in 26 unique pairings) successfully merged out of 71520 (in 16076 pairings) input.

    ## 144 paired-reads (in 37 unique pairings) successfully merged out of 57584 (in 14019 pairings) input.

    ## 2163 paired-reads (in 473 unique pairings) successfully merged out of 783560 (in 167036 pairings) input.

    ## 50 paired-reads (in 27 unique pairings) successfully merged out of 190480 (in 37335 pairings) input.

    ## 1539 paired-reads (in 348 unique pairings) successfully merged out of 381124 (in 121196 pairings) input.

    ## 1067 paired-reads (in 181 unique pairings) successfully merged out of 134817 (in 59081 pairings) input.

    ## 1011 paired-reads (in 211 unique pairings) successfully merged out of 169545 (in 71832 pairings) input.

    ## 1878 paired-reads (in 299 unique pairings) successfully merged out of 102792 (in 52930 pairings) input.

    ## 1372 paired-reads (in 318 unique pairings) successfully merged out of 196151 (in 79306 pairings) input.

    ## 2706 paired-reads (in 417 unique pairings) successfully merged out of 143486 (in 66681 pairings) input.

    ## 2494 paired-reads (in 413 unique pairings) successfully merged out of 118883 (in 49717 pairings) input.

    ## 4332 paired-reads (in 695 unique pairings) successfully merged out of 220019 (in 75840 pairings) input.

    ## 373 paired-reads (in 84 unique pairings) successfully merged out of 175779 (in 38791 pairings) input.

    ## 149 paired-reads (in 53 unique pairings) successfully merged out of 95124 (in 16651 pairings) input.

    ## 214 paired-reads (in 40 unique pairings) successfully merged out of 111856 (in 42684 pairings) input.

    ## 9141 paired-reads (in 1368 unique pairings) successfully merged out of 1425059 (in 344596 pairings) input.

    ## 12792 paired-reads (in 290 unique pairings) successfully merged out of 227655 (in 62827 pairings) input.

    ## 934 paired-reads (in 162 unique pairings) successfully merged out of 121756 (in 44928 pairings) input.

    ## 2013 paired-reads (in 334 unique pairings) successfully merged out of 238115 (in 80698 pairings) input.

    ## 1397 paired-reads (in 258 unique pairings) successfully merged out of 250058 (in 70920 pairings) input.

    ## 2695 paired-reads (in 403 unique pairings) successfully merged out of 370553 (in 104880 pairings) input.

    ## 10895 paired-reads (in 658 unique pairings) successfully merged out of 140062 (in 40370 pairings) input.

    ## 2345 paired-reads (in 256 unique pairings) successfully merged out of 121584 (in 28545 pairings) input.

    ## 4524 paired-reads (in 578 unique pairings) successfully merged out of 299072 (in 79187 pairings) input.

    ## 1846 paired-reads (in 302 unique pairings) successfully merged out of 132075 (in 40082 pairings) input.

    ## 152 paired-reads (in 28 unique pairings) successfully merged out of 187412 (in 45812 pairings) input.

    ## 372 paired-reads (in 54 unique pairings) successfully merged out of 197914 (in 65454 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                   sequence
    ## 109 CGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTTCCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGATAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTATCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGGA
    ## 112 CGCGGTAATACAGAGGGTGCAAGCGTTAATCGGATTTACTGGGCGTAAAGCGCGCGTAGGTGGCCAATTAAGTCAAATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGTTGGCTAGAGTATGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAATACTGACACTGAGGTGCGAAAGCATGGGGAGCAAACAGGA
    ## 120                                                                                 TGTATAAGAGACAGCCTACGGGTGGCTGCATTCCAAACTGGAAGTCTAGAGTGCAGGAGAGGAGAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGAGATTAGGAAGAACACCAGTGGCGAAGGCGACTCTCTGGACTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGA
    ## 133 CGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGA
    ## 190  CGCGGTAATACGTAGGGGGCTAGCGTTATCCGGAATTACTGGGCGTAAAGGGTGCGTAGGTGGTTTCTTAAGTCAGAGGTGAAAGGCTACGGCTCAACCGTAGTAAGCCTTTGAAACTGGGAAACTTGAGTGCAGGAGAGGAGAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTTGCGAAGGCGGCTCTCTGGACTGTAACTGACACTGAGGCACGAAAGCGTGGGGAGCAAACAGGA
    ## 232                                                                                 TGTATAAGAGACAGCCTACGGGTGGCTGCATCTGATACTGGTTGACTTGAGTGCAGGAGAGGAAAGCGGAATTCCTAGTGTAGCGGTGAAATGCATAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTTCTGGACTGCAACTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGGA
    ##     abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 109       173     190       3     22         0      0      2   TRUE
    ## 112       168     201       2     22         0      0      2   TRUE
    ## 120       162     315     805    102         0      0      1   TRUE
    ## 133       150     233       1     22         0      0      2   TRUE
    ## 190       121     177       4     23         0      0      2   TRUE
    ## 232       105     389     933    102         0      0      1   TRUE

``` r
#Nous avons enfin nos ASVs que nous allons stoquer dans l’objet seqtab, construire une table de variantes de séquence d'amplicon (ASV), une version à plus haute résolution de la table OTU produite par des méthodes traditionnelles.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]   37 6258

``` r
# Ce tableau contient 376258 ASV
```

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  182  183  184  185  187  188  189  190  191  192  193  194  195  196  197  198 
    ##  255    1    1    2    1    7    3   36   51   26   46    1   40    2    3    3 
    ##  199  200  201  202  204  206  207  209  210  211  213  215  216  217  222  223 
    ##    4   13  108    3    2    1    1    2    1    2    2    1    2    3    1    1 
    ##  229  230  231  232  233  234  235  236  237  238  239  240  241  244  245  247 
    ##    2    3   17   37    2    3    5    6    3   45   11    2    1    2    3    2 
    ##  248  249  254  255  257  258  259  260  261  262  263  264  265  266  268  272 
    ##    5    2    3    2    1    2    7   12  951 4179  106  198   15    2    2    2

``` r
#Suppression des chimères,cette étape vise à éliminer toutes les séquences non biologiques
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 2514 bimeras out of 6258 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   37 3744

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9065052

``` r
#Tableau de suivi,voir combien de séquences ont été éliminées à chaque étape
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##              input filtered denoisedF denoisedR merged nonchim
    ## SRR19667519 382439   283095    280494    278091   4279    3835
    ## SRR19667520 247892   188360    186150    183912   3066    2777
    ## SRR19667521 203290   173731    171993    170158   1610    1384
    ## SRR19667522 343027   297510    295327    292382   1304    1053
    ## SRR19667523 282482   247892    245339    242372   1020     782
    ## SRR19667524 249004   219441    217073    213368    315     295

``` bash
#attribuer une taxonomie aux variants de séquence
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1
```

    ## --2023-01-07 19:11:55--  https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1
    ## Resolving zenodo.org (zenodo.org)... 188.185.124.72
    ## Connecting to zenodo.org (zenodo.org)|188.185.124.72|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137283333 (131M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138.1_train_set.fa.gz?download=1.3’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 6.09M 21s
    ##     50K .......... .......... .......... .......... ..........  0% 14.5M 15s
    ##    100K .......... .......... .......... .......... ..........  0% 13.6M 13s
    ##    150K .......... .......... .......... .......... ..........  0% 14.4M 12s
    ##    200K .......... .......... .......... .......... ..........  0% 19.1M 11s
    ##    250K .......... .......... .......... .......... ..........  0% 51.0M 10s
    ##    300K .......... .......... .......... .......... ..........  0% 18.2M 9s
    ##    350K .......... .......... .......... .......... ..........  0% 52.5M 9s
    ##    400K .......... .......... .......... .......... ..........  0% 19.7M 8s
    ##    450K .......... .......... .......... .......... ..........  0% 17.9M 8s
    ##    500K .......... .......... .......... .......... ..........  0% 54.1M 8s
    ##    550K .......... .......... .......... .......... ..........  0% 60.3M 7s
    ##    600K .......... .......... .......... .......... ..........  0% 20.4M 7s
    ##    650K .......... .......... .......... .......... ..........  0% 19.4M 7s
    ##    700K .......... .......... .......... .......... ..........  0% 48.5M 7s
    ##    750K .......... .......... .......... .......... ..........  0% 24.6M 7s
    ##    800K .......... .......... .......... .......... ..........  0% 35.1M 7s
    ##    850K .......... .......... .......... .......... ..........  0% 18.5M 7s
    ##    900K .......... .......... .......... .......... ..........  0% 59.8M 6s
    ##    950K .......... .......... .......... .......... ..........  0% 23.0M 6s
    ##   1000K .......... .......... .......... .......... ..........  0% 42.7M 6s
    ##   1050K .......... .......... .......... .......... ..........  0% 19.6M 6s
    ##   1100K .......... .......... .......... .......... ..........  0% 59.2M 6s
    ##   1150K .......... .......... .......... .......... ..........  0% 18.4M 6s
    ##   1200K .......... .......... .......... .......... ..........  0% 63.8M 6s
    ##   1250K .......... .......... .......... .......... ..........  0% 21.7M 6s
    ##   1300K .......... .......... .......... .......... ..........  1% 34.4M 6s
    ##   1350K .......... .......... .......... .......... ..........  1% 21.2M 6s
    ##   1400K .......... .......... .......... .......... ..........  1% 59.8M 6s
    ##   1450K .......... .......... .......... .......... ..........  1% 41.7M 6s
    ##   1500K .......... .......... .......... .......... ..........  1% 20.5M 6s
    ##   1550K .......... .......... .......... .......... ..........  1% 34.2M 6s
    ##   1600K .......... .......... .......... .......... ..........  1% 27.4M 6s
    ##   1650K .......... .......... .......... .......... ..........  1% 36.2M 5s
    ##   1700K .......... .......... .......... .......... ..........  1% 45.5M 5s
    ##   1750K .......... .......... .......... .......... ..........  1% 25.2M 5s
    ##   1800K .......... .......... .......... .......... ..........  1% 25.1M 5s
    ##   1850K .......... .......... .......... .......... ..........  1% 67.1M 5s
    ##   1900K .......... .......... .......... .......... ..........  1% 50.2M 5s
    ##   1950K .......... .......... .......... .......... ..........  1% 19.8M 5s
    ##   2000K .......... .......... .......... .......... ..........  1% 30.9M 5s
    ##   2050K .......... .......... .......... .......... ..........  1% 43.2M 5s
    ##   2100K .......... .......... .......... .......... ..........  1% 36.6M 5s
    ##   2150K .......... .......... .......... .......... ..........  1% 33.9M 5s
    ##   2200K .......... .......... .......... .......... ..........  1% 24.5M 5s
    ##   2250K .......... .......... .......... .......... ..........  1% 52.0M 5s
    ##   2300K .......... .......... .......... .......... ..........  1% 26.2M 5s
    ##   2350K .......... .......... .......... .......... ..........  1% 24.9M 5s
    ##   2400K .......... .......... .......... .......... ..........  1% 61.7M 5s
    ##   2450K .......... .......... .......... .......... ..........  1% 20.8M 5s
    ##   2500K .......... .......... .......... .......... ..........  1% 38.5M 5s
    ##   2550K .......... .......... .......... .......... ..........  1% 37.5M 5s
    ##   2600K .......... .......... .......... .......... ..........  1% 47.9M 5s
    ##   2650K .......... .......... .......... .......... ..........  2% 24.8M 5s
    ##   2700K .......... .......... .......... .......... ..........  2% 39.7M 5s
    ##   2750K .......... .......... .......... .......... ..........  2% 19.6M 5s
    ##   2800K .......... .......... .......... .......... ..........  2% 45.6M 5s
    ##   2850K .......... .......... .......... .......... ..........  2% 61.7M 5s
    ##   2900K .......... .......... .......... .......... ..........  2% 26.6M 5s
    ##   2950K .......... .......... .......... .......... ..........  2% 26.3M 5s
    ##   3000K .......... .......... .......... .......... ..........  2% 38.0M 5s
    ##   3050K .......... .......... .......... .......... ..........  2% 68.5M 5s
    ##   3100K .......... .......... .......... .......... ..........  2% 29.6M 5s
    ##   3150K .......... .......... .......... .......... ..........  2% 27.0M 5s
    ##   3200K .......... .......... .......... .......... ..........  2% 30.5M 5s
    ##   3250K .......... .......... .......... .......... ..........  2% 45.0M 5s
    ##   3300K .......... .......... .......... .......... ..........  2% 40.9M 5s
    ##   3350K .......... .......... .......... .......... ..........  2% 28.6M 5s
    ##   3400K .......... .......... .......... .......... ..........  2% 42.5M 5s
    ##   3450K .......... .......... .......... .......... ..........  2% 28.6M 5s
    ##   3500K .......... .......... .......... .......... ..........  2% 63.0M 5s
    ##   3550K .......... .......... .......... .......... ..........  2% 29.8M 5s
    ##   3600K .......... .......... .......... .......... ..........  2% 32.0M 5s
    ##   3650K .......... .......... .......... .......... ..........  2% 37.8M 5s
    ##   3700K .......... .......... .......... .......... ..........  2% 38.8M 5s
    ##   3750K .......... .......... .......... .......... ..........  2% 57.5M 4s
    ##   3800K .......... .......... .......... .......... ..........  2% 45.2M 4s
    ##   3850K .......... .......... .......... .......... ..........  2% 36.6M 4s
    ##   3900K .......... .......... .......... .......... ..........  2% 37.2M 4s
    ##   3950K .......... .......... .......... .......... ..........  2% 33.0M 4s
    ##   4000K .......... .......... .......... .......... ..........  3% 48.3M 4s
    ##   4050K .......... .......... .......... .......... ..........  3% 50.0M 4s
    ##   4100K .......... .......... .......... .......... ..........  3% 37.7M 4s
    ##   4150K .......... .......... .......... .......... ..........  3% 31.8M 4s
    ##   4200K .......... .......... .......... .......... ..........  3% 66.2M 4s
    ##   4250K .......... .......... .......... .......... ..........  3% 33.3M 4s
    ##   4300K .......... .......... .......... .......... ..........  3% 43.6M 4s
    ##   4350K .......... .......... .......... .......... ..........  3% 33.6M 4s
    ##   4400K .......... .......... .......... .......... ..........  3% 31.7M 4s
    ##   4450K .......... .......... .......... .......... ..........  3% 37.9M 4s
    ##   4500K .......... .......... .......... .......... ..........  3% 50.2M 4s
    ##   4550K .......... .......... .......... .......... ..........  3% 40.3M 4s
    ##   4600K .......... .......... .......... .......... ..........  3% 34.6M 4s
    ##   4650K .......... .......... .......... .......... ..........  3% 61.7M 4s
    ##   4700K .......... .......... .......... .......... ..........  3% 36.8M 4s
    ##   4750K .......... .......... .......... .......... ..........  3% 40.4M 4s
    ##   4800K .......... .......... .......... .......... ..........  3% 36.0M 4s
    ##   4850K .......... .......... .......... .......... ..........  3% 31.3M 4s
    ##   4900K .......... .......... .......... .......... ..........  3% 59.2M 4s
    ##   4950K .......... .......... .......... .......... ..........  3% 46.6M 4s
    ##   5000K .......... .......... .......... .......... ..........  3% 32.1M 4s
    ##   5050K .......... .......... .......... .......... ..........  3% 41.4M 4s
    ##   5100K .......... .......... .......... .......... ..........  3% 37.6M 4s
    ##   5150K .......... .......... .......... .......... ..........  3% 40.8M 4s
    ##   5200K .......... .......... .......... .......... ..........  3% 70.6M 4s
    ##   5250K .......... .......... .......... .......... ..........  3% 32.9M 4s
    ##   5300K .......... .......... .......... .......... ..........  3% 42.9M 4s
    ##   5350K .......... .......... .......... .......... ..........  4% 34.5M 4s
    ##   5400K .......... .......... .......... .......... ..........  4% 62.8M 4s
    ##   5450K .......... .......... .......... .......... ..........  4% 42.2M 4s
    ##   5500K .......... .......... .......... .......... ..........  4% 35.9M 4s
    ##   5550K .......... .......... .......... .......... ..........  4% 33.8M 4s
    ##   5600K .......... .......... .......... .......... ..........  4% 42.0M 4s
    ##   5650K .......... .......... .......... .......... ..........  4% 51.9M 4s
    ##   5700K .......... .......... .......... .......... ..........  4% 37.4M 4s
    ##   5750K .......... .......... .......... .......... ..........  4% 31.4M 4s
    ##   5800K .......... .......... .......... .......... ..........  4% 39.2M 4s
    ##   5850K .......... .......... .......... .......... ..........  4% 41.6M 4s
    ##   5900K .......... .......... .......... .......... ..........  4% 62.2M 4s
    ##   5950K .......... .......... .......... .......... ..........  4% 30.9M 4s
    ##   6000K .......... .......... .......... .......... ..........  4% 39.4M 4s
    ##   6050K .......... .......... .......... .......... ..........  4% 50.2M 4s
    ##   6100K .......... .......... .......... .......... ..........  4% 28.8M 4s
    ##   6150K .......... .......... .......... .......... ..........  4% 53.9M 4s
    ##   6200K .......... .......... .......... .......... ..........  4% 33.8M 4s
    ##   6250K .......... .......... .......... .......... ..........  4% 44.0M 4s
    ##   6300K .......... .......... .......... .......... ..........  4% 60.5M 4s
    ##   6350K .......... .......... .......... .......... ..........  4% 34.4M 4s
    ##   6400K .......... .......... .......... .......... ..........  4% 24.2M 4s
    ##   6450K .......... .......... .......... .......... ..........  4% 52.1M 4s
    ##   6500K .......... .......... .......... .......... ..........  4% 33.9M 4s
    ##   6550K .......... .......... .......... .......... ..........  4% 50.0M 4s
    ##   6600K .......... .......... .......... .......... ..........  4% 60.8M 4s
    ##   6650K .......... .......... .......... .......... ..........  4% 27.7M 4s
    ##   6700K .......... .......... .......... .......... ..........  5% 53.4M 4s
    ##   6750K .......... .......... .......... .......... ..........  5% 34.4M 4s
    ##   6800K .......... .......... .......... .......... ..........  5% 51.4M 4s
    ##   6850K .......... .......... .......... .......... ..........  5% 52.2M 4s
    ##   6900K .......... .......... .......... .......... ..........  5% 31.9M 4s
    ##   6950K .......... .......... .......... .......... ..........  5% 55.6M 4s
    ##   7000K .......... .......... .......... .......... ..........  5% 33.0M 4s
    ##   7050K .......... .......... .......... .......... ..........  5% 41.4M 4s
    ##   7100K .......... .......... .......... .......... ..........  5% 68.5M 4s
    ##   7150K .......... .......... .......... .......... ..........  5% 24.2M 4s
    ##   7200K .......... .......... .......... .......... ..........  5% 60.5M 4s
    ##   7250K .......... .......... .......... .......... ..........  5% 36.8M 4s
    ##   7300K .......... .......... .......... .......... ..........  5% 54.4M 4s
    ##   7350K .......... .......... .......... .......... ..........  5% 48.6M 4s
    ##   7400K .......... .......... .......... .......... ..........  5% 27.1M 4s
    ##   7450K .......... .......... .......... .......... ..........  5% 56.9M 4s
    ##   7500K .......... .......... .......... .......... ..........  5% 58.9M 4s
    ##   7550K .......... .......... .......... .......... ..........  5% 27.7M 4s
    ##   7600K .......... .......... .......... .......... ..........  5% 31.2M 4s
    ##   7650K .......... .......... .......... .......... ..........  5% 53.7M 4s
    ##   7700K .......... .......... .......... .......... ..........  5% 67.4M 4s
    ##   7750K .......... .......... .......... .......... ..........  5% 45.3M 4s
    ##   7800K .......... .......... .......... .......... ..........  5% 23.7M 4s
    ##   7850K .......... .......... .......... .......... ..........  5% 60.2M 4s
    ##   7900K .......... .......... .......... .......... ..........  5% 44.6M 4s
    ##   7950K .......... .......... .......... .......... ..........  5% 33.8M 4s
    ##   8000K .......... .......... .......... .......... ..........  6% 56.3M 4s
    ##   8050K .......... .......... .......... .......... ..........  6% 35.3M 4s
    ##   8100K .......... .......... .......... .......... ..........  6% 30.5M 4s
    ##   8150K .......... .......... .......... .......... ..........  6% 51.0M 4s
    ##   8200K .......... .......... .......... .......... ..........  6% 73.4M 4s
    ##   8250K .......... .......... .......... .......... ..........  6% 29.7M 4s
    ##   8300K .......... .......... .......... .......... ..........  6% 40.8M 4s
    ##   8350K .......... .......... .......... .......... ..........  6% 62.3M 4s
    ##   8400K .......... .......... .......... .......... ..........  6% 25.9M 4s
    ##   8450K .......... .......... .......... .......... ..........  6% 54.7M 4s
    ##   8500K .......... .......... .......... .......... ..........  6% 64.7M 4s
    ##   8550K .......... .......... .......... .......... ..........  6% 41.9M 4s
    ##   8600K .......... .......... .......... .......... ..........  6% 23.1M 4s
    ##   8650K .......... .......... .......... .......... ..........  6% 52.5M 4s
    ##   8700K .......... .......... .......... .......... ..........  6% 63.8M 4s
    ##   8750K .......... .......... .......... .......... ..........  6% 28.5M 4s
    ##   8800K .......... .......... .......... .......... ..........  6% 39.9M 4s
    ##   8850K .......... .......... .......... .......... ..........  6% 50.3M 4s
    ##   8900K .......... .......... .......... .......... ..........  6% 47.2M 4s
    ##   8950K .......... .......... .......... .......... ..........  6% 28.8M 4s
    ##   9000K .......... .......... .......... .......... ..........  6% 42.2M 4s
    ##   9050K .......... .......... .......... .......... ..........  6% 70.4M 4s
    ##   9100K .......... .......... .......... .......... ..........  6% 42.7M 4s
    ##   9150K .......... .......... .......... .......... ..........  6% 27.5M 4s
    ##   9200K .......... .......... .......... .......... ..........  6% 59.1M 4s
    ##   9250K .......... .......... .......... .......... ..........  6% 32.9M 4s
    ##   9300K .......... .......... .......... .......... ..........  6% 35.3M 4s
    ##   9350K .......... .......... .......... .......... ..........  7% 46.2M 4s
    ##   9400K .......... .......... .......... .......... ..........  7% 59.4M 4s
    ##   9450K .......... .......... .......... .......... ..........  7% 28.8M 4s
    ##   9500K .......... .......... .......... .......... ..........  7% 49.6M 4s
    ##   9550K .......... .......... .......... .......... ..........  7% 52.5M 4s
    ##   9600K .......... .......... .......... .......... ..........  7% 38.8M 4s
    ##   9650K .......... .......... .......... .......... ..........  7% 35.9M 4s
    ##   9700K .......... .......... .......... .......... ..........  7% 33.0M 4s
    ##   9750K .......... .......... .......... .......... ..........  7% 54.4M 4s
    ##   9800K .......... .......... .......... .......... ..........  7% 29.4M 4s
    ##   9850K .......... .......... .......... .......... ..........  7% 57.1M 3s
    ##   9900K .......... .......... .......... .......... ..........  7% 50.6M 3s
    ##   9950K .......... .......... .......... .......... ..........  7% 43.1M 3s
    ##  10000K .......... .......... .......... .......... ..........  7% 32.8M 3s
    ##  10050K .......... .......... .......... .......... ..........  7% 57.4M 3s
    ##  10100K .......... .......... .......... .......... ..........  7% 71.2M 3s
    ##  10150K .......... .......... .......... .......... ..........  7% 23.0M 3s
    ##  10200K .......... .......... .......... .......... ..........  7% 44.6M 3s
    ##  10250K .......... .......... .......... .......... ..........  7% 34.8M 3s
    ##  10300K .......... .......... .......... .......... ..........  7% 25.6M 3s
    ##  10350K .......... .......... .......... .......... ..........  7% 67.3M 3s
    ##  10400K .......... .......... .......... .......... ..........  7% 57.9M 3s
    ##  10450K .......... .......... .......... .......... ..........  7% 61.7M 3s
    ##  10500K .......... .......... .......... .......... ..........  7% 32.5M 3s
    ##  10550K .......... .......... .......... .......... ..........  7% 65.3M 3s
    ##  10600K .......... .......... .......... .......... ..........  7% 43.6M 3s
    ##  10650K .......... .......... .......... .......... ..........  7% 24.8M 3s
    ##  10700K .......... .......... .......... .......... ..........  8% 43.6M 3s
    ##  10750K .......... .......... .......... .......... ..........  8% 28.7M 3s
    ##  10800K .......... .......... .......... .......... ..........  8% 63.1M 3s
    ##  10850K .......... .......... .......... .......... ..........  8% 58.6M 3s
    ##  10900K .......... .......... .......... .......... ..........  8% 31.3M 3s
    ##  10950K .......... .......... .......... .......... ..........  8% 55.6M 3s
    ##  11000K .......... .......... .......... .......... ..........  8% 31.0M 3s
    ##  11050K .......... .......... .......... .......... ..........  8% 54.4M 3s
    ##  11100K .......... .......... .......... .......... ..........  8% 41.3M 3s
    ##  11150K .......... .......... .......... .......... ..........  8% 30.2M 3s
    ##  11200K .......... .......... .......... .......... ..........  8% 53.9M 3s
    ##  11250K .......... .......... .......... .......... ..........  8% 33.9M 3s
    ##  11300K .......... .......... .......... .......... ..........  8% 57.6M 3s
    ##  11350K .......... .......... .......... .......... ..........  8% 46.9M 3s
    ##  11400K .......... .......... .......... .......... ..........  8% 31.3M 3s
    ##  11450K .......... .......... .......... .......... ..........  8% 49.4M 3s
    ##  11500K .......... .......... .......... .......... ..........  8% 51.7M 3s
    ##  11550K .......... .......... .......... .......... ..........  8% 27.2M 3s
    ##  11600K .......... .......... .......... .......... ..........  8% 63.5M 3s
    ##  11650K .......... .......... .......... .......... ..........  8% 28.3M 3s
    ##  11700K .......... .......... .......... .......... ..........  8% 60.2M 3s
    ##  11750K .......... .......... .......... .......... ..........  8% 31.0M 3s
    ##  11800K .......... .......... .......... .......... ..........  8% 45.8M 3s
    ##  11850K .......... .......... .......... .......... ..........  8% 35.6M 3s
    ##  11900K .......... .......... .......... .......... ..........  8% 54.9M 3s
    ##  11950K .......... .......... .......... .......... ..........  8% 24.6M 3s
    ##  12000K .......... .......... .......... .......... ..........  8% 45.9M 3s
    ##  12050K .......... .......... .......... .......... ..........  9% 72.1M 3s
    ##  12100K .......... .......... .......... .......... ..........  9% 33.7M 3s
    ##  12150K .......... .......... .......... .......... ..........  9% 58.9M 3s
    ##  12200K .......... .......... .......... .......... ..........  9% 42.8M 3s
    ##  12250K .......... .......... .......... .......... ..........  9% 31.8M 3s
    ##  12300K .......... .......... .......... .......... ..........  9% 35.3M 3s
    ##  12350K .......... .......... .......... .......... ..........  9% 35.7M 3s
    ##  12400K .......... .......... .......... .......... ..........  9% 45.9M 3s
    ##  12450K .......... .......... .......... .......... ..........  9% 38.3M 3s
    ##  12500K .......... .......... .......... .......... ..........  9% 23.7M 3s
    ##  12550K .......... .......... .......... .......... ..........  9% 44.0M 3s
    ##  12600K .......... .......... .......... .......... ..........  9% 60.7M 3s
    ##  12650K .......... .......... .......... .......... ..........  9% 40.8M 3s
    ##  12700K .......... .......... .......... .......... ..........  9% 46.0M 3s
    ##  12750K .......... .......... .......... .......... ..........  9% 41.7M 3s
    ##  12800K .......... .......... .......... .......... ..........  9% 47.1M 3s
    ##  12850K .......... .......... .......... .......... ..........  9% 42.7M 3s
    ##  12900K .......... .......... .......... .......... ..........  9% 22.8M 3s
    ##  12950K .......... .......... .......... .......... ..........  9% 29.5M 3s
    ##  13000K .......... .......... .......... .......... ..........  9% 59.5M 3s
    ##  13050K .......... .......... .......... .......... ..........  9% 45.6M 3s
    ##  13100K .......... .......... .......... .......... ..........  9% 37.5M 3s
    ##  13150K .......... .......... .......... .......... ..........  9% 36.8M 3s
    ##  13200K .......... .......... .......... .......... ..........  9% 55.2M 3s
    ##  13250K .......... .......... .......... .......... ..........  9% 44.8M 3s
    ##  13300K .......... .......... .......... .......... ..........  9% 42.7M 3s
    ##  13350K .......... .......... .......... .......... ..........  9% 43.9M 3s
    ##  13400K .......... .......... .......... .......... .......... 10% 32.5M 3s
    ##  13450K .......... .......... .......... .......... .......... 10% 48.5M 3s
    ##  13500K .......... .......... .......... .......... .......... 10% 41.4M 3s
    ##  13550K .......... .......... .......... .......... .......... 10% 27.2M 3s
    ##  13600K .......... .......... .......... .......... .......... 10% 50.2M 3s
    ##  13650K .......... .......... .......... .......... .......... 10% 66.7M 3s
    ##  13700K .......... .......... .......... .......... .......... 10% 33.5M 3s
    ##  13750K .......... .......... .......... .......... .......... 10% 33.7M 3s
    ##  13800K .......... .......... .......... .......... .......... 10% 36.6M 3s
    ##  13850K .......... .......... .......... .......... .......... 10% 58.4M 3s
    ##  13900K .......... .......... .......... .......... .......... 10% 35.5M 3s
    ##  13950K .......... .......... .......... .......... .......... 10% 38.4M 3s
    ##  14000K .......... .......... .......... .......... .......... 10% 26.5M 3s
    ##  14050K .......... .......... .......... .......... .......... 10% 48.0M 3s
    ##  14100K .......... .......... .......... .......... .......... 10% 89.3M 3s
    ##  14150K .......... .......... .......... .......... .......... 10% 47.1M 3s
    ##  14200K .......... .......... .......... .......... .......... 10% 51.0M 3s
    ##  14250K .......... .......... .......... .......... .......... 10% 39.8M 3s
    ##  14300K .......... .......... .......... .......... .......... 10% 21.9M 3s
    ##  14350K .......... .......... .......... .......... .......... 10% 41.4M 3s
    ##  14400K .......... .......... .......... .......... .......... 10% 73.6M 3s
    ##  14450K .......... .......... .......... .......... .......... 10% 82.7M 3s
    ##  14500K .......... .......... .......... .......... .......... 10% 61.9M 3s
    ##  14550K .......... .......... .......... .......... .......... 10% 52.1M 3s
    ##  14600K .......... .......... .......... .......... .......... 10% 62.1M 3s
    ##  14650K .......... .......... .......... .......... .......... 10% 75.0M 3s
    ##  14700K .......... .......... .......... .......... .......... 11% 82.3M 3s
    ##  14750K .......... .......... .......... .......... .......... 11% 63.3M 3s
    ##  14800K .......... .......... .......... .......... .......... 11% 95.0M 3s
    ##  14850K .......... .......... .......... .......... .......... 11% 70.8M 3s
    ##  14900K .......... .......... .......... .......... .......... 11% 80.1M 3s
    ##  14950K .......... .......... .......... .......... .......... 11% 56.7M 3s
    ##  15000K .......... .......... .......... .......... .......... 11% 59.0M 3s
    ##  15050K .......... .......... .......... .......... .......... 11% 56.5M 3s
    ##  15100K .......... .......... .......... .......... .......... 11% 85.3M 3s
    ##  15150K .......... .......... .......... .......... .......... 11% 48.7M 3s
    ##  15200K .......... .......... .......... .......... .......... 11% 65.6M 3s
    ##  15250K .......... .......... .......... .......... .......... 11% 57.3M 3s
    ##  15300K .......... .......... .......... .......... .......... 11% 60.6M 3s
    ##  15350K .......... .......... .......... .......... .......... 11% 64.0M 3s
    ##  15400K .......... .......... .......... .......... .......... 11% 58.7M 3s
    ##  15450K .......... .......... .......... .......... .......... 11% 22.5M 3s
    ##  15500K .......... .......... .......... .......... .......... 11% 86.7M 3s
    ##  15550K .......... .......... .......... .......... .......... 11% 71.4M 3s
    ##  15600K .......... .......... .......... .......... .......... 11% 67.1M 3s
    ##  15650K .......... .......... .......... .......... .......... 11% 67.0M 3s
    ##  15700K .......... .......... .......... .......... .......... 11% 63.1M 3s
    ##  15750K .......... .......... .......... .......... .......... 11% 60.5M 3s
    ##  15800K .......... .......... .......... .......... .......... 11% 70.1M 3s
    ##  15850K .......... .......... .......... .......... .......... 11% 53.5M 3s
    ##  15900K .......... .......... .......... .......... .......... 11% 71.8M 3s
    ##  15950K .......... .......... .......... .......... .......... 11% 65.4M 3s
    ##  16000K .......... .......... .......... .......... .......... 11% 55.6M 3s
    ##  16050K .......... .......... .......... .......... .......... 12% 77.2M 3s
    ##  16100K .......... .......... .......... .......... .......... 12% 66.0M 3s
    ##  16150K .......... .......... .......... .......... .......... 12% 62.5M 3s
    ##  16200K .......... .......... .......... .......... .......... 12% 73.0M 3s
    ##  16250K .......... .......... .......... .......... .......... 12% 66.0M 3s
    ##  16300K .......... .......... .......... .......... .......... 12% 78.0M 3s
    ##  16350K .......... .......... .......... .......... .......... 12% 42.5M 3s
    ##  16400K .......... .......... .......... .......... .......... 12% 80.0M 3s
    ##  16450K .......... .......... .......... .......... .......... 12% 76.4M 3s
    ##  16500K .......... .......... .......... .......... .......... 12% 70.7M 3s
    ##  16550K .......... .......... .......... .......... .......... 12% 66.1M 3s
    ##  16600K .......... .......... .......... .......... .......... 12% 53.6M 3s
    ##  16650K .......... .......... .......... .......... .......... 12% 70.9M 3s
    ##  16700K .......... .......... .......... .......... .......... 12% 90.7M 3s
    ##  16750K .......... .......... .......... .......... .......... 12% 53.5M 3s
    ##  16800K .......... .......... .......... .......... .......... 12% 58.5M 3s
    ##  16850K .......... .......... .......... .......... .......... 12% 55.0M 3s
    ##  16900K .......... .......... .......... .......... .......... 12% 67.0M 3s
    ##  16950K .......... .......... .......... .......... .......... 12% 72.5M 3s
    ##  17000K .......... .......... .......... .......... .......... 12% 62.6M 3s
    ##  17050K .......... .......... .......... .......... .......... 12% 53.3M 3s
    ##  17100K .......... .......... .......... .......... .......... 12% 60.7M 3s
    ##  17150K .......... .......... .......... .......... .......... 12% 60.7M 3s
    ##  17200K .......... .......... .......... .......... .......... 12% 70.5M 3s
    ##  17250K .......... .......... .......... .......... .......... 12% 55.2M 3s
    ##  17300K .......... .......... .......... .......... .......... 12% 55.0M 3s
    ##  17350K .......... .......... .......... .......... .......... 12% 61.1M 3s
    ##  17400K .......... .......... .......... .......... .......... 13% 72.4M 3s
    ##  17450K .......... .......... .......... .......... .......... 13% 68.0M 3s
    ##  17500K .......... .......... .......... .......... .......... 13% 54.5M 3s
    ##  17550K .......... .......... .......... .......... .......... 13% 65.4M 3s
    ##  17600K .......... .......... .......... .......... .......... 13% 75.2M 3s
    ##  17650K .......... .......... .......... .......... .......... 13% 86.0M 3s
    ##  17700K .......... .......... .......... .......... .......... 13% 90.0M 3s
    ##  17750K .......... .......... .......... .......... .......... 13% 68.2M 3s
    ##  17800K .......... .......... .......... .......... .......... 13% 59.3M 3s
    ##  17850K .......... .......... .......... .......... .......... 13% 53.9M 3s
    ##  17900K .......... .......... .......... .......... .......... 13% 78.2M 3s
    ##  17950K .......... .......... .......... .......... .......... 13% 68.5M 3s
    ##  18000K .......... .......... .......... .......... .......... 13% 81.9M 3s
    ##  18050K .......... .......... .......... .......... .......... 13% 83.2M 3s
    ##  18100K .......... .......... .......... .......... .......... 13%  117M 3s
    ##  18150K .......... .......... .......... .......... .......... 13% 30.1M 3s
    ##  18200K .......... .......... .......... .......... .......... 13% 83.8M 3s
    ##  18250K .......... .......... .......... .......... .......... 13% 57.2M 3s
    ##  18300K .......... .......... .......... .......... .......... 13% 92.0M 3s
    ##  18350K .......... .......... .......... .......... .......... 13% 77.1M 3s
    ##  18400K .......... .......... .......... .......... .......... 13% 92.2M 3s
    ##  18450K .......... .......... .......... .......... .......... 13%  106M 3s
    ##  18500K .......... .......... .......... .......... .......... 13% 65.4M 3s
    ##  18550K .......... .......... .......... .......... .......... 13% 61.2M 3s
    ##  18600K .......... .......... .......... .......... .......... 13% 62.1M 3s
    ##  18650K .......... .......... .......... .......... .......... 13% 93.3M 3s
    ##  18700K .......... .......... .......... .......... .......... 13% 87.6M 3s
    ##  18750K .......... .......... .......... .......... .......... 14% 79.0M 3s
    ##  18800K .......... .......... .......... .......... .......... 14% 76.0M 3s
    ##  18850K .......... .......... .......... .......... .......... 14% 80.2M 3s
    ##  18900K .......... .......... .......... .......... .......... 14% 84.1M 3s
    ##  18950K .......... .......... .......... .......... .......... 14% 77.8M 3s
    ##  19000K .......... .......... .......... .......... .......... 14% 64.6M 3s
    ##  19050K .......... .......... .......... .......... .......... 14% 61.5M 3s
    ##  19100K .......... .......... .......... .......... .......... 14% 59.7M 3s
    ##  19150K .......... .......... .......... .......... .......... 14% 58.0M 3s
    ##  19200K .......... .......... .......... .......... .......... 14% 82.6M 3s
    ##  19250K .......... .......... .......... .......... .......... 14%  100M 3s
    ##  19300K .......... .......... .......... .......... .......... 14% 83.8M 3s
    ##  19350K .......... .......... .......... .......... .......... 14%  102M 3s
    ##  19400K .......... .......... .......... .......... .......... 14% 62.5M 3s
    ##  19450K .......... .......... .......... .......... .......... 14%  100M 3s
    ##  19500K .......... .......... .......... .......... .......... 14% 38.5M 3s
    ##  19550K .......... .......... .......... .......... .......... 14% 30.3M 3s
    ##  19600K .......... .......... .......... .......... .......... 14% 36.0M 3s
    ##  19650K .......... .......... .......... .......... .......... 14% 56.6M 3s
    ##  19700K .......... .......... .......... .......... .......... 14% 58.7M 3s
    ##  19750K .......... .......... .......... .......... .......... 14% 65.0M 3s
    ##  19800K .......... .......... .......... .......... .......... 14% 70.4M 3s
    ##  19850K .......... .......... .......... .......... .......... 14% 27.6M 3s
    ##  19900K .......... .......... .......... .......... .......... 14% 77.9M 3s
    ##  19950K .......... .......... .......... .......... .......... 14% 47.2M 3s
    ##  20000K .......... .......... .......... .......... .......... 14% 68.3M 3s
    ##  20050K .......... .......... .......... .......... .......... 14% 33.1M 3s
    ##  20100K .......... .......... .......... .......... .......... 15% 37.1M 3s
    ##  20150K .......... .......... .......... .......... .......... 15% 47.9M 3s
    ##  20200K .......... .......... .......... .......... .......... 15% 54.2M 3s
    ##  20250K .......... .......... .......... .......... .......... 15% 62.6M 3s
    ##  20300K .......... .......... .......... .......... .......... 15% 71.8M 3s
    ##  20350K .......... .......... .......... .......... .......... 15% 35.1M 3s
    ##  20400K .......... .......... .......... .......... .......... 15% 43.7M 3s
    ##  20450K .......... .......... .......... .......... .......... 15% 61.0M 3s
    ##  20500K .......... .......... .......... .......... .......... 15% 56.3M 3s
    ##  20550K .......... .......... .......... .......... .......... 15% 29.7M 3s
    ##  20600K .......... .......... .......... .......... .......... 15% 56.6M 3s
    ##  20650K .......... .......... .......... .......... .......... 15% 52.6M 3s
    ##  20700K .......... .......... .......... .......... .......... 15% 30.5M 3s
    ##  20750K .......... .......... .......... .......... .......... 15% 57.0M 3s
    ##  20800K .......... .......... .......... .......... .......... 15% 82.7M 3s
    ##  20850K .......... .......... .......... .......... .......... 15% 41.9M 3s
    ##  20900K .......... .......... .......... .......... .......... 15% 63.1M 3s
    ##  20950K .......... .......... .......... .......... .......... 15% 43.5M 3s
    ##  21000K .......... .......... .......... .......... .......... 15% 44.3M 3s
    ##  21050K .......... .......... .......... .......... .......... 15% 42.0M 3s
    ##  21100K .......... .......... .......... .......... .......... 15% 55.9M 3s
    ##  21150K .......... .......... .......... .......... .......... 15% 46.0M 3s
    ##  21200K .......... .......... .......... .......... .......... 15% 30.9M 3s
    ##  21250K .......... .......... .......... .......... .......... 15% 66.9M 3s
    ##  21300K .......... .......... .......... .......... .......... 15% 36.7M 3s
    ##  21350K .......... .......... .......... .......... .......... 15% 65.9M 3s
    ##  21400K .......... .......... .......... .......... .......... 15% 38.9M 3s
    ##  21450K .......... .......... .......... .......... .......... 16% 63.3M 3s
    ##  21500K .......... .......... .......... .......... .......... 16% 63.7M 3s
    ##  21550K .......... .......... .......... .......... .......... 16% 43.1M 3s
    ##  21600K .......... .......... .......... .......... .......... 16% 35.2M 3s
    ##  21650K .......... .......... .......... .......... .......... 16% 73.2M 3s
    ##  21700K .......... .......... .......... .......... .......... 16% 38.8M 3s
    ##  21750K .......... .......... .......... .......... .......... 16% 43.6M 3s
    ##  21800K .......... .......... .......... .......... .......... 16% 32.7M 3s
    ##  21850K .......... .......... .......... .......... .......... 16% 66.4M 3s
    ##  21900K .......... .......... .......... .......... .......... 16% 79.3M 3s
    ##  21950K .......... .......... .......... .......... .......... 16% 41.6M 3s
    ##  22000K .......... .......... .......... .......... .......... 16% 47.0M 3s
    ##  22050K .......... .......... .......... .......... .......... 16% 44.0M 3s
    ##  22100K .......... .......... .......... .......... .......... 16% 72.4M 3s
    ##  22150K .......... .......... .......... .......... .......... 16% 30.6M 3s
    ##  22200K .......... .......... .......... .......... .......... 16% 63.8M 3s
    ##  22250K .......... .......... .......... .......... .......... 16% 37.8M 3s
    ##  22300K .......... .......... .......... .......... .......... 16% 58.1M 3s
    ##  22350K .......... .......... .......... .......... .......... 16% 33.8M 3s
    ##  22400K .......... .......... .......... .......... .......... 16% 40.8M 3s
    ##  22450K .......... .......... .......... .......... .......... 16% 56.4M 3s
    ##  22500K .......... .......... .......... .......... .......... 16% 53.6M 3s
    ##  22550K .......... .......... .......... .......... .......... 16% 73.7M 3s
    ##  22600K .......... .......... .......... .......... .......... 16% 37.0M 3s
    ##  22650K .......... .......... .......... .......... .......... 16% 69.0M 3s
    ##  22700K .......... .......... .......... .......... .......... 16% 29.3M 3s
    ##  22750K .......... .......... .......... .......... .......... 17% 48.2M 3s
    ##  22800K .......... .......... .......... .......... .......... 17% 66.7M 3s
    ##  22850K .......... .......... .......... .......... .......... 17% 70.9M 3s
    ##  22900K .......... .......... .......... .......... .......... 17% 30.9M 3s
    ##  22950K .......... .......... .......... .......... .......... 17% 58.0M 3s
    ##  23000K .......... .......... .......... .......... .......... 17% 90.3M 3s
    ##  23050K .......... .......... .......... .......... .......... 17% 67.0M 3s
    ##  23100K .......... .......... .......... .......... .......... 17% 27.1M 3s
    ##  23150K .......... .......... .......... .......... .......... 17% 49.8M 3s
    ##  23200K .......... .......... .......... .......... .......... 17% 32.5M 3s
    ##  23250K .......... .......... .......... .......... .......... 17% 57.8M 3s
    ##  23300K .......... .......... .......... .......... .......... 17% 74.7M 3s
    ##  23350K .......... .......... .......... .......... .......... 17% 62.0M 3s
    ##  23400K .......... .......... .......... .......... .......... 17% 71.2M 3s
    ##  23450K .......... .......... .......... .......... .......... 17% 38.3M 3s
    ##  23500K .......... .......... .......... .......... .......... 17% 64.2M 3s
    ##  23550K .......... .......... .......... .......... .......... 17% 50.0M 3s
    ##  23600K .......... .......... .......... .......... .......... 17% 36.5M 3s
    ##  23650K .......... .......... .......... .......... .......... 17% 36.4M 3s
    ##  23700K .......... .......... .......... .......... .......... 17% 33.0M 3s
    ##  23750K .......... .......... .......... .......... .......... 17% 53.1M 3s
    ##  23800K .......... .......... .......... .......... .......... 17% 85.2M 3s
    ##  23850K .......... .......... .......... .......... .......... 17% 55.7M 3s
    ##  23900K .......... .......... .......... .......... .......... 17% 37.4M 3s
    ##  23950K .......... .......... .......... .......... .......... 17% 54.0M 3s
    ##  24000K .......... .......... .......... .......... .......... 17% 77.0M 3s
    ##  24050K .......... .......... .......... .......... .......... 17% 47.7M 3s
    ##  24100K .......... .......... .......... .......... .......... 18% 41.0M 3s
    ##  24150K .......... .......... .......... .......... .......... 18% 40.6M 3s
    ##  24200K .......... .......... .......... .......... .......... 18% 67.5M 3s
    ##  24250K .......... .......... .......... .......... .......... 18% 56.8M 3s
    ##  24300K .......... .......... .......... .......... .......... 18% 35.5M 3s
    ##  24350K .......... .......... .......... .......... .......... 18% 54.8M 3s
    ##  24400K .......... .......... .......... .......... .......... 18% 33.1M 3s
    ##  24450K .......... .......... .......... .......... .......... 18% 12.9M 3s
    ##  24500K .......... .......... .......... .......... .......... 18% 64.3M 3s
    ##  24550K .......... .......... .......... .......... .......... 18% 84.3M 3s
    ##  24600K .......... .......... .......... .......... .......... 18% 63.3M 3s
    ##  24650K .......... .......... .......... .......... .......... 18% 69.2M 3s
    ##  24700K .......... .......... .......... .......... .......... 18% 78.8M 3s
    ##  24750K .......... .......... .......... .......... .......... 18% 60.8M 3s
    ##  24800K .......... .......... .......... .......... .......... 18% 91.2M 3s
    ##  24850K .......... .......... .......... .......... .......... 18% 99.4M 3s
    ##  24900K .......... .......... .......... .......... .......... 18% 78.6M 3s
    ##  24950K .......... .......... .......... .......... .......... 18% 62.6M 3s
    ##  25000K .......... .......... .......... .......... .......... 18% 82.4M 3s
    ##  25050K .......... .......... .......... .......... .......... 18% 97.1M 3s
    ##  25100K .......... .......... .......... .......... .......... 18% 73.1M 3s
    ##  25150K .......... .......... .......... .......... .......... 18% 62.9M 3s
    ##  25200K .......... .......... .......... .......... .......... 18% 92.5M 3s
    ##  25250K .......... .......... .......... .......... .......... 18% 69.5M 2s
    ##  25300K .......... .......... .......... .......... .......... 18% 69.8M 2s
    ##  25350K .......... .......... .......... .......... .......... 18% 67.4M 2s
    ##  25400K .......... .......... .......... .......... .......... 18% 65.5M 2s
    ##  25450K .......... .......... .......... .......... .......... 19% 69.5M 2s
    ##  25500K .......... .......... .......... .......... .......... 19% 72.7M 2s
    ##  25550K .......... .......... .......... .......... .......... 19% 62.4M 2s
    ##  25600K .......... .......... .......... .......... .......... 19% 58.9M 2s
    ##  25650K .......... .......... .......... .......... .......... 19% 62.6M 2s
    ##  25700K .......... .......... .......... .......... .......... 19% 86.6M 2s
    ##  25750K .......... .......... .......... .......... .......... 19% 67.2M 2s
    ##  25800K .......... .......... .......... .......... .......... 19% 64.5M 2s
    ##  25850K .......... .......... .......... .......... .......... 19% 82.9M 2s
    ##  25900K .......... .......... .......... .......... .......... 19% 63.3M 2s
    ##  25950K .......... .......... .......... .......... .......... 19% 51.9M 2s
    ##  26000K .......... .......... .......... .......... .......... 19% 59.2M 2s
    ##  26050K .......... .......... .......... .......... .......... 19% 85.3M 2s
    ##  26100K .......... .......... .......... .......... .......... 19% 67.6M 2s
    ##  26150K .......... .......... .......... .......... .......... 19% 58.1M 2s
    ##  26200K .......... .......... .......... .......... .......... 19% 94.2M 2s
    ##  26250K .......... .......... .......... .......... .......... 19% 68.5M 2s
    ##  26300K .......... .......... .......... .......... .......... 19% 70.3M 2s
    ##  26350K .......... .......... .......... .......... .......... 19% 78.8M 2s
    ##  26400K .......... .......... .......... .......... .......... 19% 66.0M 2s
    ##  26450K .......... .......... .......... .......... .......... 19% 73.6M 2s
    ##  26500K .......... .......... .......... .......... .......... 19% 55.3M 2s
    ##  26550K .......... .......... .......... .......... .......... 19% 58.2M 2s
    ##  26600K .......... .......... .......... .......... .......... 19% 70.6M 2s
    ##  26650K .......... .......... .......... .......... .......... 19% 73.6M 2s
    ##  26700K .......... .......... .......... .......... .......... 19% 64.2M 2s
    ##  26750K .......... .......... .......... .......... .......... 19% 61.5M 2s
    ##  26800K .......... .......... .......... .......... .......... 20% 77.5M 2s
    ##  26850K .......... .......... .......... .......... .......... 20% 79.8M 2s
    ##  26900K .......... .......... .......... .......... .......... 20% 70.9M 2s
    ##  26950K .......... .......... .......... .......... .......... 20% 66.3M 2s
    ##  27000K .......... .......... .......... .......... .......... 20% 65.3M 2s
    ##  27050K .......... .......... .......... .......... .......... 20% 64.1M 2s
    ##  27100K .......... .......... .......... .......... .......... 20% 75.3M 2s
    ##  27150K .......... .......... .......... .......... .......... 20% 59.5M 2s
    ##  27200K .......... .......... .......... .......... .......... 20% 60.1M 2s
    ##  27250K .......... .......... .......... .......... .......... 20% 56.7M 2s
    ##  27300K .......... .......... .......... .......... .......... 20% 71.7M 2s
    ##  27350K .......... .......... .......... .......... .......... 20% 70.4M 2s
    ##  27400K .......... .......... .......... .......... .......... 20% 71.4M 2s
    ##  27450K .......... .......... .......... .......... .......... 20% 71.3M 2s
    ##  27500K .......... .......... .......... .......... .......... 20% 69.0M 2s
    ##  27550K .......... .......... .......... .......... .......... 20% 45.8M 2s
    ##  27600K .......... .......... .......... .......... .......... 20% 67.4M 2s
    ##  27650K .......... .......... .......... .......... .......... 20% 68.3M 2s
    ##  27700K .......... .......... .......... .......... .......... 20% 53.4M 2s
    ##  27750K .......... .......... .......... .......... .......... 20% 60.3M 2s
    ##  27800K .......... .......... .......... .......... .......... 20% 59.5M 2s
    ##  27850K .......... .......... .......... .......... .......... 20% 69.6M 2s
    ##  27900K .......... .......... .......... .......... .......... 20% 82.4M 2s
    ##  27950K .......... .......... .......... .......... .......... 20% 48.2M 2s
    ##  28000K .......... .......... .......... .......... .......... 20% 87.8M 2s
    ##  28050K .......... .......... .......... .......... .......... 20% 62.1M 2s
    ##  28100K .......... .......... .......... .......... .......... 20% 77.5M 2s
    ##  28150K .......... .......... .......... .......... .......... 21% 62.3M 2s
    ##  28200K .......... .......... .......... .......... .......... 21% 59.8M 2s
    ##  28250K .......... .......... .......... .......... .......... 21% 55.9M 2s
    ##  28300K .......... .......... .......... .......... .......... 21% 56.8M 2s
    ##  28350K .......... .......... .......... .......... .......... 21% 55.6M 2s
    ##  28400K .......... .......... .......... .......... .......... 21% 76.8M 2s
    ##  28450K .......... .......... .......... .......... .......... 21% 68.7M 2s
    ##  28500K .......... .......... .......... .......... .......... 21% 79.1M 2s
    ##  28550K .......... .......... .......... .......... .......... 21% 56.6M 2s
    ##  28600K .......... .......... .......... .......... .......... 21% 61.3M 2s
    ##  28650K .......... .......... .......... .......... .......... 21% 58.3M 2s
    ##  28700K .......... .......... .......... .......... .......... 21% 65.6M 2s
    ##  28750K .......... .......... .......... .......... .......... 21% 54.0M 2s
    ##  28800K .......... .......... .......... .......... .......... 21% 66.5M 2s
    ##  28850K .......... .......... .......... .......... .......... 21% 61.4M 2s
    ##  28900K .......... .......... .......... .......... .......... 21% 63.9M 2s
    ##  28950K .......... .......... .......... .......... .......... 21% 49.5M 2s
    ##  29000K .......... .......... .......... .......... .......... 21% 76.7M 2s
    ##  29050K .......... .......... .......... .......... .......... 21% 81.0M 2s
    ##  29100K .......... .......... .......... .......... .......... 21% 58.8M 2s
    ##  29150K .......... .......... .......... .......... .......... 21% 54.5M 2s
    ##  29200K .......... .......... .......... .......... .......... 21% 77.3M 2s
    ##  29250K .......... .......... .......... .......... .......... 21% 53.7M 2s
    ##  29300K .......... .......... .......... .......... .......... 21% 62.7M 2s
    ##  29350K .......... .......... .......... .......... .......... 21% 63.0M 2s
    ##  29400K .......... .......... .......... .......... .......... 21% 81.7M 2s
    ##  29450K .......... .......... .......... .......... .......... 22% 76.0M 2s
    ##  29500K .......... .......... .......... .......... .......... 22% 67.4M 2s
    ##  29550K .......... .......... .......... .......... .......... 22% 55.2M 2s
    ##  29600K .......... .......... .......... .......... .......... 22% 64.5M 2s
    ##  29650K .......... .......... .......... .......... .......... 22% 78.9M 2s
    ##  29700K .......... .......... .......... .......... .......... 22% 72.2M 2s
    ##  29750K .......... .......... .......... .......... .......... 22% 57.9M 2s
    ##  29800K .......... .......... .......... .......... .......... 22% 60.1M 2s
    ##  29850K .......... .......... .......... .......... .......... 22% 72.2M 2s
    ##  29900K .......... .......... .......... .......... .......... 22% 62.0M 2s
    ##  29950K .......... .......... .......... .......... .......... 22% 61.5M 2s
    ##  30000K .......... .......... .......... .......... .......... 22% 55.9M 2s
    ##  30050K .......... .......... .......... .......... .......... 22% 61.1M 2s
    ##  30100K .......... .......... .......... .......... .......... 22% 64.9M 2s
    ##  30150K .......... .......... .......... .......... .......... 22% 54.7M 2s
    ##  30200K .......... .......... .......... .......... .......... 22% 65.3M 2s
    ##  30250K .......... .......... .......... .......... .......... 22% 56.3M 2s
    ##  30300K .......... .......... .......... .......... .......... 22% 45.2M 2s
    ##  30350K .......... .......... .......... .......... .......... 22% 58.9M 2s
    ##  30400K .......... .......... .......... .......... .......... 22% 67.1M 2s
    ##  30450K .......... .......... .......... .......... .......... 22% 65.8M 2s
    ##  30500K .......... .......... .......... .......... .......... 22% 71.0M 2s
    ##  30550K .......... .......... .......... .......... .......... 22% 69.4M 2s
    ##  30600K .......... .......... .......... .......... .......... 22% 63.2M 2s
    ##  30650K .......... .......... .......... .......... .......... 22% 86.9M 2s
    ##  30700K .......... .......... .......... .......... .......... 22% 66.9M 2s
    ##  30750K .......... .......... .......... .......... .......... 22% 67.6M 2s
    ##  30800K .......... .......... .......... .......... .......... 23% 69.4M 2s
    ##  30850K .......... .......... .......... .......... .......... 23% 69.9M 2s
    ##  30900K .......... .......... .......... .......... .......... 23% 70.0M 2s
    ##  30950K .......... .......... .......... .......... .......... 23% 64.1M 2s
    ##  31000K .......... .......... .......... .......... .......... 23% 69.3M 2s
    ##  31050K .......... .......... .......... .......... .......... 23% 75.2M 2s
    ##  31100K .......... .......... .......... .......... .......... 23% 74.2M 2s
    ##  31150K .......... .......... .......... .......... .......... 23% 46.7M 2s
    ##  31200K .......... .......... .......... .......... .......... 23% 61.2M 2s
    ##  31250K .......... .......... .......... .......... .......... 23% 57.7M 2s
    ##  31300K .......... .......... .......... .......... .......... 23% 75.6M 2s
    ##  31350K .......... .......... .......... .......... .......... 23% 67.0M 2s
    ##  31400K .......... .......... .......... .......... .......... 23% 60.9M 2s
    ##  31450K .......... .......... .......... .......... .......... 23% 62.8M 2s
    ##  31500K .......... .......... .......... .......... .......... 23% 96.6M 2s
    ##  31550K .......... .......... .......... .......... .......... 23% 49.2M 2s
    ##  31600K .......... .......... .......... .......... .......... 23% 57.5M 2s
    ##  31650K .......... .......... .......... .......... .......... 23% 70.3M 2s
    ##  31700K .......... .......... .......... .......... .......... 23% 80.7M 2s
    ##  31750K .......... .......... .......... .......... .......... 23% 54.0M 2s
    ##  31800K .......... .......... .......... .......... .......... 23% 61.7M 2s
    ##  31850K .......... .......... .......... .......... .......... 23% 69.5M 2s
    ##  31900K .......... .......... .......... .......... .......... 23% 82.1M 2s
    ##  31950K .......... .......... .......... .......... .......... 23% 28.2M 2s
    ##  32000K .......... .......... .......... .......... .......... 23% 58.6M 2s
    ##  32050K .......... .......... .......... .......... .......... 23% 32.3M 2s
    ##  32100K .......... .......... .......... .......... .......... 23% 54.6M 2s
    ##  32150K .......... .......... .......... .......... .......... 24% 73.8M 2s
    ##  32200K .......... .......... .......... .......... .......... 24% 57.0M 2s
    ##  32250K .......... .......... .......... .......... .......... 24% 63.3M 2s
    ##  32300K .......... .......... .......... .......... .......... 24% 61.1M 2s
    ##  32350K .......... .......... .......... .......... .......... 24% 47.2M 2s
    ##  32400K .......... .......... .......... .......... .......... 24% 70.4M 2s
    ##  32450K .......... .......... .......... .......... .......... 24% 61.5M 2s
    ##  32500K .......... .......... .......... .......... .......... 24% 51.7M 2s
    ##  32550K .......... .......... .......... .......... .......... 24% 48.2M 2s
    ##  32600K .......... .......... .......... .......... .......... 24% 82.4M 2s
    ##  32650K .......... .......... .......... .......... .......... 24% 53.5M 2s
    ##  32700K .......... .......... .......... .......... .......... 24% 73.4M 2s
    ##  32750K .......... .......... .......... .......... .......... 24% 49.4M 2s
    ##  32800K .......... .......... .......... .......... .......... 24% 61.1M 2s
    ##  32850K .......... .......... .......... .......... .......... 24% 66.4M 2s
    ##  32900K .......... .......... .......... .......... .......... 24% 71.9M 2s
    ##  32950K .......... .......... .......... .......... .......... 24% 64.6M 2s
    ##  33000K .......... .......... .......... .......... .......... 24% 73.5M 2s
    ##  33050K .......... .......... .......... .......... .......... 24% 70.6M 2s
    ##  33100K .......... .......... .......... .......... .......... 24% 70.0M 2s
    ##  33150K .......... .......... .......... .......... .......... 24% 48.5M 2s
    ##  33200K .......... .......... .......... .......... .......... 24% 51.6M 2s
    ##  33250K .......... .......... .......... .......... .......... 24% 54.8M 2s
    ##  33300K .......... .......... .......... .......... .......... 24% 81.7M 2s
    ##  33350K .......... .......... .......... .......... .......... 24% 67.6M 2s
    ##  33400K .......... .......... .......... .......... .......... 24% 76.0M 2s
    ##  33450K .......... .......... .......... .......... .......... 24% 63.8M 2s
    ##  33500K .......... .......... .......... .......... .......... 25% 63.7M 2s
    ##  33550K .......... .......... .......... .......... .......... 25% 56.2M 2s
    ##  33600K .......... .......... .......... .......... .......... 25% 73.9M 2s
    ##  33650K .......... .......... .......... .......... .......... 25% 54.8M 2s
    ##  33700K .......... .......... .......... .......... .......... 25% 60.5M 2s
    ##  33750K .......... .......... .......... .......... .......... 25% 53.1M 2s
    ##  33800K .......... .......... .......... .......... .......... 25% 70.8M 2s
    ##  33850K .......... .......... .......... .......... .......... 25% 59.2M 2s
    ##  33900K .......... .......... .......... .......... .......... 25% 69.2M 2s
    ##  33950K .......... .......... .......... .......... .......... 25% 68.9M 2s
    ##  34000K .......... .......... .......... .......... .......... 25% 64.5M 2s
    ##  34050K .......... .......... .......... .......... .......... 25% 57.9M 2s
    ##  34100K .......... .......... .......... .......... .......... 25% 65.9M 2s
    ##  34150K .......... .......... .......... .......... .......... 25% 60.1M 2s
    ##  34200K .......... .......... .......... .......... .......... 25% 86.4M 2s
    ##  34250K .......... .......... .......... .......... .......... 25% 47.4M 2s
    ##  34300K .......... .......... .......... .......... .......... 25% 54.1M 2s
    ##  34350K .......... .......... .......... .......... .......... 25% 44.0M 2s
    ##  34400K .......... .......... .......... .......... .......... 25% 96.5M 2s
    ##  34450K .......... .......... .......... .......... .......... 25% 61.9M 2s
    ##  34500K .......... .......... .......... .......... .......... 25% 60.4M 2s
    ##  34550K .......... .......... .......... .......... .......... 25% 67.9M 2s
    ##  34600K .......... .......... .......... .......... .......... 25% 58.4M 2s
    ##  34650K .......... .......... .......... .......... .......... 25% 77.2M 2s
    ##  34700K .......... .......... .......... .......... .......... 25% 70.6M 2s
    ##  34750K .......... .......... .......... .......... .......... 25% 60.0M 2s
    ##  34800K .......... .......... .......... .......... .......... 25% 73.4M 2s
    ##  34850K .......... .......... .......... .......... .......... 26% 78.7M 2s
    ##  34900K .......... .......... .......... .......... .......... 26% 81.4M 2s
    ##  34950K .......... .......... .......... .......... .......... 26% 53.6M 2s
    ##  35000K .......... .......... .......... .......... .......... 26% 71.0M 2s
    ##  35050K .......... .......... .......... .......... .......... 26% 92.3M 2s
    ##  35100K .......... .......... .......... .......... .......... 26% 57.6M 2s
    ##  35150K .......... .......... .......... .......... .......... 26% 52.8M 2s
    ##  35200K .......... .......... .......... .......... .......... 26% 68.8M 2s
    ##  35250K .......... .......... .......... .......... .......... 26% 71.0M 2s
    ##  35300K .......... .......... .......... .......... .......... 26% 52.1M 2s
    ##  35350K .......... .......... .......... .......... .......... 26% 47.6M 2s
    ##  35400K .......... .......... .......... .......... .......... 26% 96.4M 2s
    ##  35450K .......... .......... .......... .......... .......... 26% 76.1M 2s
    ##  35500K .......... .......... .......... .......... .......... 26% 58.0M 2s
    ##  35550K .......... .......... .......... .......... .......... 26% 72.2M 2s
    ##  35600K .......... .......... .......... .......... .......... 26% 47.7M 2s
    ##  35650K .......... .......... .......... .......... .......... 26% 65.7M 2s
    ##  35700K .......... .......... .......... .......... .......... 26% 78.3M 2s
    ##  35750K .......... .......... .......... .......... .......... 26% 69.4M 2s
    ##  35800K .......... .......... .......... .......... .......... 26% 74.2M 2s
    ##  35850K .......... .......... .......... .......... .......... 26% 73.1M 2s
    ##  35900K .......... .......... .......... .......... .......... 26% 71.5M 2s
    ##  35950K .......... .......... .......... .......... .......... 26% 57.6M 2s
    ##  36000K .......... .......... .......... .......... .......... 26% 69.1M 2s
    ##  36050K .......... .......... .......... .......... .......... 26% 85.8M 2s
    ##  36100K .......... .......... .......... .......... .......... 26% 53.3M 2s
    ##  36150K .......... .......... .......... .......... .......... 27% 60.3M 2s
    ##  36200K .......... .......... .......... .......... .......... 27% 71.8M 2s
    ##  36250K .......... .......... .......... .......... .......... 27% 86.1M 2s
    ##  36300K .......... .......... .......... .......... .......... 27% 77.3M 2s
    ##  36350K .......... .......... .......... .......... .......... 27% 62.7M 2s
    ##  36400K .......... .......... .......... .......... .......... 27% 52.0M 2s
    ##  36450K .......... .......... .......... .......... .......... 27% 49.4M 2s
    ##  36500K .......... .......... .......... .......... .......... 27% 8.43M 2s
    ##  36550K .......... .......... .......... .......... .......... 27% 12.5M 2s
    ##  36600K .......... .......... .......... .......... .......... 27% 13.4M 2s
    ##  36650K .......... .......... .......... .......... .......... 27% 11.4M 2s
    ##  36700K .......... .......... .......... .......... .......... 27% 12.0M 2s
    ##  36750K .......... .......... .......... .......... .......... 27% 12.8M 2s
    ##  36800K .......... .......... .......... .......... .......... 27% 14.8M 2s
    ##  36850K .......... .......... .......... .......... .......... 27% 18.9M 2s
    ##  36900K .......... .......... .......... .......... .......... 27% 22.2M 2s
    ##  36950K .......... .......... .......... .......... .......... 27% 16.6M 2s
    ##  37000K .......... .......... .......... .......... .......... 27% 19.3M 2s
    ##  37050K .......... .......... .......... .......... .......... 27% 19.3M 2s
    ##  37100K .......... .......... .......... .......... .......... 27% 22.1M 2s
    ##  37150K .......... .......... .......... .......... .......... 27% 19.2M 2s
    ##  37200K .......... .......... .......... .......... .......... 27% 19.5M 2s
    ##  37250K .......... .......... .......... .......... .......... 27% 20.8M 2s
    ##  37300K .......... .......... .......... .......... .......... 27% 21.8M 2s
    ##  37350K .......... .......... .......... .......... .......... 27% 16.4M 2s
    ##  37400K .......... .......... .......... .......... .......... 27% 19.2M 2s
    ##  37450K .......... .......... .......... .......... .......... 27% 18.7M 2s
    ##  37500K .......... .......... .......... .......... .......... 28% 20.2M 2s
    ##  37550K .......... .......... .......... .......... .......... 28% 14.3M 2s
    ##  37600K .......... .......... .......... .......... .......... 28% 17.4M 2s
    ##  37650K .......... .......... .......... .......... .......... 28% 20.4M 2s
    ##  37700K .......... .......... .......... .......... .......... 28% 17.6M 2s
    ##  37750K .......... .......... .......... .......... .......... 28% 16.6M 2s
    ##  37800K .......... .......... .......... .......... .......... 28% 22.7M 2s
    ##  37850K .......... .......... .......... .......... .......... 28% 19.6M 2s
    ##  37900K .......... .......... .......... .......... .......... 28% 19.0M 2s
    ##  37950K .......... .......... .......... .......... .......... 28% 18.4M 2s
    ##  38000K .......... .......... .......... .......... .......... 28% 18.9M 2s
    ##  38050K .......... .......... .......... .......... .......... 28% 18.9M 2s
    ##  38100K .......... .......... .......... .......... .......... 28% 18.2M 2s
    ##  38150K .......... .......... .......... .......... .......... 28% 18.1M 2s
    ##  38200K .......... .......... .......... .......... .......... 28% 22.5M 2s
    ##  38250K .......... .......... .......... .......... .......... 28% 19.1M 2s
    ##  38300K .......... .......... .......... .......... .......... 28% 19.2M 2s
    ##  38350K .......... .......... .......... .......... .......... 28% 16.5M 2s
    ##  38400K .......... .......... .......... .......... .......... 28% 17.7M 2s
    ##  38450K .......... .......... .......... .......... .......... 28% 19.6M 2s
    ##  38500K .......... .......... .......... .......... .......... 28% 19.9M 2s
    ##  38550K .......... .......... .......... .......... .......... 28% 14.9M 2s
    ##  38600K .......... .......... .......... .......... .......... 28% 19.4M 2s
    ##  38650K .......... .......... .......... .......... .......... 28% 19.8M 2s
    ##  38700K .......... .......... .......... .......... .......... 28% 16.3M 2s
    ##  38750K .......... .......... .......... .......... .......... 28% 20.5M 2s
    ##  38800K .......... .......... .......... .......... .......... 28% 18.3M 2s
    ##  38850K .......... .......... .......... .......... .......... 29% 16.8M 2s
    ##  38900K .......... .......... .......... .......... .......... 29% 19.0M 2s
    ##  38950K .......... .......... .......... .......... .......... 29% 14.9M 2s
    ##  39000K .......... .......... .......... .......... .......... 29% 17.6M 2s
    ##  39050K .......... .......... .......... .......... .......... 29% 19.8M 2s
    ##  39100K .......... .......... .......... .......... .......... 29% 21.6M 2s
    ##  39150K .......... .......... .......... .......... .......... 29% 12.6M 2s
    ##  39200K .......... .......... .......... .......... .......... 29% 31.1M 2s
    ##  39250K .......... .......... .......... .......... .......... 29% 35.6M 2s
    ##  39300K .......... .......... .......... .......... .......... 29% 28.5M 2s
    ##  39350K .......... .......... .......... .......... .......... 29% 27.0M 2s
    ##  39400K .......... .......... .......... .......... .......... 29% 31.0M 2s
    ##  39450K .......... .......... .......... .......... .......... 29% 26.9M 2s
    ##  39500K .......... .......... .......... .......... .......... 29% 26.2M 2s
    ##  39550K .......... .......... .......... .......... .......... 29% 23.9M 2s
    ##  39600K .......... .......... .......... .......... .......... 29% 29.1M 2s
    ##  39650K .......... .......... .......... .......... .......... 29% 27.6M 2s
    ##  39700K .......... .......... .......... .......... .......... 29% 30.9M 2s
    ##  39750K .......... .......... .......... .......... .......... 29% 26.3M 2s
    ##  39800K .......... .......... .......... .......... .......... 29% 27.0M 2s
    ##  39850K .......... .......... .......... .......... .......... 29% 37.1M 2s
    ##  39900K .......... .......... .......... .......... .......... 29% 27.2M 2s
    ##  39950K .......... .......... .......... .......... .......... 29% 26.1M 2s
    ##  40000K .......... .......... .......... .......... .......... 29% 31.7M 2s
    ##  40050K .......... .......... .......... .......... .......... 29% 28.8M 2s
    ##  40100K .......... .......... .......... .......... .......... 29% 29.8M 2s
    ##  40150K .......... .......... .......... .......... .......... 29% 20.9M 2s
    ##  40200K .......... .......... .......... .......... .......... 30% 32.9M 2s
    ##  40250K .......... .......... .......... .......... .......... 30% 27.3M 2s
    ##  40300K .......... .......... .......... .......... .......... 30% 35.8M 2s
    ##  40350K .......... .......... .......... .......... .......... 30% 22.1M 2s
    ##  40400K .......... .......... .......... .......... .......... 30% 33.0M 2s
    ##  40450K .......... .......... .......... .......... .......... 30% 27.2M 2s
    ##  40500K .......... .......... .......... .......... .......... 30% 27.8M 2s
    ##  40550K .......... .......... .......... .......... .......... 30% 30.0M 2s
    ##  40600K .......... .......... .......... .......... .......... 30% 27.7M 2s
    ##  40650K .......... .......... .......... .......... .......... 30% 27.2M 2s
    ##  40700K .......... .......... .......... .......... .......... 30% 27.1M 2s
    ##  40750K .......... .......... .......... .......... .......... 30% 24.9M 2s
    ##  40800K .......... .......... .......... .......... .......... 30% 37.5M 2s
    ##  40850K .......... .......... .......... .......... .......... 30% 32.8M 2s
    ##  40900K .......... .......... .......... .......... .......... 30% 40.3M 2s
    ##  40950K .......... .......... .......... .......... .......... 30% 23.4M 2s
    ##  41000K .......... .......... .......... .......... .......... 30% 38.3M 2s
    ##  41050K .......... .......... .......... .......... .......... 30% 25.6M 2s
    ##  41100K .......... .......... .......... .......... .......... 30% 37.9M 2s
    ##  41150K .......... .......... .......... .......... .......... 30% 26.8M 2s
    ##  41200K .......... .......... .......... .......... .......... 30% 37.0M 2s
    ##  41250K .......... .......... .......... .......... .......... 30% 29.5M 2s
    ##  41300K .......... .......... .......... .......... .......... 30% 34.8M 2s
    ##  41350K .......... .......... .......... .......... .......... 30% 26.3M 2s
    ##  41400K .......... .......... .......... .......... .......... 30% 40.3M 2s
    ##  41450K .......... .......... .......... .......... .......... 30% 32.0M 2s
    ##  41500K .......... .......... .......... .......... .......... 30% 30.9M 2s
    ##  41550K .......... .......... .......... .......... .......... 31% 22.4M 2s
    ##  41600K .......... .......... .......... .......... .......... 31% 29.7M 2s
    ##  41650K .......... .......... .......... .......... .......... 31% 29.3M 2s
    ##  41700K .......... .......... .......... .......... .......... 31% 29.6M 2s
    ##  41750K .......... .......... .......... .......... .......... 31% 26.0M 2s
    ##  41800K .......... .......... .......... .......... .......... 31% 30.6M 2s
    ##  41850K .......... .......... .......... .......... .......... 31% 30.4M 2s
    ##  41900K .......... .......... .......... .......... .......... 31% 29.3M 2s
    ##  41950K .......... .......... .......... .......... .......... 31% 27.1M 2s
    ##  42000K .......... .......... .......... .......... .......... 31% 33.3M 2s
    ##  42050K .......... .......... .......... .......... .......... 31% 33.4M 2s
    ##  42100K .......... .......... .......... .......... .......... 31% 22.6M 2s
    ##  42150K .......... .......... .......... .......... .......... 31% 26.8M 2s
    ##  42200K .......... .......... .......... .......... .......... 31% 28.3M 2s
    ##  42250K .......... .......... .......... .......... .......... 31% 36.1M 2s
    ##  42300K .......... .......... .......... .......... .......... 31% 28.9M 2s
    ##  42350K .......... .......... .......... .......... .......... 31% 33.0M 2s
    ##  42400K .......... .......... .......... .......... .......... 31% 26.6M 2s
    ##  42450K .......... .......... .......... .......... .......... 31% 38.6M 2s
    ##  42500K .......... .......... .......... .......... .......... 31% 28.9M 2s
    ##  42550K .......... .......... .......... .......... .......... 31% 32.6M 2s
    ##  42600K .......... .......... .......... .......... .......... 31% 26.9M 2s
    ##  42650K .......... .......... .......... .......... .......... 31% 33.9M 2s
    ##  42700K .......... .......... .......... .......... .......... 31% 24.7M 2s
    ##  42750K .......... .......... .......... .......... .......... 31% 31.2M 2s
    ##  42800K .......... .......... .......... .......... .......... 31% 29.1M 2s
    ##  42850K .......... .......... .......... .......... .......... 31% 39.3M 2s
    ##  42900K .......... .......... .......... .......... .......... 32% 28.5M 2s
    ##  42950K .......... .......... .......... .......... .......... 32% 31.8M 2s
    ##  43000K .......... .......... .......... .......... .......... 32% 28.9M 2s
    ##  43050K .......... .......... .......... .......... .......... 32% 38.5M 2s
    ##  43100K .......... .......... .......... .......... .......... 32% 25.6M 2s
    ##  43150K .......... .......... .......... .......... .......... 32% 31.5M 2s
    ##  43200K .......... .......... .......... .......... .......... 32% 29.4M 2s
    ##  43250K .......... .......... .......... .......... .......... 32% 32.1M 2s
    ##  43300K .......... .......... .......... .......... .......... 32% 30.5M 2s
    ##  43350K .......... .......... .......... .......... .......... 32% 22.9M 2s
    ##  43400K .......... .......... .......... .......... .......... 32% 32.6M 2s
    ##  43450K .......... .......... .......... .......... .......... 32% 28.8M 2s
    ##  43500K .......... .......... .......... .......... .......... 32% 35.1M 2s
    ##  43550K .......... .......... .......... .......... .......... 32% 27.8M 2s
    ##  43600K .......... .......... .......... .......... .......... 32% 30.6M 2s
    ##  43650K .......... .......... .......... .......... .......... 32% 31.4M 2s
    ##  43700K .......... .......... .......... .......... .......... 32% 30.5M 2s
    ##  43750K .......... .......... .......... .......... .......... 32% 25.7M 2s
    ##  43800K .......... .......... .......... .......... .......... 32% 27.7M 2s
    ##  43850K .......... .......... .......... .......... .......... 32% 28.5M 2s
    ##  43900K .......... .......... .......... .......... .......... 32% 32.5M 2s
    ##  43950K .......... .......... .......... .......... .......... 32% 25.6M 2s
    ##  44000K .......... .......... .......... .......... .......... 32% 33.2M 2s
    ##  44050K .......... .......... .......... .......... .......... 32% 29.7M 2s
    ##  44100K .......... .......... .......... .......... .......... 32% 37.9M 2s
    ##  44150K .......... .......... .......... .......... .......... 32% 28.2M 2s
    ##  44200K .......... .......... .......... .......... .......... 33% 38.5M 2s
    ##  44250K .......... .......... .......... .......... .......... 33% 27.3M 2s
    ##  44300K .......... .......... .......... .......... .......... 33% 37.9M 2s
    ##  44350K .......... .......... .......... .......... .......... 33% 21.8M 2s
    ##  44400K .......... .......... .......... .......... .......... 33% 34.9M 2s
    ##  44450K .......... .......... .......... .......... .......... 33% 36.4M 2s
    ##  44500K .......... .......... .......... .......... .......... 33% 27.4M 2s
    ##  44550K .......... .......... .......... .......... .......... 33% 33.7M 2s
    ##  44600K .......... .......... .......... .......... .......... 33% 25.9M 2s
    ##  44650K .......... .......... .......... .......... .......... 33% 39.8M 2s
    ##  44700K .......... .......... .......... .......... .......... 33% 26.6M 2s
    ##  44750K .......... .......... .......... .......... .......... 33% 33.4M 2s
    ##  44800K .......... .......... .......... .......... .......... 33% 25.6M 2s
    ##  44850K .......... .......... .......... .......... .......... 33% 39.8M 2s
    ##  44900K .......... .......... .......... .......... .......... 33% 29.4M 2s
    ##  44950K .......... .......... .......... .......... .......... 33% 33.4M 2s
    ##  45000K .......... .......... .......... .......... .......... 33% 27.3M 2s
    ##  45050K .......... .......... .......... .......... .......... 33% 36.0M 2s
    ##  45100K .......... .......... .......... .......... .......... 33% 28.8M 2s
    ##  45150K .......... .......... .......... .......... .......... 33% 31.8M 2s
    ##  45200K .......... .......... .......... .......... .......... 33% 29.7M 2s
    ##  45250K .......... .......... .......... .......... .......... 33% 30.6M 2s
    ##  45300K .......... .......... .......... .......... .......... 33% 30.3M 2s
    ##  45350K .......... .......... .......... .......... .......... 33% 28.9M 2s
    ##  45400K .......... .......... .......... .......... .......... 33% 30.2M 2s
    ##  45450K .......... .......... .......... .......... .......... 33% 27.6M 2s
    ##  45500K .......... .......... .......... .......... .......... 33% 35.4M 2s
    ##  45550K .......... .......... .......... .......... .......... 34% 25.1M 2s
    ##  45600K .......... .......... .......... .......... .......... 34% 31.2M 2s
    ##  45650K .......... .......... .......... .......... .......... 34% 29.1M 2s
    ##  45700K .......... .......... .......... .......... .......... 34% 28.7M 2s
    ##  45750K .......... .......... .......... .......... .......... 34% 26.5M 2s
    ##  45800K .......... .......... .......... .......... .......... 34% 26.4M 2s
    ##  45850K .......... .......... .......... .......... .......... 34% 39.4M 2s
    ##  45900K .......... .......... .......... .......... .......... 34% 29.0M 2s
    ##  45950K .......... .......... .......... .......... .......... 34% 24.8M 2s
    ##  46000K .......... .......... .......... .......... .......... 34% 31.9M 2s
    ##  46050K .......... .......... .......... .......... .......... 34% 32.6M 2s
    ##  46100K .......... .......... .......... .......... .......... 34% 31.1M 2s
    ##  46150K .......... .......... .......... .......... .......... 34% 26.5M 2s
    ##  46200K .......... .......... .......... .......... .......... 34% 32.7M 2s
    ##  46250K .......... .......... .......... .......... .......... 34% 29.4M 2s
    ##  46300K .......... .......... .......... .......... .......... 34% 29.9M 2s
    ##  46350K .......... .......... .......... .......... .......... 34% 24.7M 2s
    ##  46400K .......... .......... .......... .......... .......... 34% 30.6M 2s
    ##  46450K .......... .......... .......... .......... .......... 34% 32.5M 2s
    ##  46500K .......... .......... .......... .......... .......... 34% 35.6M 2s
    ##  46550K .......... .......... .......... .......... .......... 34% 23.1M 2s
    ##  46600K .......... .......... .......... .......... .......... 34% 33.9M 2s
    ##  46650K .......... .......... .......... .......... .......... 34% 29.8M 2s
    ##  46700K .......... .......... .......... .......... .......... 34% 32.9M 2s
    ##  46750K .......... .......... .......... .......... .......... 34% 24.2M 2s
    ##  46800K .......... .......... .......... .......... .......... 34% 25.7M 2s
    ##  46850K .......... .......... .......... .......... .......... 34% 34.6M 2s
    ##  46900K .......... .......... .......... .......... .......... 35% 5.06M 2s
    ##  46950K .......... .......... .......... .......... .......... 35% 32.5M 2s
    ##  47000K .......... .......... .......... .......... .......... 35% 39.6M 2s
    ##  47050K .......... .......... .......... .......... .......... 35% 31.3M 2s
    ##  47100K .......... .......... .......... .......... .......... 35% 37.0M 2s
    ##  47150K .......... .......... .......... .......... .......... 35% 28.0M 2s
    ##  47200K .......... .......... .......... .......... .......... 35% 32.7M 2s
    ##  47250K .......... .......... .......... .......... .......... 35% 31.7M 2s
    ##  47300K .......... .......... .......... .......... .......... 35% 36.9M 2s
    ##  47350K .......... .......... .......... .......... .......... 35% 33.3M 2s
    ##  47400K .......... .......... .......... .......... .......... 35% 27.5M 2s
    ##  47450K .......... .......... .......... .......... .......... 35% 37.2M 2s
    ##  47500K .......... .......... .......... .......... .......... 35% 26.9M 2s
    ##  47550K .......... .......... .......... .......... .......... 35% 36.8M 2s
    ##  47600K .......... .......... .......... .......... .......... 35% 69.9M 2s
    ##  47650K .......... .......... .......... .......... .......... 35% 74.0M 2s
    ##  47700K .......... .......... .......... .......... .......... 35% 73.8M 2s
    ##  47750K .......... .......... .......... .......... .......... 35% 76.8M 2s
    ##  47800K .......... .......... .......... .......... .......... 35% 76.1M 2s
    ##  47850K .......... .......... .......... .......... .......... 35% 89.0M 2s
    ##  47900K .......... .......... .......... .......... .......... 35% 77.9M 2s
    ##  47950K .......... .......... .......... .......... .......... 35% 64.0M 2s
    ##  48000K .......... .......... .......... .......... .......... 35% 97.8M 2s
    ##  48050K .......... .......... .......... .......... .......... 35% 76.0M 2s
    ##  48100K .......... .......... .......... .......... .......... 35% 76.3M 2s
    ##  48150K .......... .......... .......... .......... .......... 35% 61.0M 2s
    ##  48200K .......... .......... .......... .......... .......... 35% 78.7M 2s
    ##  48250K .......... .......... .......... .......... .......... 36% 80.7M 2s
    ##  48300K .......... .......... .......... .......... .......... 36% 87.3M 2s
    ##  48350K .......... .......... .......... .......... .......... 36% 59.3M 2s
    ##  48400K .......... .......... .......... .......... .......... 36% 53.1M 2s
    ##  48450K .......... .......... .......... .......... .......... 36% 63.4M 2s
    ##  48500K .......... .......... .......... .......... .......... 36% 61.8M 2s
    ##  48550K .......... .......... .......... .......... .......... 36% 58.3M 2s
    ##  48600K .......... .......... .......... .......... .......... 36% 71.2M 2s
    ##  48650K .......... .......... .......... .......... .......... 36% 73.4M 2s
    ##  48700K .......... .......... .......... .......... .......... 36% 65.7M 2s
    ##  48750K .......... .......... .......... .......... .......... 36% 65.6M 2s
    ##  48800K .......... .......... .......... .......... .......... 36% 69.0M 2s
    ##  48850K .......... .......... .......... .......... .......... 36% 63.5M 2s
    ##  48900K .......... .......... .......... .......... .......... 36% 67.4M 2s
    ##  48950K .......... .......... .......... .......... .......... 36% 62.3M 2s
    ##  49000K .......... .......... .......... .......... .......... 36% 80.5M 2s
    ##  49050K .......... .......... .......... .......... .......... 36% 76.1M 2s
    ##  49100K .......... .......... .......... .......... .......... 36% 63.0M 2s
    ##  49150K .......... .......... .......... .......... .......... 36% 49.7M 2s
    ##  49200K .......... .......... .......... .......... .......... 36% 83.0M 2s
    ##  49250K .......... .......... .......... .......... .......... 36% 55.7M 2s
    ##  49300K .......... .......... .......... .......... .......... 36% 63.4M 2s
    ##  49350K .......... .......... .......... .......... .......... 36% 51.1M 2s
    ##  49400K .......... .......... .......... .......... .......... 36% 58.4M 2s
    ##  49450K .......... .......... .......... .......... .......... 36% 63.4M 2s
    ##  49500K .......... .......... .......... .......... .......... 36% 83.0M 2s
    ##  49550K .......... .......... .......... .......... .......... 36% 63.0M 2s
    ##  49600K .......... .......... .......... .......... .......... 37% 70.8M 2s
    ##  49650K .......... .......... .......... .......... .......... 37% 67.3M 2s
    ##  49700K .......... .......... .......... .......... .......... 37% 78.4M 2s
    ##  49750K .......... .......... .......... .......... .......... 37% 63.7M 2s
    ##  49800K .......... .......... .......... .......... .......... 37% 67.0M 2s
    ##  49850K .......... .......... .......... .......... .......... 37% 68.5M 2s
    ##  49900K .......... .......... .......... .......... .......... 37% 87.6M 2s
    ##  49950K .......... .......... .......... .......... .......... 37% 55.4M 2s
    ##  50000K .......... .......... .......... .......... .......... 37% 64.9M 2s
    ##  50050K .......... .......... .......... .......... .......... 37% 52.4M 2s
    ##  50100K .......... .......... .......... .......... .......... 37% 72.5M 2s
    ##  50150K .......... .......... .......... .......... .......... 37% 70.4M 2s
    ##  50200K .......... .......... .......... .......... .......... 37% 60.5M 2s
    ##  50250K .......... .......... .......... .......... .......... 37% 59.2M 2s
    ##  50300K .......... .......... .......... .......... .......... 37% 76.0M 2s
    ##  50350K .......... .......... .......... .......... .......... 37% 61.5M 2s
    ##  50400K .......... .......... .......... .......... .......... 37% 61.5M 2s
    ##  50450K .......... .......... .......... .......... .......... 37% 66.9M 2s
    ##  50500K .......... .......... .......... .......... .......... 37% 90.7M 2s
    ##  50550K .......... .......... .......... .......... .......... 37% 56.7M 2s
    ##  50600K .......... .......... .......... .......... .......... 37% 57.8M 2s
    ##  50650K .......... .......... .......... .......... .......... 37% 60.3M 2s
    ##  50700K .......... .......... .......... .......... .......... 37% 84.7M 2s
    ##  50750K .......... .......... .......... .......... .......... 37% 62.7M 2s
    ##  50800K .......... .......... .......... .......... .......... 37% 72.2M 2s
    ##  50850K .......... .......... .......... .......... .......... 37% 65.2M 2s
    ##  50900K .......... .......... .......... .......... .......... 38% 57.9M 2s
    ##  50950K .......... .......... .......... .......... .......... 38% 62.5M 2s
    ##  51000K .......... .......... .......... .......... .......... 38% 96.5M 2s
    ##  51050K .......... .......... .......... .......... .......... 38% 73.6M 2s
    ##  51100K .......... .......... .......... .......... .......... 38% 56.5M 2s
    ##  51150K .......... .......... .......... .......... .......... 38% 58.3M 2s
    ##  51200K .......... .......... .......... .......... .......... 38% 70.2M 2s
    ##  51250K .......... .......... .......... .......... .......... 38% 67.4M 2s
    ##  51300K .......... .......... .......... .......... .......... 38% 65.4M 2s
    ##  51350K .......... .......... .......... .......... .......... 38% 62.0M 2s
    ##  51400K .......... .......... .......... .......... .......... 38% 55.9M 2s
    ##  51450K .......... .......... .......... .......... .......... 38% 71.1M 2s
    ##  51500K .......... .......... .......... .......... .......... 38% 71.1M 2s
    ##  51550K .......... .......... .......... .......... .......... 38% 60.7M 2s
    ##  51600K .......... .......... .......... .......... .......... 38% 78.6M 2s
    ##  51650K .......... .......... .......... .......... .......... 38% 77.7M 2s
    ##  51700K .......... .......... .......... .......... .......... 38% 84.6M 2s
    ##  51750K .......... .......... .......... .......... .......... 38% 63.1M 2s
    ##  51800K .......... .......... .......... .......... .......... 38% 77.8M 2s
    ##  51850K .......... .......... .......... .......... .......... 38% 75.9M 2s
    ##  51900K .......... .......... .......... .......... .......... 38% 59.9M 2s
    ##  51950K .......... .......... .......... .......... .......... 38% 46.1M 2s
    ##  52000K .......... .......... .......... .......... .......... 38% 72.9M 2s
    ##  52050K .......... .......... .......... .......... .......... 38% 71.6M 2s
    ##  52100K .......... .......... .......... .......... .......... 38% 74.6M 2s
    ##  52150K .......... .......... .......... .......... .......... 38% 72.3M 2s
    ##  52200K .......... .......... .......... .......... .......... 38% 91.5M 2s
    ##  52250K .......... .......... .......... .......... .......... 39% 89.4M 2s
    ##  52300K .......... .......... .......... .......... .......... 39% 82.1M 2s
    ##  52350K .......... .......... .......... .......... .......... 39% 61.9M 2s
    ##  52400K .......... .......... .......... .......... .......... 39% 59.8M 2s
    ##  52450K .......... .......... .......... .......... .......... 39% 67.9M 2s
    ##  52500K .......... .......... .......... .......... .......... 39% 93.1M 2s
    ##  52550K .......... .......... .......... .......... .......... 39% 78.5M 2s
    ##  52600K .......... .......... .......... .......... .......... 39% 71.3M 2s
    ##  52650K .......... .......... .......... .......... .......... 39% 62.4M 2s
    ##  52700K .......... .......... .......... .......... .......... 39% 68.1M 2s
    ##  52750K .......... .......... .......... .......... .......... 39% 66.4M 2s
    ##  52800K .......... .......... .......... .......... .......... 39% 60.3M 2s
    ##  52850K .......... .......... .......... .......... .......... 39% 62.0M 2s
    ##  52900K .......... .......... .......... .......... .......... 39% 66.5M 2s
    ##  52950K .......... .......... .......... .......... .......... 39% 63.5M 2s
    ##  53000K .......... .......... .......... .......... .......... 39% 79.4M 2s
    ##  53050K .......... .......... .......... .......... .......... 39% 95.8M 2s
    ##  53100K .......... .......... .......... .......... .......... 39% 73.9M 2s
    ##  53150K .......... .......... .......... .......... .......... 39% 65.1M 2s
    ##  53200K .......... .......... .......... .......... .......... 39% 73.1M 2s
    ##  53250K .......... .......... .......... .......... .......... 39% 64.4M 2s
    ##  53300K .......... .......... .......... .......... .......... 39% 71.3M 2s
    ##  53350K .......... .......... .......... .......... .......... 39% 52.5M 2s
    ##  53400K .......... .......... .......... .......... .......... 39% 78.5M 2s
    ##  53450K .......... .......... .......... .......... .......... 39% 64.5M 2s
    ##  53500K .......... .......... .......... .......... .......... 39% 57.6M 2s
    ##  53550K .......... .......... .......... .......... .......... 39% 59.8M 2s
    ##  53600K .......... .......... .......... .......... .......... 40% 61.4M 2s
    ##  53650K .......... .......... .......... .......... .......... 40% 82.3M 2s
    ##  53700K .......... .......... .......... .......... .......... 40% 80.5M 2s
    ##  53750K .......... .......... .......... .......... .......... 40% 56.8M 2s
    ##  53800K .......... .......... .......... .......... .......... 40% 85.9M 2s
    ##  53850K .......... .......... .......... .......... .......... 40% 80.9M 2s
    ##  53900K .......... .......... .......... .......... .......... 40% 66.5M 2s
    ##  53950K .......... .......... .......... .......... .......... 40% 54.9M 2s
    ##  54000K .......... .......... .......... .......... .......... 40% 88.3M 2s
    ##  54050K .......... .......... .......... .......... .......... 40% 79.4M 2s
    ##  54100K .......... .......... .......... .......... .......... 40% 71.8M 2s
    ##  54150K .......... .......... .......... .......... .......... 40% 66.6M 2s
    ##  54200K .......... .......... .......... .......... .......... 40% 80.4M 2s
    ##  54250K .......... .......... .......... .......... .......... 40% 61.8M 2s
    ##  54300K .......... .......... .......... .......... .......... 40% 73.6M 2s
    ##  54350K .......... .......... .......... .......... .......... 40% 55.9M 2s
    ##  54400K .......... .......... .......... .......... .......... 40% 64.3M 2s
    ##  54450K .......... .......... .......... .......... .......... 40% 73.3M 2s
    ##  54500K .......... .......... .......... .......... .......... 40% 62.7M 2s
    ##  54550K .......... .......... .......... .......... .......... 40% 69.4M 2s
    ##  54600K .......... .......... .......... .......... .......... 40% 80.9M 2s
    ##  54650K .......... .......... .......... .......... .......... 40% 64.4M 2s
    ##  54700K .......... .......... .......... .......... .......... 40% 48.5M 2s
    ##  54750K .......... .......... .......... .......... .......... 40% 65.9M 2s
    ##  54800K .......... .......... .......... .......... .......... 40% 85.0M 2s
    ##  54850K .......... .......... .......... .......... .......... 40% 82.3M 2s
    ##  54900K .......... .......... .......... .......... .......... 40% 77.8M 2s
    ##  54950K .......... .......... .......... .......... .......... 41% 68.3M 2s
    ##  55000K .......... .......... .......... .......... .......... 41% 71.4M 2s
    ##  55050K .......... .......... .......... .......... .......... 41% 59.7M 2s
    ##  55100K .......... .......... .......... .......... .......... 41% 62.7M 2s
    ##  55150K .......... .......... .......... .......... .......... 41% 60.0M 2s
    ##  55200K .......... .......... .......... .......... .......... 41% 76.3M 2s
    ##  55250K .......... .......... .......... .......... .......... 41% 75.6M 2s
    ##  55300K .......... .......... .......... .......... .......... 41% 59.3M 2s
    ##  55350K .......... .......... .......... .......... .......... 41% 60.3M 2s
    ##  55400K .......... .......... .......... .......... .......... 41% 58.2M 2s
    ##  55450K .......... .......... .......... .......... .......... 41% 71.4M 2s
    ##  55500K .......... .......... .......... .......... .......... 41% 71.9M 2s
    ##  55550K .......... .......... .......... .......... .......... 41% 56.3M 2s
    ##  55600K .......... .......... .......... .......... .......... 41% 82.1M 2s
    ##  55650K .......... .......... .......... .......... .......... 41% 70.2M 2s
    ##  55700K .......... .......... .......... .......... .......... 41% 67.5M 2s
    ##  55750K .......... .......... .......... .......... .......... 41% 66.0M 2s
    ##  55800K .......... .......... .......... .......... .......... 41% 75.5M 2s
    ##  55850K .......... .......... .......... .......... .......... 41% 67.4M 2s
    ##  55900K .......... .......... .......... .......... .......... 41% 62.2M 2s
    ##  55950K .......... .......... .......... .......... .......... 41% 57.0M 2s
    ##  56000K .......... .......... .......... .......... .......... 41% 77.2M 2s
    ##  56050K .......... .......... .......... .......... .......... 41% 74.1M 2s
    ##  56100K .......... .......... .......... .......... .......... 41% 75.8M 2s
    ##  56150K .......... .......... .......... .......... .......... 41% 67.8M 2s
    ##  56200K .......... .......... .......... .......... .......... 41% 77.1M 2s
    ##  56250K .......... .......... .......... .......... .......... 41% 71.3M 2s
    ##  56300K .......... .......... .......... .......... .......... 42% 86.9M 2s
    ##  56350K .......... .......... .......... .......... .......... 42% 69.5M 2s
    ##  56400K .......... .......... .......... .......... .......... 42% 55.7M 2s
    ##  56450K .......... .......... .......... .......... .......... 42% 62.0M 2s
    ##  56500K .......... .......... .......... .......... .......... 42% 63.2M 2s
    ##  56550K .......... .......... .......... .......... .......... 42% 61.9M 2s
    ##  56600K .......... .......... .......... .......... .......... 42% 58.4M 2s
    ##  56650K .......... .......... .......... .......... .......... 42% 61.4M 2s
    ##  56700K .......... .......... .......... .......... .......... 42% 63.2M 2s
    ##  56750K .......... .......... .......... .......... .......... 42% 65.4M 2s
    ##  56800K .......... .......... .......... .......... .......... 42% 75.4M 2s
    ##  56850K .......... .......... .......... .......... .......... 42% 71.1M 2s
    ##  56900K .......... .......... .......... .......... .......... 42% 64.0M 2s
    ##  56950K .......... .......... .......... .......... .......... 42% 62.2M 2s
    ##  57000K .......... .......... .......... .......... .......... 42% 73.1M 2s
    ##  57050K .......... .......... .......... .......... .......... 42% 58.1M 2s
    ##  57100K .......... .......... .......... .......... .......... 42% 63.8M 2s
    ##  57150K .......... .......... .......... .......... .......... 42% 68.2M 2s
    ##  57200K .......... .......... .......... .......... .......... 42% 68.3M 2s
    ##  57250K .......... .......... .......... .......... .......... 42% 68.0M 2s
    ##  57300K .......... .......... .......... .......... .......... 42% 63.7M 2s
    ##  57350K .......... .......... .......... .......... .......... 42% 72.8M 2s
    ##  57400K .......... .......... .......... .......... .......... 42% 89.1M 2s
    ##  57450K .......... .......... .......... .......... .......... 42% 88.7M 2s
    ##  57500K .......... .......... .......... .......... .......... 42% 22.8M 2s
    ##  57550K .......... .......... .......... .......... .......... 42% 61.5M 2s
    ##  57600K .......... .......... .......... .......... .......... 43% 70.9M 2s
    ##  57650K .......... .......... .......... .......... .......... 43% 74.2M 2s
    ##  57700K .......... .......... .......... .......... .......... 43% 60.5M 2s
    ##  57750K .......... .......... .......... .......... .......... 43% 54.5M 2s
    ##  57800K .......... .......... .......... .......... .......... 43% 82.4M 2s
    ##  57850K .......... .......... .......... .......... .......... 43% 60.3M 2s
    ##  57900K .......... .......... .......... .......... .......... 43% 66.1M 2s
    ##  57950K .......... .......... .......... .......... .......... 43% 62.2M 2s
    ##  58000K .......... .......... .......... .......... .......... 43% 57.6M 2s
    ##  58050K .......... .......... .......... .......... .......... 43% 69.8M 2s
    ##  58100K .......... .......... .......... .......... .......... 43% 62.3M 2s
    ##  58150K .......... .......... .......... .......... .......... 43% 63.7M 2s
    ##  58200K .......... .......... .......... .......... .......... 43% 57.7M 2s
    ##  58250K .......... .......... .......... .......... .......... 43% 60.0M 2s
    ##  58300K .......... .......... .......... .......... .......... 43% 68.2M 2s
    ##  58350K .......... .......... .......... .......... .......... 43% 70.5M 2s
    ##  58400K .......... .......... .......... .......... .......... 43% 71.5M 2s
    ##  58450K .......... .......... .......... .......... .......... 43% 78.2M 2s
    ##  58500K .......... .......... .......... .......... .......... 43% 60.3M 2s
    ##  58550K .......... .......... .......... .......... .......... 43% 52.6M 2s
    ##  58600K .......... .......... .......... .......... .......... 43% 62.6M 2s
    ##  58650K .......... .......... .......... .......... .......... 43% 81.3M 2s
    ##  58700K .......... .......... .......... .......... .......... 43% 58.6M 2s
    ##  58750K .......... .......... .......... .......... .......... 43% 60.9M 2s
    ##  58800K .......... .......... .......... .......... .......... 43% 62.7M 2s
    ##  58850K .......... .......... .......... .......... .......... 43% 76.3M 2s
    ##  58900K .......... .......... .......... .......... .......... 43% 74.5M 2s
    ##  58950K .......... .......... .......... .......... .......... 44% 68.2M 2s
    ##  59000K .......... .......... .......... .......... .......... 44% 82.3M 2s
    ##  59050K .......... .......... .......... .......... .......... 44% 68.5M 2s
    ##  59100K .......... .......... .......... .......... .......... 44% 84.6M 2s
    ##  59150K .......... .......... .......... .......... .......... 44% 65.8M 2s
    ##  59200K .......... .......... .......... .......... .......... 44% 66.8M 2s
    ##  59250K .......... .......... .......... .......... .......... 44%  101M 2s
    ##  59300K .......... .......... .......... .......... .......... 44% 84.1M 2s
    ##  59350K .......... .......... .......... .......... .......... 44% 58.0M 2s
    ##  59400K .......... .......... .......... .......... .......... 44% 67.4M 2s
    ##  59450K .......... .......... .......... .......... .......... 44% 64.2M 2s
    ##  59500K .......... .......... .......... .......... .......... 44% 64.1M 2s
    ##  59550K .......... .......... .......... .......... .......... 44% 70.7M 2s
    ##  59600K .......... .......... .......... .......... .......... 44% 69.2M 2s
    ##  59650K .......... .......... .......... .......... .......... 44% 63.4M 2s
    ##  59700K .......... .......... .......... .......... .......... 44% 71.3M 2s
    ##  59750K .......... .......... .......... .......... .......... 44% 67.3M 2s
    ##  59800K .......... .......... .......... .......... .......... 44% 78.8M 2s
    ##  59850K .......... .......... .......... .......... .......... 44% 76.5M 2s
    ##  59900K .......... .......... .......... .......... .......... 44% 63.9M 2s
    ##  59950K .......... .......... .......... .......... .......... 44% 62.8M 2s
    ##  60000K .......... .......... .......... .......... .......... 44% 77.2M 2s
    ##  60050K .......... .......... .......... .......... .......... 44% 58.9M 2s
    ##  60100K .......... .......... .......... .......... .......... 44% 59.4M 2s
    ##  60150K .......... .......... .......... .......... .......... 44% 55.3M 2s
    ##  60200K .......... .......... .......... .......... .......... 44% 66.6M 2s
    ##  60250K .......... .......... .......... .......... .......... 44% 60.7M 2s
    ##  60300K .......... .......... .......... .......... .......... 45% 68.3M 2s
    ##  60350K .......... .......... .......... .......... .......... 45% 53.2M 2s
    ##  60400K .......... .......... .......... .......... .......... 45% 61.4M 2s
    ##  60450K .......... .......... .......... .......... .......... 45% 60.3M 2s
    ##  60500K .......... .......... .......... .......... .......... 45% 58.3M 2s
    ##  60550K .......... .......... .......... .......... .......... 45% 71.9M 2s
    ##  60600K .......... .......... .......... .......... .......... 45% 77.5M 2s
    ##  60650K .......... .......... .......... .......... .......... 45% 75.0M 2s
    ##  60700K .......... .......... .......... .......... .......... 45% 63.0M 2s
    ##  60750K .......... .......... .......... .......... .......... 45% 46.5M 2s
    ##  60800K .......... .......... .......... .......... .......... 45% 84.5M 2s
    ##  60850K .......... .......... .......... .......... .......... 45% 82.5M 2s
    ##  60900K .......... .......... .......... .......... .......... 45% 80.7M 2s
    ##  60950K .......... .......... .......... .......... .......... 45% 68.8M 2s
    ##  61000K .......... .......... .......... .......... .......... 45% 95.5M 2s
    ##  61050K .......... .......... .......... .......... .......... 45% 83.6M 2s
    ##  61100K .......... .......... .......... .......... .......... 45% 71.9M 2s
    ##  61150K .......... .......... .......... .......... .......... 45% 54.3M 2s
    ##  61200K .......... .......... .......... .......... .......... 45% 70.8M 2s
    ##  61250K .......... .......... .......... .......... .......... 45% 75.8M 2s
    ##  61300K .......... .......... .......... .......... .......... 45% 66.3M 2s
    ##  61350K .......... .......... .......... .......... .......... 45% 77.8M 2s
    ##  61400K .......... .......... .......... .......... .......... 45% 83.0M 2s
    ##  61450K .......... .......... .......... .......... .......... 45% 77.4M 2s
    ##  61500K .......... .......... .......... .......... .......... 45%  104M 2s
    ##  61550K .......... .......... .......... .......... .......... 45% 70.1M 2s
    ##  61600K .......... .......... .......... .......... .......... 45% 86.0M 2s
    ##  61650K .......... .......... .......... .......... .......... 46% 77.3M 2s
    ##  61700K .......... .......... .......... .......... .......... 46% 73.1M 2s
    ##  61750K .......... .......... .......... .......... .......... 46% 89.7M 2s
    ##  61800K .......... .......... .......... .......... .......... 46% 79.5M 2s
    ##  61850K .......... .......... .......... .......... .......... 46% 86.9M 2s
    ##  61900K .......... .......... .......... .......... .......... 46% 95.4M 2s
    ##  61950K .......... .......... .......... .......... .......... 46% 70.2M 2s
    ##  62000K .......... .......... .......... .......... .......... 46% 86.7M 2s
    ##  62050K .......... .......... .......... .......... .......... 46% 73.4M 2s
    ##  62100K .......... .......... .......... .......... .......... 46% 83.3M 2s
    ##  62150K .......... .......... .......... .......... .......... 46% 61.9M 2s
    ##  62200K .......... .......... .......... .......... .......... 46% 69.3M 2s
    ##  62250K .......... .......... .......... .......... .......... 46% 92.3M 2s
    ##  62300K .......... .......... .......... .......... .......... 46% 62.5M 2s
    ##  62350K .......... .......... .......... .......... .......... 46% 66.9M 2s
    ##  62400K .......... .......... .......... .......... .......... 46% 81.7M 2s
    ##  62450K .......... .......... .......... .......... .......... 46% 85.7M 2s
    ##  62500K .......... .......... .......... .......... .......... 46% 86.7M 2s
    ##  62550K .......... .......... .......... .......... .......... 46% 23.1M 2s
    ##  62600K .......... .......... .......... .......... .......... 46% 78.5M 2s
    ##  62650K .......... .......... .......... .......... .......... 46% 81.5M 2s
    ##  62700K .......... .......... .......... .......... .......... 46% 92.1M 2s
    ##  62750K .......... .......... .......... .......... .......... 46% 69.0M 2s
    ##  62800K .......... .......... .......... .......... .......... 46% 69.2M 2s
    ##  62850K .......... .......... .......... .......... .......... 46% 76.9M 2s
    ##  62900K .......... .......... .......... .......... .......... 46% 96.8M 2s
    ##  62950K .......... .......... .......... .......... .......... 46% 70.1M 2s
    ##  63000K .......... .......... .......... .......... .......... 47% 72.8M 2s
    ##  63050K .......... .......... .......... .......... .......... 47% 76.1M 2s
    ##  63100K .......... .......... .......... .......... .......... 47% 86.9M 2s
    ##  63150K .......... .......... .......... .......... .......... 47% 71.1M 2s
    ##  63200K .......... .......... .......... .......... .......... 47% 73.6M 2s
    ##  63250K .......... .......... .......... .......... .......... 47% 83.7M 2s
    ##  63300K .......... .......... .......... .......... .......... 47% 75.5M 2s
    ##  63350K .......... .......... .......... .......... .......... 47% 75.0M 2s
    ##  63400K .......... .......... .......... .......... .......... 47% 55.0M 2s
    ##  63450K .......... .......... .......... .......... .......... 47% 27.6M 2s
    ##  63500K .......... .......... .......... .......... .......... 47% 91.9M 2s
    ##  63550K .......... .......... .......... .......... .......... 47% 61.8M 2s
    ##  63600K .......... .......... .......... .......... .......... 47% 73.3M 2s
    ##  63650K .......... .......... .......... .......... .......... 47% 77.0M 2s
    ##  63700K .......... .......... .......... .......... .......... 47% 73.9M 2s
    ##  63750K .......... .......... .......... .......... .......... 47% 59.1M 2s
    ##  63800K .......... .......... .......... .......... .......... 47% 63.4M 2s
    ##  63850K .......... .......... .......... .......... .......... 47% 85.8M 2s
    ##  63900K .......... .......... .......... .......... .......... 47% 69.7M 2s
    ##  63950K .......... .......... .......... .......... .......... 47% 66.3M 2s
    ##  64000K .......... .......... .......... .......... .......... 47% 86.3M 2s
    ##  64050K .......... .......... .......... .......... .......... 47% 70.6M 2s
    ##  64100K .......... .......... .......... .......... .......... 47% 72.1M 2s
    ##  64150K .......... .......... .......... .......... .......... 47% 68.2M 2s
    ##  64200K .......... .......... .......... .......... .......... 47% 68.4M 2s
    ##  64250K .......... .......... .......... .......... .......... 47% 79.9M 2s
    ##  64300K .......... .......... .......... .......... .......... 47% 91.0M 2s
    ##  64350K .......... .......... .......... .......... .......... 48% 62.7M 2s
    ##  64400K .......... .......... .......... .......... .......... 48% 89.1M 2s
    ##  64450K .......... .......... .......... .......... .......... 48% 82.8M 2s
    ##  64500K .......... .......... .......... .......... .......... 48% 82.1M 2s
    ##  64550K .......... .......... .......... .......... .......... 48% 64.7M 2s
    ##  64600K .......... .......... .......... .......... .......... 48% 79.3M 2s
    ##  64650K .......... .......... .......... .......... .......... 48% 75.2M 2s
    ##  64700K .......... .......... .......... .......... .......... 48% 76.9M 2s
    ##  64750K .......... .......... .......... .......... .......... 48% 63.8M 2s
    ##  64800K .......... .......... .......... .......... .......... 48% 66.9M 2s
    ##  64850K .......... .......... .......... .......... .......... 48% 86.3M 2s
    ##  64900K .......... .......... .......... .......... .......... 48% 79.3M 2s
    ##  64950K .......... .......... .......... .......... .......... 48% 73.8M 2s
    ##  65000K .......... .......... .......... .......... .......... 48% 83.8M 2s
    ##  65050K .......... .......... .......... .......... .......... 48% 72.9M 2s
    ##  65100K .......... .......... .......... .......... .......... 48% 71.1M 2s
    ##  65150K .......... .......... .......... .......... .......... 48% 71.8M 2s
    ##  65200K .......... .......... .......... .......... .......... 48% 74.6M 2s
    ##  65250K .......... .......... .......... .......... .......... 48% 76.7M 2s
    ##  65300K .......... .......... .......... .......... .......... 48% 80.3M 2s
    ##  65350K .......... .......... .......... .......... .......... 48% 56.7M 2s
    ##  65400K .......... .......... .......... .......... .......... 48% 69.9M 2s
    ##  65450K .......... .......... .......... .......... .......... 48% 79.6M 2s
    ##  65500K .......... .......... .......... .......... .......... 48% 78.3M 2s
    ##  65550K .......... .......... .......... .......... .......... 48% 56.5M 2s
    ##  65600K .......... .......... .......... .......... .......... 48% 74.6M 2s
    ##  65650K .......... .......... .......... .......... .......... 49% 67.0M 2s
    ##  65700K .......... .......... .......... .......... .......... 49% 78.5M 2s
    ##  65750K .......... .......... .......... .......... .......... 49% 66.8M 2s
    ##  65800K .......... .......... .......... .......... .......... 49% 71.3M 2s
    ##  65850K .......... .......... .......... .......... .......... 49% 75.1M 2s
    ##  65900K .......... .......... .......... .......... .......... 49% 64.8M 1s
    ##  65950K .......... .......... .......... .......... .......... 49% 61.8M 1s
    ##  66000K .......... .......... .......... .......... .......... 49% 84.5M 1s
    ##  66050K .......... .......... .......... .......... .......... 49% 60.6M 1s
    ##  66100K .......... .......... .......... .......... .......... 49% 76.9M 1s
    ##  66150K .......... .......... .......... .......... .......... 49% 71.8M 1s
    ##  66200K .......... .......... .......... .......... .......... 49% 77.5M 1s
    ##  66250K .......... .......... .......... .......... .......... 49% 79.3M 1s
    ##  66300K .......... .......... .......... .......... .......... 49% 71.3M 1s
    ##  66350K .......... .......... .......... .......... .......... 49% 56.2M 1s
    ##  66400K .......... .......... .......... .......... .......... 49% 78.1M 1s
    ##  66450K .......... .......... .......... .......... .......... 49% 72.7M 1s
    ##  66500K .......... .......... .......... .......... .......... 49% 86.3M 1s
    ##  66550K .......... .......... .......... .......... .......... 49% 67.5M 1s
    ##  66600K .......... .......... .......... .......... .......... 49% 74.7M 1s
    ##  66650K .......... .......... .......... .......... .......... 49% 85.0M 1s
    ##  66700K .......... .......... .......... .......... .......... 49% 75.2M 1s
    ##  66750K .......... .......... .......... .......... .......... 49% 64.5M 1s
    ##  66800K .......... .......... .......... .......... .......... 49% 80.7M 1s
    ##  66850K .......... .......... .......... .......... .......... 49% 71.5M 1s
    ##  66900K .......... .......... .......... .......... .......... 49% 76.5M 1s
    ##  66950K .......... .......... .......... .......... .......... 49% 70.0M 1s
    ##  67000K .......... .......... .......... .......... .......... 50% 69.1M 1s
    ##  67050K .......... .......... .......... .......... .......... 50% 88.1M 1s
    ##  67100K .......... .......... .......... .......... .......... 50% 87.3M 1s
    ##  67150K .......... .......... .......... .......... .......... 50% 57.1M 1s
    ##  67200K .......... .......... .......... .......... .......... 50% 76.6M 1s
    ##  67250K .......... .......... .......... .......... .......... 50% 71.8M 1s
    ##  67300K .......... .......... .......... .......... .......... 50% 94.6M 1s
    ##  67350K .......... .......... .......... .......... .......... 50% 79.7M 1s
    ##  67400K .......... .......... .......... .......... .......... 50% 75.2M 1s
    ##  67450K .......... .......... .......... .......... .......... 50% 74.7M 1s
    ##  67500K .......... .......... .......... .......... .......... 50% 72.2M 1s
    ##  67550K .......... .......... .......... .......... .......... 50% 71.1M 1s
    ##  67600K .......... .......... .......... .......... .......... 50% 88.6M 1s
    ##  67650K .......... .......... .......... .......... .......... 50% 86.4M 1s
    ##  67700K .......... .......... .......... .......... .......... 50% 82.2M 1s
    ##  67750K .......... .......... .......... .......... .......... 50% 77.3M 1s
    ##  67800K .......... .......... .......... .......... .......... 50% 65.9M 1s
    ##  67850K .......... .......... .......... .......... .......... 50% 69.2M 1s
    ##  67900K .......... .......... .......... .......... .......... 50% 82.6M 1s
    ##  67950K .......... .......... .......... .......... .......... 50% 70.0M 1s
    ##  68000K .......... .......... .......... .......... .......... 50% 76.0M 1s
    ##  68050K .......... .......... .......... .......... .......... 50% 67.2M 1s
    ##  68100K .......... .......... .......... .......... .......... 50% 77.7M 1s
    ##  68150K .......... .......... .......... .......... .......... 50% 78.3M 1s
    ##  68200K .......... .......... .......... .......... .......... 50% 78.6M 1s
    ##  68250K .......... .......... .......... .......... .......... 50% 69.1M 1s
    ##  68300K .......... .......... .......... .......... .......... 50% 75.1M 1s
    ##  68350K .......... .......... .......... .......... .......... 51% 44.9M 1s
    ##  68400K .......... .......... .......... .......... .......... 51% 75.7M 1s
    ##  68450K .......... .......... .......... .......... .......... 51% 81.5M 1s
    ##  68500K .......... .......... .......... .......... .......... 51% 67.8M 1s
    ##  68550K .......... .......... .......... .......... .......... 51% 58.2M 1s
    ##  68600K .......... .......... .......... .......... .......... 51% 57.7M 1s
    ##  68650K .......... .......... .......... .......... .......... 51% 74.7M 1s
    ##  68700K .......... .......... .......... .......... .......... 51% 69.7M 1s
    ##  68750K .......... .......... .......... .......... .......... 51% 45.7M 1s
    ##  68800K .......... .......... .......... .......... .......... 51% 35.6M 1s
    ##  68850K .......... .......... .......... .......... .......... 51% 90.8M 1s
    ##  68900K .......... .......... .......... .......... .......... 51% 76.2M 1s
    ##  68950K .......... .......... .......... .......... .......... 51% 51.5M 1s
    ##  69000K .......... .......... .......... .......... .......... 51% 63.5M 1s
    ##  69050K .......... .......... .......... .......... .......... 51% 58.4M 1s
    ##  69100K .......... .......... .......... .......... .......... 51% 71.5M 1s
    ##  69150K .......... .......... .......... .......... .......... 51% 68.6M 1s
    ##  69200K .......... .......... .......... .......... .......... 51% 70.5M 1s
    ##  69250K .......... .......... .......... .......... .......... 51% 97.8M 1s
    ##  69300K .......... .......... .......... .......... .......... 51% 99.5M 1s
    ##  69350K .......... .......... .......... .......... .......... 51% 56.0M 1s
    ##  69400K .......... .......... .......... .......... .......... 51% 77.2M 1s
    ##  69450K .......... .......... .......... .......... .......... 51% 60.2M 1s
    ##  69500K .......... .......... .......... .......... .......... 51% 65.6M 1s
    ##  69550K .......... .......... .......... .......... .......... 51% 61.1M 1s
    ##  69600K .......... .......... .......... .......... .......... 51% 86.0M 1s
    ##  69650K .......... .......... .......... .......... .......... 51% 61.4M 1s
    ##  69700K .......... .......... .......... .......... .......... 52% 84.7M 1s
    ##  69750K .......... .......... .......... .......... .......... 52% 60.8M 1s
    ##  69800K .......... .......... .......... .......... .......... 52% 79.2M 1s
    ##  69850K .......... .......... .......... .......... .......... 52% 72.4M 1s
    ##  69900K .......... .......... .......... .......... .......... 52% 58.9M 1s
    ##  69950K .......... .......... .......... .......... .......... 52% 50.2M 1s
    ##  70000K .......... .......... .......... .......... .......... 52% 81.1M 1s
    ##  70050K .......... .......... .......... .......... .......... 52%  114M 1s
    ##  70100K .......... .......... .......... .......... .......... 52% 93.5M 1s
    ##  70150K .......... .......... .......... .......... .......... 52% 77.2M 1s
    ##  70200K .......... .......... .......... .......... .......... 52% 60.4M 1s
    ##  70250K .......... .......... .......... .......... .......... 52% 57.9M 1s
    ##  70300K .......... .......... .......... .......... .......... 52% 68.7M 1s
    ##  70350K .......... .......... .......... .......... .......... 52% 51.8M 1s
    ##  70400K .......... .......... .......... .......... .......... 52% 75.2M 1s
    ##  70450K .......... .......... .......... .......... .......... 52% 68.0M 1s
    ##  70500K .......... .......... .......... .......... .......... 52% 68.3M 1s
    ##  70550K .......... .......... .......... .......... .......... 52% 62.0M 1s
    ##  70600K .......... .......... .......... .......... .......... 52% 70.0M 1s
    ##  70650K .......... .......... .......... .......... .......... 52% 59.2M 1s
    ##  70700K .......... .......... .......... .......... .......... 52% 66.1M 1s
    ##  70750K .......... .......... .......... .......... .......... 52% 63.3M 1s
    ##  70800K .......... .......... .......... .......... .......... 52% 53.4M 1s
    ##  70850K .......... .......... .......... .......... .......... 52% 70.6M 1s
    ##  70900K .......... .......... .......... .......... .......... 52% 63.4M 1s
    ##  70950K .......... .......... .......... .......... .......... 52% 71.1M 1s
    ##  71000K .......... .......... .......... .......... .......... 52% 85.7M 1s
    ##  71050K .......... .......... .......... .......... .......... 53% 73.7M 1s
    ##  71100K .......... .......... .......... .......... .......... 53% 72.9M 1s
    ##  71150K .......... .......... .......... .......... .......... 53% 67.4M 1s
    ##  71200K .......... .......... .......... .......... .......... 53% 77.4M 1s
    ##  71250K .......... .......... .......... .......... .......... 53% 57.8M 1s
    ##  71300K .......... .......... .......... .......... .......... 53% 70.8M 1s
    ##  71350K .......... .......... .......... .......... .......... 53% 61.7M 1s
    ##  71400K .......... .......... .......... .......... .......... 53% 59.2M 1s
    ##  71450K .......... .......... .......... .......... .......... 53% 77.0M 1s
    ##  71500K .......... .......... .......... .......... .......... 53% 70.7M 1s
    ##  71550K .......... .......... .......... .......... .......... 53% 58.0M 1s
    ##  71600K .......... .......... .......... .......... .......... 53% 69.9M 1s
    ##  71650K .......... .......... .......... .......... .......... 53% 64.7M 1s
    ##  71700K .......... .......... .......... .......... .......... 53% 63.5M 1s
    ##  71750K .......... .......... .......... .......... .......... 53% 56.4M 1s
    ##  71800K .......... .......... .......... .......... .......... 53% 73.6M 1s
    ##  71850K .......... .......... .......... .......... .......... 53% 73.0M 1s
    ##  71900K .......... .......... .......... .......... .......... 53% 79.8M 1s
    ##  71950K .......... .......... .......... .......... .......... 53% 50.5M 1s
    ##  72000K .......... .......... .......... .......... .......... 53% 81.0M 1s
    ##  72050K .......... .......... .......... .......... .......... 53% 69.7M 1s
    ##  72100K .......... .......... .......... .......... .......... 53% 85.9M 1s
    ##  72150K .......... .......... .......... .......... .......... 53% 60.4M 1s
    ##  72200K .......... .......... .......... .......... .......... 53% 71.7M 1s
    ##  72250K .......... .......... .......... .......... .......... 53% 79.0M 1s
    ##  72300K .......... .......... .......... .......... .......... 53% 64.2M 1s
    ##  72350K .......... .......... .......... .......... .......... 54% 65.0M 1s
    ##  72400K .......... .......... .......... .......... .......... 54% 69.7M 1s
    ##  72450K .......... .......... .......... .......... .......... 54% 68.8M 1s
    ##  72500K .......... .......... .......... .......... .......... 54% 63.1M 1s
    ##  72550K .......... .......... .......... .......... .......... 54% 61.2M 1s
    ##  72600K .......... .......... .......... .......... .......... 54% 78.7M 1s
    ##  72650K .......... .......... .......... .......... .......... 54% 73.9M 1s
    ##  72700K .......... .......... .......... .......... .......... 54% 73.1M 1s
    ##  72750K .......... .......... .......... .......... .......... 54% 55.6M 1s
    ##  72800K .......... .......... .......... .......... .......... 54% 69.6M 1s
    ##  72850K .......... .......... .......... .......... .......... 54% 58.7M 1s
    ##  72900K .......... .......... .......... .......... .......... 54% 75.8M 1s
    ##  72950K .......... .......... .......... .......... .......... 54% 57.0M 1s
    ##  73000K .......... .......... .......... .......... .......... 54% 78.4M 1s
    ##  73050K .......... .......... .......... .......... .......... 54% 66.6M 1s
    ##  73100K .......... .......... .......... .......... .......... 54% 64.0M 1s
    ##  73150K .......... .......... .......... .......... .......... 54% 60.5M 1s
    ##  73200K .......... .......... .......... .......... .......... 54% 89.2M 1s
    ##  73250K .......... .......... .......... .......... .......... 54% 66.8M 1s
    ##  73300K .......... .......... .......... .......... .......... 54% 59.4M 1s
    ##  73350K .......... .......... .......... .......... .......... 54% 51.7M 1s
    ##  73400K .......... .......... .......... .......... .......... 54% 68.0M 1s
    ##  73450K .......... .......... .......... .......... .......... 54% 65.7M 1s
    ##  73500K .......... .......... .......... .......... .......... 54% 72.3M 1s
    ##  73550K .......... .......... .......... .......... .......... 54% 55.2M 1s
    ##  73600K .......... .......... .......... .......... .......... 54% 64.9M 1s
    ##  73650K .......... .......... .......... .......... .......... 54% 65.3M 1s
    ##  73700K .......... .......... .......... .......... .......... 55% 77.5M 1s
    ##  73750K .......... .......... .......... .......... .......... 55% 61.5M 1s
    ##  73800K .......... .......... .......... .......... .......... 55% 60.8M 1s
    ##  73850K .......... .......... .......... .......... .......... 55% 73.9M 1s
    ##  73900K .......... .......... .......... .......... .......... 55% 64.5M 1s
    ##  73950K .......... .......... .......... .......... .......... 55% 69.4M 1s
    ##  74000K .......... .......... .......... .......... .......... 55% 77.2M 1s
    ##  74050K .......... .......... .......... .......... .......... 55% 82.7M 1s
    ##  74100K .......... .......... .......... .......... .......... 55% 68.6M 1s
    ##  74150K .......... .......... .......... .......... .......... 55% 65.7M 1s
    ##  74200K .......... .......... .......... .......... .......... 55% 67.4M 1s
    ##  74250K .......... .......... .......... .......... .......... 55% 65.7M 1s
    ##  74300K .......... .......... .......... .......... .......... 55% 80.4M 1s
    ##  74350K .......... .......... .......... .......... .......... 55% 56.3M 1s
    ##  74400K .......... .......... .......... .......... .......... 55% 66.1M 1s
    ##  74450K .......... .......... .......... .......... .......... 55% 99.1M 1s
    ##  74500K .......... .......... .......... .......... .......... 55% 58.5M 1s
    ##  74550K .......... .......... .......... .......... .......... 55% 60.8M 1s
    ##  74600K .......... .......... .......... .......... .......... 55% 64.4M 1s
    ##  74650K .......... .......... .......... .......... .......... 55% 73.1M 1s
    ##  74700K .......... .......... .......... .......... .......... 55% 58.2M 1s
    ##  74750K .......... .......... .......... .......... .......... 55% 53.4M 1s
    ##  74800K .......... .......... .......... .......... .......... 55% 78.4M 1s
    ##  74850K .......... .......... .......... .......... .......... 55% 77.4M 1s
    ##  74900K .......... .......... .......... .......... .......... 55% 65.0M 1s
    ##  74950K .......... .......... .......... .......... .......... 55% 65.7M 1s
    ##  75000K .......... .......... .......... .......... .......... 55% 79.5M 1s
    ##  75050K .......... .......... .......... .......... .......... 56% 73.6M 1s
    ##  75100K .......... .......... .......... .......... .......... 56% 64.2M 1s
    ##  75150K .......... .......... .......... .......... .......... 56% 61.8M 1s
    ##  75200K .......... .......... .......... .......... .......... 56% 75.0M 1s
    ##  75250K .......... .......... .......... .......... .......... 56% 57.5M 1s
    ##  75300K .......... .......... .......... .......... .......... 56% 69.1M 1s
    ##  75350K .......... .......... .......... .......... .......... 56% 66.8M 1s
    ##  75400K .......... .......... .......... .......... .......... 56% 73.5M 1s
    ##  75450K .......... .......... .......... .......... .......... 56% 70.3M 1s
    ##  75500K .......... .......... .......... .......... .......... 56% 81.6M 1s
    ##  75550K .......... .......... .......... .......... .......... 56% 61.8M 1s
    ##  75600K .......... .......... .......... .......... .......... 56% 70.3M 1s
    ##  75650K .......... .......... .......... .......... .......... 56% 57.5M 1s
    ##  75700K .......... .......... .......... .......... .......... 56% 69.4M 1s
    ##  75750K .......... .......... .......... .......... .......... 56% 62.6M 1s
    ##  75800K .......... .......... .......... .......... .......... 56% 67.1M 1s
    ##  75850K .......... .......... .......... .......... .......... 56% 58.1M 1s
    ##  75900K .......... .......... .......... .......... .......... 56% 56.6M 1s
    ##  75950K .......... .......... .......... .......... .......... 56% 65.4M 1s
    ##  76000K .......... .......... .......... .......... .......... 56% 55.6M 1s
    ##  76050K .......... .......... .......... .......... .......... 56% 54.8M 1s
    ##  76100K .......... .......... .......... .......... .......... 56% 68.6M 1s
    ##  76150K .......... .......... .......... .......... .......... 56% 56.7M 1s
    ##  76200K .......... .......... .......... .......... .......... 56% 72.6M 1s
    ##  76250K .......... .......... .......... .......... .......... 56% 64.4M 1s
    ##  76300K .......... .......... .......... .......... .......... 56% 61.4M 1s
    ##  76350K .......... .......... .......... .......... .......... 56% 50.4M 1s
    ##  76400K .......... .......... .......... .......... .......... 57% 60.7M 1s
    ##  76450K .......... .......... .......... .......... .......... 57% 75.1M 1s
    ##  76500K .......... .......... .......... .......... .......... 57% 74.1M 1s
    ##  76550K .......... .......... .......... .......... .......... 57% 65.1M 1s
    ##  76600K .......... .......... .......... .......... .......... 57% 68.4M 1s
    ##  76650K .......... .......... .......... .......... .......... 57%  101M 1s
    ##  76700K .......... .......... .......... .......... .......... 57% 81.5M 1s
    ##  76750K .......... .......... .......... .......... .......... 57% 63.2M 1s
    ##  76800K .......... .......... .......... .......... .......... 57% 66.4M 1s
    ##  76850K .......... .......... .......... .......... .......... 57% 69.5M 1s
    ##  76900K .......... .......... .......... .......... .......... 57% 59.2M 1s
    ##  76950K .......... .......... .......... .......... .......... 57% 72.4M 1s
    ##  77000K .......... .......... .......... .......... .......... 57% 73.6M 1s
    ##  77050K .......... .......... .......... .......... .......... 57% 61.0M 1s
    ##  77100K .......... .......... .......... .......... .......... 57% 82.0M 1s
    ##  77150K .......... .......... .......... .......... .......... 57% 52.5M 1s
    ##  77200K .......... .......... .......... .......... .......... 57% 78.8M 1s
    ##  77250K .......... .......... .......... .......... .......... 57% 92.1M 1s
    ##  77300K .......... .......... .......... .......... .......... 57% 59.2M 1s
    ##  77350K .......... .......... .......... .......... .......... 57% 58.9M 1s
    ##  77400K .......... .......... .......... .......... .......... 57% 71.4M 1s
    ##  77450K .......... .......... .......... .......... .......... 57% 79.1M 1s
    ##  77500K .......... .......... .......... .......... .......... 57% 70.1M 1s
    ##  77550K .......... .......... .......... .......... .......... 57% 62.9M 1s
    ##  77600K .......... .......... .......... .......... .......... 57% 59.6M 1s
    ##  77650K .......... .......... .......... .......... .......... 57% 77.0M 1s
    ##  77700K .......... .......... .......... .......... .......... 57% 71.6M 1s
    ##  77750K .......... .......... .......... .......... .......... 58% 54.7M 1s
    ##  77800K .......... .......... .......... .......... .......... 58% 60.1M 1s
    ##  77850K .......... .......... .......... .......... .......... 58% 71.0M 1s
    ##  77900K .......... .......... .......... .......... .......... 58% 64.8M 1s
    ##  77950K .......... .......... .......... .......... .......... 58% 61.1M 1s
    ##  78000K .......... .......... .......... .......... .......... 58% 68.8M 1s
    ##  78050K .......... .......... .......... .......... .......... 58% 64.5M 1s
    ##  78100K .......... .......... .......... .......... .......... 58% 63.2M 1s
    ##  78150K .......... .......... .......... .......... .......... 58% 59.4M 1s
    ##  78200K .......... .......... .......... .......... .......... 58% 69.5M 1s
    ##  78250K .......... .......... .......... .......... .......... 58% 78.4M 1s
    ##  78300K .......... .......... .......... .......... .......... 58% 66.5M 1s
    ##  78350K .......... .......... .......... .......... .......... 58% 62.8M 1s
    ##  78400K .......... .......... .......... .......... .......... 58% 81.0M 1s
    ##  78450K .......... .......... .......... .......... .......... 58% 73.9M 1s
    ##  78500K .......... .......... .......... .......... .......... 58% 60.5M 1s
    ##  78550K .......... .......... .......... .......... .......... 58% 52.1M 1s
    ##  78600K .......... .......... .......... .......... .......... 58% 57.3M 1s
    ##  78650K .......... .......... .......... .......... .......... 58% 57.6M 1s
    ##  78700K .......... .......... .......... .......... .......... 58% 73.6M 1s
    ##  78750K .......... .......... .......... .......... .......... 58% 66.2M 1s
    ##  78800K .......... .......... .......... .......... .......... 58% 72.6M 1s
    ##  78850K .......... .......... .......... .......... .......... 58% 64.5M 1s
    ##  78900K .......... .......... .......... .......... .......... 58% 95.4M 1s
    ##  78950K .......... .......... .......... .......... .......... 58% 56.5M 1s
    ##  79000K .......... .......... .......... .......... .......... 58% 72.3M 1s
    ##  79050K .......... .......... .......... .......... .......... 59% 68.7M 1s
    ##  79100K .......... .......... .......... .......... .......... 59% 78.0M 1s
    ##  79150K .......... .......... .......... .......... .......... 59% 65.6M 1s
    ##  79200K .......... .......... .......... .......... .......... 59% 69.6M 1s
    ##  79250K .......... .......... .......... .......... .......... 59% 70.4M 1s
    ##  79300K .......... .......... .......... .......... .......... 59% 75.0M 1s
    ##  79350K .......... .......... .......... .......... .......... 59% 63.4M 1s
    ##  79400K .......... .......... .......... .......... .......... 59% 63.0M 1s
    ##  79450K .......... .......... .......... .......... .......... 59% 64.0M 1s
    ##  79500K .......... .......... .......... .......... .......... 59% 61.9M 1s
    ##  79550K .......... .......... .......... .......... .......... 59% 57.9M 1s
    ##  79600K .......... .......... .......... .......... .......... 59% 62.9M 1s
    ##  79650K .......... .......... .......... .......... .......... 59% 55.2M 1s
    ##  79700K .......... .......... .......... .......... .......... 59% 76.1M 1s
    ##  79750K .......... .......... .......... .......... .......... 59% 83.9M 1s
    ##  79800K .......... .......... .......... .......... .......... 59% 67.8M 1s
    ##  79850K .......... .......... .......... .......... .......... 59% 56.6M 1s
    ##  79900K .......... .......... .......... .......... .......... 59% 74.5M 1s
    ##  79950K .......... .......... .......... .......... .......... 59% 51.5M 1s
    ##  80000K .......... .......... .......... .......... .......... 59% 75.2M 1s
    ##  80050K .......... .......... .......... .......... .......... 59% 86.1M 1s
    ##  80100K .......... .......... .......... .......... .......... 59% 69.9M 1s
    ##  80150K .......... .......... .......... .......... .......... 59% 49.8M 1s
    ##  80200K .......... .......... .......... .......... .......... 59% 78.8M 1s
    ##  80250K .......... .......... .......... .......... .......... 59% 77.3M 1s
    ##  80300K .......... .......... .......... .......... .......... 59% 94.6M 1s
    ##  80350K .......... .......... .......... .......... .......... 59% 49.5M 1s
    ##  80400K .......... .......... .......... .......... .......... 60% 62.1M 1s
    ##  80450K .......... .......... .......... .......... .......... 60% 74.9M 1s
    ##  80500K .......... .......... .......... .......... .......... 60% 75.3M 1s
    ##  80550K .......... .......... .......... .......... .......... 60% 63.9M 1s
    ##  80600K .......... .......... .......... .......... .......... 60% 55.7M 1s
    ##  80650K .......... .......... .......... .......... .......... 60% 55.3M 1s
    ##  80700K .......... .......... .......... .......... .......... 60% 63.1M 1s
    ##  80750K .......... .......... .......... .......... .......... 60% 61.9M 1s
    ##  80800K .......... .......... .......... .......... .......... 60% 67.9M 1s
    ##  80850K .......... .......... .......... .......... .......... 60% 57.9M 1s
    ##  80900K .......... .......... .......... .......... .......... 60% 67.1M 1s
    ##  80950K .......... .......... .......... .......... .......... 60% 61.3M 1s
    ##  81000K .......... .......... .......... .......... .......... 60% 58.7M 1s
    ##  81050K .......... .......... .......... .......... .......... 60% 86.0M 1s
    ##  81100K .......... .......... .......... .......... .......... 60% 62.3M 1s
    ##  81150K .......... .......... .......... .......... .......... 60% 54.1M 1s
    ##  81200K .......... .......... .......... .......... .......... 60% 59.9M 1s
    ##  81250K .......... .......... .......... .......... .......... 60% 65.5M 1s
    ##  81300K .......... .......... .......... .......... .......... 60% 67.2M 1s
    ##  81350K .......... .......... .......... .......... .......... 60% 58.3M 1s
    ##  81400K .......... .......... .......... .......... .......... 60% 75.0M 1s
    ##  81450K .......... .......... .......... .......... .......... 60% 66.4M 1s
    ##  81500K .......... .......... .......... .......... .......... 60% 62.9M 1s
    ##  81550K .......... .......... .......... .......... .......... 60% 49.6M 1s
    ##  81600K .......... .......... .......... .......... .......... 60% 60.7M 1s
    ##  81650K .......... .......... .......... .......... .......... 60% 89.2M 1s
    ##  81700K .......... .......... .......... .......... .......... 60% 97.4M 1s
    ##  81750K .......... .......... .......... .......... .......... 61% 59.1M 1s
    ##  81800K .......... .......... .......... .......... .......... 61% 73.9M 1s
    ##  81850K .......... .......... .......... .......... .......... 61% 98.1M 1s
    ##  81900K .......... .......... .......... .......... .......... 61% 79.4M 1s
    ##  81950K .......... .......... .......... .......... .......... 61% 63.3M 1s
    ##  82000K .......... .......... .......... .......... .......... 61% 67.7M 1s
    ##  82050K .......... .......... .......... .......... .......... 61% 58.1M 1s
    ##  82100K .......... .......... .......... .......... .......... 61% 70.8M 1s
    ##  82150K .......... .......... .......... .......... .......... 61% 83.8M 1s
    ##  82200K .......... .......... .......... .......... .......... 61% 71.7M 1s
    ##  82250K .......... .......... .......... .......... .......... 61% 72.7M 1s
    ##  82300K .......... .......... .......... .......... .......... 61% 76.9M 1s
    ##  82350K .......... .......... .......... .......... .......... 61% 48.8M 1s
    ##  82400K .......... .......... .......... .......... .......... 61% 56.4M 1s
    ##  82450K .......... .......... .......... .......... .......... 61% 71.1M 1s
    ##  82500K .......... .......... .......... .......... .......... 61% 68.6M 1s
    ##  82550K .......... .......... .......... .......... .......... 61% 57.6M 1s
    ##  82600K .......... .......... .......... .......... .......... 61% 66.9M 1s
    ##  82650K .......... .......... .......... .......... .......... 61% 67.4M 1s
    ##  82700K .......... .......... .......... .......... .......... 61% 58.1M 1s
    ##  82750K .......... .......... .......... .......... .......... 61% 58.2M 1s
    ##  82800K .......... .......... .......... .......... .......... 61% 77.8M 1s
    ##  82850K .......... .......... .......... .......... .......... 61% 64.1M 1s
    ##  82900K .......... .......... .......... .......... .......... 61% 98.1M 1s
    ##  82950K .......... .......... .......... .......... .......... 61% 64.6M 1s
    ##  83000K .......... .......... .......... .......... .......... 61% 95.9M 1s
    ##  83050K .......... .......... .......... .......... .......... 61% 83.7M 1s
    ##  83100K .......... .......... .......... .......... .......... 62% 63.3M 1s
    ##  83150K .......... .......... .......... .......... .......... 62% 64.1M 1s
    ##  83200K .......... .......... .......... .......... .......... 62% 72.6M 1s
    ##  83250K .......... .......... .......... .......... .......... 62% 77.7M 1s
    ##  83300K .......... .......... .......... .......... .......... 62% 72.8M 1s
    ##  83350K .......... .......... .......... .......... .......... 62% 53.3M 1s
    ##  83400K .......... .......... .......... .......... .......... 62% 73.6M 1s
    ##  83450K .......... .......... .......... .......... .......... 62% 76.1M 1s
    ##  83500K .......... .......... .......... .......... .......... 62% 68.1M 1s
    ##  83550K .......... .......... .......... .......... .......... 62% 69.5M 1s
    ##  83600K .......... .......... .......... .......... .......... 62% 75.4M 1s
    ##  83650K .......... .......... .......... .......... .......... 62% 71.3M 1s
    ##  83700K .......... .......... .......... .......... .......... 62% 29.8M 1s
    ##  83750K .......... .......... .......... .......... .......... 62% 56.3M 1s
    ##  83800K .......... .......... .......... .......... .......... 62% 26.0M 1s
    ##  83850K .......... .......... .......... .......... .......... 62% 65.7M 1s
    ##  83900K .......... .......... .......... .......... .......... 62% 41.3M 1s
    ##  83950K .......... .......... .......... .......... .......... 62% 55.1M 1s
    ##  84000K .......... .......... .......... .......... .......... 62% 49.9M 1s
    ##  84050K .......... .......... .......... .......... .......... 62% 30.5M 1s
    ##  84100K .......... .......... .......... .......... .......... 62% 53.9M 1s
    ##  84150K .......... .......... .......... .......... .......... 62% 47.3M 1s
    ##  84200K .......... .......... .......... .......... .......... 62% 52.7M 1s
    ##  84250K .......... .......... .......... .......... .......... 62% 38.3M 1s
    ##  84300K .......... .......... .......... .......... .......... 62% 59.3M 1s
    ##  84350K .......... .......... .......... .......... .......... 62% 24.6M 1s
    ##  84400K .......... .......... .......... .......... .......... 62% 38.5M 1s
    ##  84450K .......... .......... .......... .......... .......... 63% 66.7M 1s
    ##  84500K .......... .......... .......... .......... .......... 63% 59.3M 1s
    ##  84550K .......... .......... .......... .......... .......... 63% 38.2M 1s
    ##  84600K .......... .......... .......... .......... .......... 63% 25.4M 1s
    ##  84650K .......... .......... .......... .......... .......... 63% 63.3M 1s
    ##  84700K .......... .......... .......... .......... .......... 63% 69.3M 1s
    ##  84750K .......... .......... .......... .......... .......... 63% 31.8M 1s
    ##  84800K .......... .......... .......... .......... .......... 63% 80.3M 1s
    ##  84850K .......... .......... .......... .......... .......... 63% 63.0M 1s
    ##  84900K .......... .......... .......... .......... .......... 63% 53.7M 1s
    ##  84950K .......... .......... .......... .......... .......... 63% 64.9M 1s
    ##  85000K .......... .......... .......... .......... .......... 63% 60.0M 1s
    ##  85050K .......... .......... .......... .......... .......... 63% 73.4M 1s
    ##  85100K .......... .......... .......... .......... .......... 63% 66.6M 1s
    ##  85150K .......... .......... .......... .......... .......... 63% 62.7M 1s
    ##  85200K .......... .......... .......... .......... .......... 63% 72.9M 1s
    ##  85250K .......... .......... .......... .......... .......... 63% 75.1M 1s
    ##  85300K .......... .......... .......... .......... .......... 63% 82.6M 1s
    ##  85350K .......... .......... .......... .......... .......... 63% 66.8M 1s
    ##  85400K .......... .......... .......... .......... .......... 63% 76.2M 1s
    ##  85450K .......... .......... .......... .......... .......... 63% 59.7M 1s
    ##  85500K .......... .......... .......... .......... .......... 63% 77.6M 1s
    ##  85550K .......... .......... .......... .......... .......... 63% 59.2M 1s
    ##  85600K .......... .......... .......... .......... .......... 63% 66.2M 1s
    ##  85650K .......... .......... .......... .......... .......... 63% 76.4M 1s
    ##  85700K .......... .......... .......... .......... .......... 63% 72.5M 1s
    ##  85750K .......... .......... .......... .......... .......... 63% 53.4M 1s
    ##  85800K .......... .......... .......... .......... .......... 64% 71.7M 1s
    ##  85850K .......... .......... .......... .......... .......... 64% 83.2M 1s
    ##  85900K .......... .......... .......... .......... .......... 64% 63.9M 1s
    ##  85950K .......... .......... .......... .......... .......... 64% 64.2M 1s
    ##  86000K .......... .......... .......... .......... .......... 64% 66.5M 1s
    ##  86050K .......... .......... .......... .......... .......... 64% 59.5M 1s
    ##  86100K .......... .......... .......... .......... .......... 64% 70.8M 1s
    ##  86150K .......... .......... .......... .......... .......... 64% 52.7M 1s
    ##  86200K .......... .......... .......... .......... .......... 64% 59.1M 1s
    ##  86250K .......... .......... .......... .......... .......... 64% 61.3M 1s
    ##  86300K .......... .......... .......... .......... .......... 64% 68.2M 1s
    ##  86350K .......... .......... .......... .......... .......... 64% 65.8M 1s
    ##  86400K .......... .......... .......... .......... .......... 64% 70.9M 1s
    ##  86450K .......... .......... .......... .......... .......... 64% 65.8M 1s
    ##  86500K .......... .......... .......... .......... .......... 64% 62.2M 1s
    ##  86550K .......... .......... .......... .......... .......... 64% 83.0M 1s
    ##  86600K .......... .......... .......... .......... .......... 64%  111M 1s
    ##  86650K .......... .......... .......... .......... .......... 64% 80.8M 1s
    ##  86700K .......... .......... .......... .......... .......... 64% 79.8M 1s
    ##  86750K .......... .......... .......... .......... .......... 64% 68.8M 1s
    ##  86800K .......... .......... .......... .......... .......... 64% 84.1M 1s
    ##  86850K .......... .......... .......... .......... .......... 64% 93.6M 1s
    ##  86900K .......... .......... .......... .......... .......... 64% 98.8M 1s
    ##  86950K .......... .......... .......... .......... .......... 64% 60.6M 1s
    ##  87000K .......... .......... .......... .......... .......... 64% 61.8M 1s
    ##  87050K .......... .......... .......... .......... .......... 64% 71.9M 1s
    ##  87100K .......... .......... .......... .......... .......... 65% 70.5M 1s
    ##  87150K .......... .......... .......... .......... .......... 65% 64.9M 1s
    ##  87200K .......... .......... .......... .......... .......... 65% 76.0M 1s
    ##  87250K .......... .......... .......... .......... .......... 65% 49.0M 1s
    ##  87300K .......... .......... .......... .......... .......... 65% 67.6M 1s
    ##  87350K .......... .......... .......... .......... .......... 65% 67.1M 1s
    ##  87400K .......... .......... .......... .......... .......... 65% 67.3M 1s
    ##  87450K .......... .......... .......... .......... .......... 65% 53.8M 1s
    ##  87500K .......... .......... .......... .......... .......... 65% 63.0M 1s
    ##  87550K .......... .......... .......... .......... .......... 65% 53.8M 1s
    ##  87600K .......... .......... .......... .......... .......... 65% 72.6M 1s
    ##  87650K .......... .......... .......... .......... .......... 65% 62.8M 1s
    ##  87700K .......... .......... .......... .......... .......... 65% 71.5M 1s
    ##  87750K .......... .......... .......... .......... .......... 65% 48.9M 1s
    ##  87800K .......... .......... .......... .......... .......... 65% 55.0M 1s
    ##  87850K .......... .......... .......... .......... .......... 65% 57.9M 1s
    ##  87900K .......... .......... .......... .......... .......... 65% 62.0M 1s
    ##  87950K .......... .......... .......... .......... .......... 65% 63.7M 1s
    ##  88000K .......... .......... .......... .......... .......... 65% 87.1M 1s
    ##  88050K .......... .......... .......... .......... .......... 65% 68.6M 1s
    ##  88100K .......... .......... .......... .......... .......... 65% 72.0M 1s
    ##  88150K .......... .......... .......... .......... .......... 65% 58.3M 1s
    ##  88200K .......... .......... .......... .......... .......... 65% 65.2M 1s
    ##  88250K .......... .......... .......... .......... .......... 65% 55.4M 1s
    ##  88300K .......... .......... .......... .......... .......... 65% 67.4M 1s
    ##  88350K .......... .......... .......... .......... .......... 65% 59.5M 1s
    ##  88400K .......... .......... .......... .......... .......... 65% 74.0M 1s
    ##  88450K .......... .......... .......... .......... .......... 66% 57.6M 1s
    ##  88500K .......... .......... .......... .......... .......... 66% 68.7M 1s
    ##  88550K .......... .......... .......... .......... .......... 66% 61.5M 1s
    ##  88600K .......... .......... .......... .......... .......... 66% 74.2M 1s
    ##  88650K .......... .......... .......... .......... .......... 66% 69.5M 1s
    ##  88700K .......... .......... .......... .......... .......... 66% 58.0M 1s
    ##  88750K .......... .......... .......... .......... .......... 66% 56.3M 1s
    ##  88800K .......... .......... .......... .......... .......... 66% 80.0M 1s
    ##  88850K .......... .......... .......... .......... .......... 66% 68.5M 1s
    ##  88900K .......... .......... .......... .......... .......... 66% 58.3M 1s
    ##  88950K .......... .......... .......... .......... .......... 66% 52.5M 1s
    ##  89000K .......... .......... .......... .......... .......... 66% 69.8M 1s
    ##  89050K .......... .......... .......... .......... .......... 66% 62.6M 1s
    ##  89100K .......... .......... .......... .......... .......... 66% 59.8M 1s
    ##  89150K .......... .......... .......... .......... .......... 66% 50.1M 1s
    ##  89200K .......... .......... .......... .......... .......... 66% 66.7M 1s
    ##  89250K .......... .......... .......... .......... .......... 66% 76.0M 1s
    ##  89300K .......... .......... .......... .......... .......... 66% 66.7M 1s
    ##  89350K .......... .......... .......... .......... .......... 66% 51.8M 1s
    ##  89400K .......... .......... .......... .......... .......... 66% 57.3M 1s
    ##  89450K .......... .......... .......... .......... .......... 66% 69.0M 1s
    ##  89500K .......... .......... .......... .......... .......... 66% 67.0M 1s
    ##  89550K .......... .......... .......... .......... .......... 66% 58.6M 1s
    ##  89600K .......... .......... .......... .......... .......... 66% 60.5M 1s
    ##  89650K .......... .......... .......... .......... .......... 66% 74.2M 1s
    ##  89700K .......... .......... .......... .......... .......... 66% 75.7M 1s
    ##  89750K .......... .......... .......... .......... .......... 66% 64.6M 1s
    ##  89800K .......... .......... .......... .......... .......... 67% 80.6M 1s
    ##  89850K .......... .......... .......... .......... .......... 67% 74.1M 1s
    ##  89900K .......... .......... .......... .......... .......... 67% 66.4M 1s
    ##  89950K .......... .......... .......... .......... .......... 67% 47.6M 1s
    ##  90000K .......... .......... .......... .......... .......... 67% 54.4M 1s
    ##  90050K .......... .......... .......... .......... .......... 67% 67.1M 1s
    ##  90100K .......... .......... .......... .......... .......... 67% 66.2M 1s
    ##  90150K .......... .......... .......... .......... .......... 67% 58.4M 1s
    ##  90200K .......... .......... .......... .......... .......... 67% 77.8M 1s
    ##  90250K .......... .......... .......... .......... .......... 67% 75.2M 1s
    ##  90300K .......... .......... .......... .......... .......... 67% 68.4M 1s
    ##  90350K .......... .......... .......... .......... .......... 67% 54.3M 1s
    ##  90400K .......... .......... .......... .......... .......... 67% 62.1M 1s
    ##  90450K .......... .......... .......... .......... .......... 67% 60.2M 1s
    ##  90500K .......... .......... .......... .......... .......... 67% 75.0M 1s
    ##  90550K .......... .......... .......... .......... .......... 67% 59.8M 1s
    ##  90600K .......... .......... .......... .......... .......... 67% 74.7M 1s
    ##  90650K .......... .......... .......... .......... .......... 67% 80.9M 1s
    ##  90700K .......... .......... .......... .......... .......... 67% 72.9M 1s
    ##  90750K .......... .......... .......... .......... .......... 67% 53.7M 1s
    ##  90800K .......... .......... .......... .......... .......... 67% 54.2M 1s
    ##  90850K .......... .......... .......... .......... .......... 67% 69.2M 1s
    ##  90900K .......... .......... .......... .......... .......... 67% 69.9M 1s
    ##  90950K .......... .......... .......... .......... .......... 67% 71.1M 1s
    ##  91000K .......... .......... .......... .......... .......... 67% 77.3M 1s
    ##  91050K .......... .......... .......... .......... .......... 67% 64.4M 1s
    ##  91100K .......... .......... .......... .......... .......... 67% 62.2M 1s
    ##  91150K .......... .......... .......... .......... .......... 68% 61.4M 1s
    ##  91200K .......... .......... .......... .......... .......... 68% 61.7M 1s
    ##  91250K .......... .......... .......... .......... .......... 68% 57.2M 1s
    ##  91300K .......... .......... .......... .......... .......... 68% 62.7M 1s
    ##  91350K .......... .......... .......... .......... .......... 68% 66.6M 1s
    ##  91400K .......... .......... .......... .......... .......... 68% 69.4M 1s
    ##  91450K .......... .......... .......... .......... .......... 68% 63.8M 1s
    ##  91500K .......... .......... .......... .......... .......... 68% 67.6M 1s
    ##  91550K .......... .......... .......... .......... .......... 68% 58.0M 1s
    ##  91600K .......... .......... .......... .......... .......... 68% 58.0M 1s
    ##  91650K .......... .......... .......... .......... .......... 68% 80.7M 1s
    ##  91700K .......... .......... .......... .......... .......... 68% 70.0M 1s
    ##  91750K .......... .......... .......... .......... .......... 68% 74.2M 1s
    ##  91800K .......... .......... .......... .......... .......... 68% 86.5M 1s
    ##  91850K .......... .......... .......... .......... .......... 68% 53.0M 1s
    ##  91900K .......... .......... .......... .......... .......... 68% 59.4M 1s
    ##  91950K .......... .......... .......... .......... .......... 68% 50.8M 1s
    ##  92000K .......... .......... .......... .......... .......... 68% 73.6M 1s
    ##  92050K .......... .......... .......... .......... .......... 68% 66.9M 1s
    ##  92100K .......... .......... .......... .......... .......... 68% 61.8M 1s
    ##  92150K .......... .......... .......... .......... .......... 68% 58.0M 1s
    ##  92200K .......... .......... .......... .......... .......... 68% 72.0M 1s
    ##  92250K .......... .......... .......... .......... .......... 68% 78.4M 1s
    ##  92300K .......... .......... .......... .......... .......... 68% 69.1M 1s
    ##  92350K .......... .......... .......... .......... .......... 68% 61.4M 1s
    ##  92400K .......... .......... .......... .......... .......... 68% 83.5M 1s
    ##  92450K .......... .......... .......... .......... .......... 68% 68.9M 1s
    ##  92500K .......... .......... .......... .......... .......... 69% 73.3M 1s
    ##  92550K .......... .......... .......... .......... .......... 69% 49.9M 1s
    ##  92600K .......... .......... .......... .......... .......... 69% 55.9M 1s
    ##  92650K .......... .......... .......... .......... .......... 69% 67.6M 1s
    ##  92700K .......... .......... .......... .......... .......... 69% 71.0M 1s
    ##  92750K .......... .......... .......... .......... .......... 69% 62.8M 1s
    ##  92800K .......... .......... .......... .......... .......... 69% 81.7M 1s
    ##  92850K .......... .......... .......... .......... .......... 69% 67.4M 1s
    ##  92900K .......... .......... .......... .......... .......... 69% 83.3M 1s
    ##  92950K .......... .......... .......... .......... .......... 69% 65.4M 1s
    ##  93000K .......... .......... .......... .......... .......... 69% 59.9M 1s
    ##  93050K .......... .......... .......... .......... .......... 69% 71.8M 1s
    ##  93100K .......... .......... .......... .......... .......... 69% 66.5M 1s
    ##  93150K .......... .......... .......... .......... .......... 69% 66.8M 1s
    ##  93200K .......... .......... .......... .......... .......... 69% 75.3M 1s
    ##  93250K .......... .......... .......... .......... .......... 69% 62.3M 1s
    ##  93300K .......... .......... .......... .......... .......... 69% 79.4M 1s
    ##  93350K .......... .......... .......... .......... .......... 69% 61.0M 1s
    ##  93400K .......... .......... .......... .......... .......... 69% 73.6M 1s
    ##  93450K .......... .......... .......... .......... .......... 69% 63.3M 1s
    ##  93500K .......... .......... .......... .......... .......... 69% 60.4M 1s
    ##  93550K .......... .......... .......... .......... .......... 69% 56.6M 1s
    ##  93600K .......... .......... .......... .......... .......... 69% 68.1M 1s
    ##  93650K .......... .......... .......... .......... .......... 69% 68.3M 1s
    ##  93700K .......... .......... .......... .......... .......... 69% 83.2M 1s
    ##  93750K .......... .......... .......... .......... .......... 69% 70.5M 1s
    ##  93800K .......... .......... .......... .......... .......... 70% 75.9M 1s
    ##  93850K .......... .......... .......... .......... .......... 70% 58.4M 1s
    ##  93900K .......... .......... .......... .......... .......... 70% 82.8M 1s
    ##  93950K .......... .......... .......... .......... .......... 70% 53.0M 1s
    ##  94000K .......... .......... .......... .......... .......... 70% 71.2M 1s
    ##  94050K .......... .......... .......... .......... .......... 70% 75.0M 1s
    ##  94100K .......... .......... .......... .......... .......... 70% 57.6M 1s
    ##  94150K .......... .......... .......... .......... .......... 70% 67.3M 1s
    ##  94200K .......... .......... .......... .......... .......... 70% 74.6M 1s
    ##  94250K .......... .......... .......... .......... .......... 70% 85.5M 1s
    ##  94300K .......... .......... .......... .......... .......... 70% 75.8M 1s
    ##  94350K .......... .......... .......... .......... .......... 70% 58.3M 1s
    ##  94400K .......... .......... .......... .......... .......... 70% 67.3M 1s
    ##  94450K .......... .......... .......... .......... .......... 70% 67.6M 1s
    ##  94500K .......... .......... .......... .......... .......... 70% 77.8M 1s
    ##  94550K .......... .......... .......... .......... .......... 70% 64.9M 1s
    ##  94600K .......... .......... .......... .......... .......... 70% 71.2M 1s
    ##  94650K .......... .......... .......... .......... .......... 70% 60.6M 1s
    ##  94700K .......... .......... .......... .......... .......... 70% 56.1M 1s
    ##  94750K .......... .......... .......... .......... .......... 70% 47.0M 1s
    ##  94800K .......... .......... .......... .......... .......... 70% 56.3M 1s
    ##  94850K .......... .......... .......... .......... .......... 70% 64.4M 1s
    ##  94900K .......... .......... .......... .......... .......... 70% 78.2M 1s
    ##  94950K .......... .......... .......... .......... .......... 70% 79.3M 1s
    ##  95000K .......... .......... .......... .......... .......... 70% 86.9M 1s
    ##  95050K .......... .......... .......... .......... .......... 70% 69.7M 1s
    ##  95100K .......... .......... .......... .......... .......... 70% 73.5M 1s
    ##  95150K .......... .......... .......... .......... .......... 71% 57.1M 1s
    ##  95200K .......... .......... .......... .......... .......... 71% 55.2M 1s
    ##  95250K .......... .......... .......... .......... .......... 71% 76.4M 1s
    ##  95300K .......... .......... .......... .......... .......... 71% 71.7M 1s
    ##  95350K .......... .......... .......... .......... .......... 71% 54.7M 1s
    ##  95400K .......... .......... .......... .......... .......... 71% 61.6M 1s
    ##  95450K .......... .......... .......... .......... .......... 71% 80.3M 1s
    ##  95500K .......... .......... .......... .......... .......... 71% 81.2M 1s
    ##  95550K .......... .......... .......... .......... .......... 71% 64.2M 1s
    ##  95600K .......... .......... .......... .......... .......... 71% 69.3M 1s
    ##  95650K .......... .......... .......... .......... .......... 71% 78.0M 1s
    ##  95700K .......... .......... .......... .......... .......... 71% 55.0M 1s
    ##  95750K .......... .......... .......... .......... .......... 71% 70.1M 1s
    ##  95800K .......... .......... .......... .......... .......... 71% 66.0M 1s
    ##  95850K .......... .......... .......... .......... .......... 71% 75.6M 1s
    ##  95900K .......... .......... .......... .......... .......... 71% 67.4M 1s
    ##  95950K .......... .......... .......... .......... .......... 71% 53.2M 1s
    ##  96000K .......... .......... .......... .......... .......... 71% 63.1M 1s
    ##  96050K .......... .......... .......... .......... .......... 71% 72.7M 1s
    ##  96100K .......... .......... .......... .......... .......... 71% 55.6M 1s
    ##  96150K .......... .......... .......... .......... .......... 71% 58.3M 1s
    ##  96200K .......... .......... .......... .......... .......... 71% 70.3M 1s
    ##  96250K .......... .......... .......... .......... .......... 71% 28.2M 1s
    ##  96300K .......... .......... .......... .......... .......... 71%  121M 1s
    ##  96350K .......... .......... .......... .......... .......... 71%  129M 1s
    ##  96400K .......... .......... .......... .......... .......... 71% 92.7M 1s
    ##  96450K .......... .......... .......... .......... .......... 71%  114M 1s
    ##  96500K .......... .......... .......... .......... .......... 72%  113M 1s
    ##  96550K .......... .......... .......... .......... .......... 72%  129M 1s
    ##  96600K .......... .......... .......... .......... .......... 72%  141M 1s
    ##  96650K .......... .......... .......... .......... .......... 72% 95.4M 1s
    ##  96700K .......... .......... .......... .......... .......... 72% 97.6M 1s
    ##  96750K .......... .......... .......... .......... .......... 72% 84.1M 1s
    ##  96800K .......... .......... .......... .......... .......... 72% 85.5M 1s
    ##  96850K .......... .......... .......... .......... .......... 72%  118M 1s
    ##  96900K .......... .......... .......... .......... .......... 72%  103M 1s
    ##  96950K .......... .......... .......... .......... .......... 72% 80.9M 1s
    ##  97000K .......... .......... .......... .......... .......... 72%  119M 1s
    ##  97050K .......... .......... .......... .......... .......... 72%  119M 1s
    ##  97100K .......... .......... .......... .......... .......... 72%  140M 1s
    ##  97150K .......... .......... .......... .......... .......... 72%  106M 1s
    ##  97200K .......... .......... .......... .......... .......... 72%  102M 1s
    ##  97250K .......... .......... .......... .......... .......... 72%  120M 1s
    ##  97300K .......... .......... .......... .......... .......... 72%  129M 1s
    ##  97350K .......... .......... .......... .......... .......... 72% 94.5M 1s
    ##  97400K .......... .......... .......... .......... .......... 72% 95.4M 1s
    ##  97450K .......... .......... .......... .......... .......... 72%  116M 1s
    ##  97500K .......... .......... .......... .......... .......... 72%  105M 1s
    ##  97550K .......... .......... .......... .......... .......... 72% 77.5M 1s
    ##  97600K .......... .......... .......... .......... .......... 72% 85.2M 1s
    ##  97650K .......... .......... .......... .......... .......... 72%  101M 1s
    ##  97700K .......... .......... .......... .......... .......... 72%  148M 1s
    ##  97750K .......... .......... .......... .......... .......... 72% 83.3M 1s
    ##  97800K .......... .......... .......... .......... .......... 72%  127M 1s
    ##  97850K .......... .......... .......... .......... .......... 73% 99.3M 1s
    ##  97900K .......... .......... .......... .......... .......... 73%  135M 1s
    ##  97950K .......... .......... .......... .......... .......... 73% 83.3M 1s
    ##  98000K .......... .......... .......... .......... .......... 73% 90.1M 1s
    ##  98050K .......... .......... .......... .......... .......... 73%  108M 1s
    ##  98100K .......... .......... .......... .......... .......... 73%  113M 1s
    ##  98150K .......... .......... .......... .......... .......... 73% 99.7M 1s
    ##  98200K .......... .......... .......... .......... .......... 73%  109M 1s
    ##  98250K .......... .......... .......... .......... .......... 73%  110M 1s
    ##  98300K .......... .......... .......... .......... .......... 73%  118M 1s
    ##  98350K .......... .......... .......... .......... .......... 73% 96.8M 1s
    ##  98400K .......... .......... .......... .......... .......... 73%  127M 1s
    ##  98450K .......... .......... .......... .......... .......... 73%  154M 1s
    ##  98500K .......... .......... .......... .......... .......... 73% 91.3M 1s
    ##  98550K .......... .......... .......... .......... .......... 73%  100M 1s
    ##  98600K .......... .......... .......... .......... .......... 73% 93.5M 1s
    ##  98650K .......... .......... .......... .......... .......... 73% 87.3M 1s
    ##  98700K .......... .......... .......... .......... .......... 73%  101M 1s
    ##  98750K .......... .......... .......... .......... .......... 73% 86.9M 1s
    ##  98800K .......... .......... .......... .......... .......... 73% 73.5M 1s
    ##  98850K .......... .......... .......... .......... .......... 73% 59.1M 1s
    ##  98900K .......... .......... .......... .......... .......... 73%  116M 1s
    ##  98950K .......... .......... .......... .......... .......... 73%  105M 1s
    ##  99000K .......... .......... .......... .......... .......... 73% 81.7M 1s
    ##  99050K .......... .......... .......... .......... .......... 73% 87.1M 1s
    ##  99100K .......... .......... .......... .......... .......... 73%  105M 1s
    ##  99150K .......... .......... .......... .......... .......... 73% 71.5M 1s
    ##  99200K .......... .......... .......... .......... .......... 74% 84.4M 1s
    ##  99250K .......... .......... .......... .......... .......... 74% 96.4M 1s
    ##  99300K .......... .......... .......... .......... .......... 74%  106M 1s
    ##  99350K .......... .......... .......... .......... .......... 74% 76.5M 1s
    ##  99400K .......... .......... .......... .......... .......... 74% 80.3M 1s
    ##  99450K .......... .......... .......... .......... .......... 74% 82.8M 1s
    ##  99500K .......... .......... .......... .......... .......... 74% 97.9M 1s
    ##  99550K .......... .......... .......... .......... .......... 74% 64.1M 1s
    ##  99600K .......... .......... .......... .......... .......... 74%  118M 1s
    ##  99650K .......... .......... .......... .......... .......... 74% 66.3M 1s
    ##  99700K .......... .......... .......... .......... .......... 74% 91.0M 1s
    ##  99750K .......... .......... .......... .......... .......... 74% 79.0M 1s
    ##  99800K .......... .......... .......... .......... .......... 74%  119M 1s
    ##  99850K .......... .......... .......... .......... .......... 74%  113M 1s
    ##  99900K .......... .......... .......... .......... .......... 74% 60.4M 1s
    ##  99950K .......... .......... .......... .......... .......... 74% 67.7M 1s
    ## 100000K .......... .......... .......... .......... .......... 74%  102M 1s
    ## 100050K .......... .......... .......... .......... .......... 74% 97.6M 1s
    ## 100100K .......... .......... .......... .......... .......... 74% 87.1M 1s
    ## 100150K .......... .......... .......... .......... .......... 74% 69.0M 1s
    ## 100200K .......... .......... .......... .......... .......... 74%  114M 1s
    ## 100250K .......... .......... .......... .......... .......... 74% 76.4M 1s
    ## 100300K .......... .......... .......... .......... .......... 74%  100M 1s
    ## 100350K .......... .......... .......... .......... .......... 74% 69.4M 1s
    ## 100400K .......... .......... .......... .......... .......... 74% 50.5M 1s
    ## 100450K .......... .......... .......... .......... .......... 74%  118M 1s
    ## 100500K .......... .......... .......... .......... .......... 75%  107M 1s
    ## 100550K .......... .......... .......... .......... .......... 75%  107M 1s
    ## 100600K .......... .......... .......... .......... .......... 75%  101M 1s
    ## 100650K .......... .......... .......... .......... .......... 75% 74.5M 1s
    ## 100700K .......... .......... .......... .......... .......... 75%  112M 1s
    ## 100750K .......... .......... .......... .......... .......... 75% 71.7M 1s
    ## 100800K .......... .......... .......... .......... .......... 75% 90.9M 1s
    ## 100850K .......... .......... .......... .......... .......... 75% 87.5M 1s
    ## 100900K .......... .......... .......... .......... .......... 75% 59.4M 1s
    ## 100950K .......... .......... .......... .......... .......... 75% 98.5M 1s
    ## 101000K .......... .......... .......... .......... .......... 75%  123M 1s
    ## 101050K .......... .......... .......... .......... .......... 75%  113M 1s
    ## 101100K .......... .......... .......... .......... .......... 75%  101M 1s
    ## 101150K .......... .......... .......... .......... .......... 75% 44.9M 1s
    ## 101200K .......... .......... .......... .......... .......... 75% 95.0M 1s
    ## 101250K .......... .......... .......... .......... .......... 75% 98.9M 1s
    ## 101300K .......... .......... .......... .......... .......... 75%  117M 1s
    ## 101350K .......... .......... .......... .......... .......... 75% 88.8M 1s
    ## 101400K .......... .......... .......... .......... .......... 75% 71.1M 1s
    ## 101450K .......... .......... .......... .......... .......... 75% 98.4M 1s
    ## 101500K .......... .......... .......... .......... .......... 75% 70.6M 1s
    ## 101550K .......... .......... .......... .......... .......... 75% 66.7M 1s
    ## 101600K .......... .......... .......... .......... .......... 75%  103M 1s
    ## 101650K .......... .......... .......... .......... .......... 75%  102M 1s
    ## 101700K .......... .......... .......... .......... .......... 75%  116M 1s
    ## 101750K .......... .......... .......... .......... .......... 75% 46.3M 1s
    ## 101800K .......... .......... .......... .......... .......... 75% 87.4M 1s
    ## 101850K .......... .......... .......... .......... .......... 76% 97.2M 1s
    ## 101900K .......... .......... .......... .......... .......... 76%  117M 1s
    ## 101950K .......... .......... .......... .......... .......... 76%  104M 1s
    ## 102000K .......... .......... .......... .......... .......... 76% 95.2M 1s
    ## 102050K .......... .......... .......... .......... .......... 76% 56.6M 1s
    ## 102100K .......... .......... .......... .......... .......... 76% 98.8M 1s
    ## 102150K .......... .......... .......... .......... .......... 76% 82.5M 1s
    ## 102200K .......... .......... .......... .......... .......... 76%  114M 1s
    ## 102250K .......... .......... .......... .......... .......... 76%  118M 1s
    ## 102300K .......... .......... .......... .......... .......... 76% 67.3M 1s
    ## 102350K .......... .......... .......... .......... .......... 76% 47.3M 1s
    ## 102400K .......... .......... .......... .......... .......... 76% 95.8M 1s
    ## 102450K .......... .......... .......... .......... .......... 76%  119M 1s
    ## 102500K .......... .......... .......... .......... .......... 76% 98.9M 1s
    ## 102550K .......... .......... .......... .......... .......... 76% 95.0M 1s
    ## 102600K .......... .......... .......... .......... .......... 76% 81.5M 1s
    ## 102650K .......... .......... .......... .......... .......... 76% 84.0M 1s
    ## 102700K .......... .......... .......... .......... .......... 76% 86.3M 1s
    ## 102750K .......... .......... .......... .......... .......... 76% 91.6M 1s
    ## 102800K .......... .......... .......... .......... .......... 76%  114M 1s
    ## 102850K .......... .......... .......... .......... .......... 76%  118M 1s
    ## 102900K .......... .......... .......... .......... .......... 76% 80.1M 1s
    ## 102950K .......... .......... .......... .......... .......... 76% 37.9M 1s
    ## 103000K .......... .......... .......... .......... .......... 76%  102M 1s
    ## 103050K .......... .......... .......... .......... .......... 76% 92.0M 1s
    ## 103100K .......... .......... .......... .......... .......... 76%  106M 1s
    ## 103150K .......... .......... .......... .......... .......... 76%  112M 1s
    ## 103200K .......... .......... .......... .......... .......... 77% 75.7M 1s
    ## 103250K .......... .......... .......... .......... .......... 77%  114M 1s
    ## 103300K .......... .......... .......... .......... .......... 77% 84.4M 1s
    ## 103350K .......... .......... .......... .......... .......... 77%  119M 1s
    ## 103400K .......... .......... .......... .......... .......... 77%  103M 1s
    ## 103450K .......... .......... .......... .......... .......... 77% 71.0M 1s
    ## 103500K .......... .......... .......... .......... .......... 77% 78.9M 1s
    ## 103550K .......... .......... .......... .......... .......... 77% 88.9M 1s
    ## 103600K .......... .......... .......... .......... .......... 77%  108M 1s
    ## 103650K .......... .......... .......... .......... .......... 77% 80.1M 1s
    ## 103700K .......... .......... .......... .......... .......... 77% 84.9M 1s
    ## 103750K .......... .......... .......... .......... .......... 77% 73.3M 1s
    ## 103800K .......... .......... .......... .......... .......... 77% 63.1M 1s
    ## 103850K .......... .......... .......... .......... .......... 77%  119M 1s
    ## 103900K .......... .......... .......... .......... .......... 77% 81.1M 1s
    ## 103950K .......... .......... .......... .......... .......... 77% 80.6M 1s
    ## 104000K .......... .......... .......... .......... .......... 77% 97.4M 1s
    ## 104050K .......... .......... .......... .......... .......... 77% 73.7M 1s
    ## 104100K .......... .......... .......... .......... .......... 77% 88.5M 1s
    ## 104150K .......... .......... .......... .......... .......... 77% 85.8M 1s
    ## 104200K .......... .......... .......... .......... .......... 77% 97.7M 1s
    ## 104250K .......... .......... .......... .......... .......... 77%  105M 1s
    ## 104300K .......... .......... .......... .......... .......... 77% 89.0M 1s
    ## 104350K .......... .......... .......... .......... .......... 77% 53.0M 1s
    ## 104400K .......... .......... .......... .......... .......... 77% 86.9M 1s
    ## 104450K .......... .......... .......... .......... .......... 77% 79.5M 1s
    ## 104500K .......... .......... .......... .......... .......... 77%  101M 1s
    ## 104550K .......... .......... .......... .......... .......... 78%  101M 1s
    ## 104600K .......... .......... .......... .......... .......... 78% 50.6M 1s
    ## 104650K .......... .......... .......... .......... .......... 78% 89.7M 1s
    ## 104700K .......... .......... .......... .......... .......... 78% 95.8M 1s
    ## 104750K .......... .......... .......... .......... .......... 78% 97.7M 1s
    ## 104800K .......... .......... .......... .......... .......... 78%  115M 1s
    ## 104850K .......... .......... .......... .......... .......... 78%  107M 1s
    ## 104900K .......... .......... .......... .......... .......... 78% 72.3M 1s
    ## 104950K .......... .......... .......... .......... .......... 78% 90.6M 1s
    ## 105000K .......... .......... .......... .......... .......... 78% 70.9M 1s
    ## 105050K .......... .......... .......... .......... .......... 78% 93.9M 1s
    ## 105100K .......... .......... .......... .......... .......... 78% 90.7M 1s
    ## 105150K .......... .......... .......... .......... .......... 78% 85.7M 1s
    ## 105200K .......... .......... .......... .......... .......... 78% 59.9M 1s
    ## 105250K .......... .......... .......... .......... .......... 78% 80.7M 1s
    ## 105300K .......... .......... .......... .......... .......... 78% 89.8M 1s
    ## 105350K .......... .......... .......... .......... .......... 78% 94.8M 1s
    ## 105400K .......... .......... .......... .......... .......... 78%  123M 1s
    ## 105450K .......... .......... .......... .......... .......... 78%  101M 1s
    ## 105500K .......... .......... .......... .......... .......... 78% 72.5M 1s
    ## 105550K .......... .......... .......... .......... .......... 78% 78.0M 1s
    ## 105600K .......... .......... .......... .......... .......... 78% 76.8M 1s
    ## 105650K .......... .......... .......... .......... .......... 78% 92.0M 1s
    ## 105700K .......... .......... .......... .......... .......... 78%  121M 1s
    ## 105750K .......... .......... .......... .......... .......... 78% 57.6M 1s
    ## 105800K .......... .......... .......... .......... .......... 78%  102M 1s
    ## 105850K .......... .......... .......... .......... .......... 78% 78.9M 1s
    ## 105900K .......... .......... .......... .......... .......... 79% 93.2M 1s
    ## 105950K .......... .......... .......... .......... .......... 79% 90.5M 1s
    ## 106000K .......... .......... .......... .......... .......... 79%  106M 1s
    ## 106050K .......... .......... .......... .......... .......... 79%  117M 1s
    ## 106100K .......... .......... .......... .......... .......... 79% 51.1M 1s
    ## 106150K .......... .......... .......... .......... .......... 79% 84.6M 1s
    ## 106200K .......... .......... .......... .......... .......... 79%  107M 1s
    ## 106250K .......... .......... .......... .......... .......... 79%  114M 1s
    ## 106300K .......... .......... .......... .......... .......... 79% 77.7M 1s
    ## 106350K .......... .......... .......... .......... .......... 79% 49.8M 1s
    ## 106400K .......... .......... .......... .......... .......... 79%  119M 1s
    ## 106450K .......... .......... .......... .......... .......... 79%  110M 1s
    ## 106500K .......... .......... .......... .......... .......... 79% 97.1M 1s
    ## 106550K .......... .......... .......... .......... .......... 79% 89.9M 1s
    ## 106600K .......... .......... .......... .......... .......... 79%  116M 1s
    ## 106650K .......... .......... .......... .......... .......... 79% 91.4M 1s
    ## 106700K .......... .......... .......... .......... .......... 79% 38.8M 1s
    ## 106750K .......... .......... .......... .......... .......... 79% 74.7M 1s
    ## 106800K .......... .......... .......... .......... .......... 79%  117M 1s
    ## 106850K .......... .......... .......... .......... .......... 79%  102M 1s
    ## 106900K .......... .......... .......... .......... .......... 79%  125M 1s
    ## 106950K .......... .......... .......... .......... .......... 79% 83.1M 1s
    ## 107000K .......... .......... .......... .......... .......... 79% 81.8M 1s
    ## 107050K .......... .......... .......... .......... .......... 79%  105M 1s
    ## 107100K .......... .......... .......... .......... .......... 79%  141M 1s
    ## 107150K .......... .......... .......... .......... .......... 79% 86.1M 1s
    ## 107200K .......... .......... .......... .......... .......... 79% 88.2M 1s
    ## 107250K .......... .......... .......... .......... .......... 80% 55.9M 1s
    ## 107300K .......... .......... .......... .......... .......... 80% 99.1M 1s
    ## 107350K .......... .......... .......... .......... .......... 80%  105M 1s
    ## 107400K .......... .......... .......... .......... .......... 80% 93.4M 1s
    ## 107450K .......... .......... .......... .......... .......... 80% 76.7M 1s
    ## 107500K .......... .......... .......... .......... .......... 80%  105M 1s
    ## 107550K .......... .......... .......... .......... .......... 80% 63.0M 0s
    ## 107600K .......... .......... .......... .......... .......... 80% 97.3M 0s
    ## 107650K .......... .......... .......... .......... .......... 80% 91.6M 0s
    ## 107700K .......... .......... .......... .......... .......... 80% 85.2M 0s
    ## 107750K .......... .......... .......... .......... .......... 80% 82.7M 0s
    ## 107800K .......... .......... .......... .......... .......... 80%  110M 0s
    ## 107850K .......... .......... .......... .......... .......... 80% 78.8M 0s
    ## 107900K .......... .......... .......... .......... .......... 80% 98.9M 0s
    ## 107950K .......... .......... .......... .......... .......... 80% 69.8M 0s
    ## 108000K .......... .......... .......... .......... .......... 80%  110M 0s
    ## 108050K .......... .......... .......... .......... .......... 80% 99.8M 0s
    ## 108100K .......... .......... .......... .......... .......... 80% 75.3M 0s
    ## 108150K .......... .......... .......... .......... .......... 80% 63.8M 0s
    ## 108200K .......... .......... .......... .......... .......... 80%  100M 0s
    ## 108250K .......... .......... .......... .......... .......... 80%  101M 0s
    ## 108300K .......... .......... .......... .......... .......... 80%  103M 0s
    ## 108350K .......... .......... .......... .......... .......... 80% 60.2M 0s
    ## 108400K .......... .......... .......... .......... .......... 80%  108M 0s
    ## 108450K .......... .......... .......... .......... .......... 80% 58.9M 0s
    ## 108500K .......... .......... .......... .......... .......... 80% 84.9M 0s
    ## 108550K .......... .......... .......... .......... .......... 81% 91.0M 0s
    ## 108600K .......... .......... .......... .......... .......... 81% 90.0M 0s
    ## 108650K .......... .......... .......... .......... .......... 81%  100M 0s
    ## 108700K .......... .......... .......... .......... .......... 81% 53.8M 0s
    ## 108750K .......... .......... .......... .......... .......... 81% 66.6M 0s
    ## 108800K .......... .......... .......... .......... .......... 81% 86.3M 0s
    ## 108850K .......... .......... .......... .......... .......... 81%  110M 0s
    ## 108900K .......... .......... .......... .......... .......... 81% 95.1M 0s
    ## 108950K .......... .......... .......... .......... .......... 81% 87.4M 0s
    ## 109000K .......... .......... .......... .......... .......... 81%  110M 0s
    ## 109050K .......... .......... .......... .......... .......... 81% 85.4M 0s
    ## 109100K .......... .......... .......... .......... .......... 81% 94.3M 0s
    ## 109150K .......... .......... .......... .......... .......... 81% 79.3M 0s
    ## 109200K .......... .......... .......... .......... .......... 81%  109M 0s
    ## 109250K .......... .......... .......... .......... .......... 81%  106M 0s
    ## 109300K .......... .......... .......... .......... .......... 81% 69.7M 0s
    ## 109350K .......... .......... .......... .......... .......... 81% 68.9M 0s
    ## 109400K .......... .......... .......... .......... .......... 81%  116M 0s
    ## 109450K .......... .......... .......... .......... .......... 81%  100M 0s
    ## 109500K .......... .......... .......... .......... .......... 81%  106M 0s
    ## 109550K .......... .......... .......... .......... .......... 81%  104M 0s
    ## 109600K .......... .......... .......... .......... .......... 81% 46.1M 0s
    ## 109650K .......... .......... .......... .......... .......... 81% 92.4M 0s
    ## 109700K .......... .......... .......... .......... .......... 81% 92.0M 0s
    ## 109750K .......... .......... .......... .......... .......... 81% 79.5M 0s
    ## 109800K .......... .......... .......... .......... .......... 81% 81.2M 0s
    ## 109850K .......... .......... .......... .......... .......... 81%  120M 0s
    ## 109900K .......... .......... .......... .......... .......... 82% 70.9M 0s
    ## 109950K .......... .......... .......... .......... .......... 82% 86.4M 0s
    ## 110000K .......... .......... .......... .......... .......... 82%  101M 0s
    ## 110050K .......... .......... .......... .......... .......... 82% 95.0M 0s
    ## 110100K .......... .......... .......... .......... .......... 82% 91.8M 0s
    ## 110150K .......... .......... .......... .......... .......... 82%  100M 0s
    ## 110200K .......... .......... .......... .......... .......... 82% 55.4M 0s
    ## 110250K .......... .......... .......... .......... .......... 82%  108M 0s
    ## 110300K .......... .......... .......... .......... .......... 82% 93.6M 0s
    ## 110350K .......... .......... .......... .......... .......... 82% 83.8M 0s
    ## 110400K .......... .......... .......... .......... .......... 82% 87.4M 0s
    ## 110450K .......... .......... .......... .......... .......... 82% 70.8M 0s
    ## 110500K .......... .......... .......... .......... .......... 82% 89.1M 0s
    ## 110550K .......... .......... .......... .......... .......... 82% 85.1M 0s
    ## 110600K .......... .......... .......... .......... .......... 82% 92.5M 0s
    ## 110650K .......... .......... .......... .......... .......... 82% 88.4M 0s
    ## 110700K .......... .......... .......... .......... .......... 82%  138M 0s
    ## 110750K .......... .......... .......... .......... .......... 82% 45.6M 0s
    ## 110800K .......... .......... .......... .......... .......... 82%  116M 0s
    ## 110850K .......... .......... .......... .......... .......... 82%  116M 0s
    ## 110900K .......... .......... .......... .......... .......... 82%  105M 0s
    ## 110950K .......... .......... .......... .......... .......... 82%  109M 0s
    ## 111000K .......... .......... .......... .......... .......... 82% 53.5M 0s
    ## 111050K .......... .......... .......... .......... .......... 82%  100M 0s
    ## 111100K .......... .......... .......... .......... .......... 82%  108M 0s
    ## 111150K .......... .......... .......... .......... .......... 82% 96.6M 0s
    ## 111200K .......... .......... .......... .......... .......... 82%  121M 0s
    ## 111250K .......... .......... .......... .......... .......... 83%  110M 0s
    ## 111300K .......... .......... .......... .......... .......... 83% 43.0M 0s
    ## 111350K .......... .......... .......... .......... .......... 83% 89.1M 0s
    ## 111400K .......... .......... .......... .......... .......... 83%  111M 0s
    ## 111450K .......... .......... .......... .......... .......... 83% 89.6M 0s
    ## 111500K .......... .......... .......... .......... .......... 83%  110M 0s
    ## 111550K .......... .......... .......... .......... .......... 83% 73.9M 0s
    ## 111600K .......... .......... .......... .......... .......... 83% 60.5M 0s
    ## 111650K .......... .......... .......... .......... .......... 83% 85.6M 0s
    ## 111700K .......... .......... .......... .......... .......... 83% 88.8M 0s
    ## 111750K .......... .......... .......... .......... .......... 83% 84.9M 0s
    ## 111800K .......... .......... .......... .......... .......... 83% 93.1M 0s
    ## 111850K .......... .......... .......... .......... .......... 83%  148M 0s
    ## 111900K .......... .......... .......... .......... .......... 83% 51.5M 0s
    ## 111950K .......... .......... .......... .......... .......... 83% 80.1M 0s
    ## 112000K .......... .......... .......... .......... .......... 83% 87.9M 0s
    ## 112050K .......... .......... .......... .......... .......... 83%  105M 0s
    ## 112100K .......... .......... .......... .......... .......... 83%  138M 0s
    ## 112150K .......... .......... .......... .......... .......... 83% 66.6M 0s
    ## 112200K .......... .......... .......... .......... .......... 83% 96.1M 0s
    ## 112250K .......... .......... .......... .......... .......... 83%  103M 0s
    ## 112300K .......... .......... .......... .......... .......... 83%  116M 0s
    ## 112350K .......... .......... .......... .......... .......... 83% 85.4M 0s
    ## 112400K .......... .......... .......... .......... .......... 83% 92.8M 0s
    ## 112450K .......... .......... .......... .......... .......... 83% 56.3M 0s
    ## 112500K .......... .......... .......... .......... .......... 83% 72.4M 0s
    ## 112550K .......... .......... .......... .......... .......... 83% 82.7M 0s
    ## 112600K .......... .......... .......... .......... .......... 84% 90.4M 0s
    ## 112650K .......... .......... .......... .......... .......... 84%  115M 0s
    ## 112700K .......... .......... .......... .......... .......... 84% 98.6M 0s
    ## 112750K .......... .......... .......... .......... .......... 84% 63.1M 0s
    ## 112800K .......... .......... .......... .......... .......... 84% 73.1M 0s
    ## 112850K .......... .......... .......... .......... .......... 84%  107M 0s
    ## 112900K .......... .......... .......... .......... .......... 84% 87.4M 0s
    ## 112950K .......... .......... .......... .......... .......... 84% 90.9M 0s
    ## 113000K .......... .......... .......... .......... .......... 84%  139M 0s
    ## 113050K .......... .......... .......... .......... .......... 84% 57.0M 0s
    ## 113100K .......... .......... .......... .......... .......... 84%  100M 0s
    ## 113150K .......... .......... .......... .......... .......... 84% 68.0M 0s
    ## 113200K .......... .......... .......... .......... .......... 84%  118M 0s
    ## 113250K .......... .......... .......... .......... .......... 84%  124M 0s
    ## 113300K .......... .......... .......... .......... .......... 84% 80.0M 0s
    ## 113350K .......... .......... .......... .......... .......... 84% 51.7M 0s
    ## 113400K .......... .......... .......... .......... .......... 84%  122M 0s
    ## 113450K .......... .......... .......... .......... .......... 84%  101M 0s
    ## 113500K .......... .......... .......... .......... .......... 84%  106M 0s
    ## 113550K .......... .......... .......... .......... .......... 84%  123M 0s
    ## 113600K .......... .......... .......... .......... .......... 84% 65.9M 0s
    ## 113650K .......... .......... .......... .......... .......... 84% 73.2M 0s
    ## 113700K .......... .......... .......... .......... .......... 84% 80.7M 0s
    ## 113750K .......... .......... .......... .......... .......... 84% 89.9M 0s
    ## 113800K .......... .......... .......... .......... .......... 84%  120M 0s
    ## 113850K .......... .......... .......... .......... .......... 84%  141M 0s
    ## 113900K .......... .......... .......... .......... .......... 84% 54.8M 0s
    ## 113950K .......... .......... .......... .......... .......... 85% 54.1M 0s
    ## 114000K .......... .......... .......... .......... .......... 85% 81.0M 0s
    ## 114050K .......... .......... .......... .......... .......... 85%  105M 0s
    ## 114100K .......... .......... .......... .......... .......... 85%  153M 0s
    ## 114150K .......... .......... .......... .......... .......... 85%  113M 0s
    ## 114200K .......... .......... .......... .......... .......... 85% 61.9M 0s
    ## 114250K .......... .......... .......... .......... .......... 85% 86.5M 0s
    ## 114300K .......... .......... .......... .......... .......... 85% 88.5M 0s
    ## 114350K .......... .......... .......... .......... .......... 85%  119M 0s
    ## 114400K .......... .......... .......... .......... .......... 85%  123M 0s
    ## 114450K .......... .......... .......... .......... .......... 85% 51.4M 0s
    ## 114500K .......... .......... .......... .......... .......... 85% 52.5M 0s
    ## 114550K .......... .......... .......... .......... .......... 85% 73.5M 0s
    ## 114600K .......... .......... .......... .......... .......... 85% 88.0M 0s
    ## 114650K .......... .......... .......... .......... .......... 85%  112M 0s
    ## 114700K .......... .......... .......... .......... .......... 85%  136M 0s
    ## 114750K .......... .......... .......... .......... .......... 85% 71.9M 0s
    ## 114800K .......... .......... .......... .......... .......... 85% 95.5M 0s
    ## 114850K .......... .......... .......... .......... .......... 85% 94.7M 0s
    ## 114900K .......... .......... .......... .......... .......... 85%  101M 0s
    ## 114950K .......... .......... .......... .......... .......... 85%  122M 0s
    ## 115000K .......... .......... .......... .......... .......... 85% 62.0M 0s
    ## 115050K .......... .......... .......... .......... .......... 85% 75.2M 0s
    ## 115100K .......... .......... .......... .......... .......... 85% 83.1M 0s
    ## 115150K .......... .......... .......... .......... .......... 85% 69.6M 0s
    ## 115200K .......... .......... .......... .......... .......... 85% 98.7M 0s
    ## 115250K .......... .......... .......... .......... .......... 86%  143M 0s
    ## 115300K .......... .......... .......... .......... .......... 86% 54.6M 0s
    ## 115350K .......... .......... .......... .......... .......... 86% 70.8M 0s
    ## 115400K .......... .......... .......... .......... .......... 86% 84.4M 0s
    ## 115450K .......... .......... .......... .......... .......... 86% 95.3M 0s
    ## 115500K .......... .......... .......... .......... .......... 86%  114M 0s
    ## 115550K .......... .......... .......... .......... .......... 86% 95.5M 0s
    ## 115600K .......... .......... .......... .......... .......... 86%  101M 0s
    ## 115650K .......... .......... .......... .......... .......... 86% 66.2M 0s
    ## 115700K .......... .......... .......... .......... .......... 86% 92.0M 0s
    ## 115750K .......... .......... .......... .......... .......... 86% 89.6M 0s
    ## 115800K .......... .......... .......... .......... .......... 86%  130M 0s
    ## 115850K .......... .......... .......... .......... .......... 86%  143M 0s
    ## 115900K .......... .......... .......... .......... .......... 86% 48.9M 0s
    ## 115950K .......... .......... .......... .......... .......... 86% 79.7M 0s
    ## 116000K .......... .......... .......... .......... .......... 86%  106M 0s
    ## 116050K .......... .......... .......... .......... .......... 86%  129M 0s
    ## 116100K .......... .......... .......... .......... .......... 86%  146M 0s
    ## 116150K .......... .......... .......... .......... .......... 86% 47.8M 0s
    ## 116200K .......... .......... .......... .......... .......... 86%  108M 0s
    ## 116250K .......... .......... .......... .......... .......... 86% 64.6M 0s
    ## 116300K .......... .......... .......... .......... .......... 86% 80.0M 0s
    ## 116350K .......... .......... .......... .......... .......... 86% 79.5M 0s
    ## 116400K .......... .......... .......... .......... .......... 86%  107M 0s
    ## 116450K .......... .......... .......... .......... .......... 86% 98.2M 0s
    ## 116500K .......... .......... .......... .......... .......... 86% 79.5M 0s
    ## 116550K .......... .......... .......... .......... .......... 86% 93.2M 0s
    ## 116600K .......... .......... .......... .......... .......... 87%  102M 0s
    ## 116650K .......... .......... .......... .......... .......... 87%  145M 0s
    ## 116700K .......... .......... .......... .......... .......... 87%  126M 0s
    ## 116750K .......... .......... .......... .......... .......... 87% 55.1M 0s
    ## 116800K .......... .......... .......... .......... .......... 87% 56.4M 0s
    ## 116850K .......... .......... .......... .......... .......... 87%  100M 0s
    ## 116900K .......... .......... .......... .......... .......... 87%  103M 0s
    ## 116950K .......... .......... .......... .......... .......... 87%  118M 0s
    ## 117000K .......... .......... .......... .......... .......... 87%  148M 0s
    ## 117050K .......... .......... .......... .......... .......... 87% 40.6M 0s
    ## 117100K .......... .......... .......... .......... .......... 87% 98.5M 0s
    ## 117150K .......... .......... .......... .......... .......... 87% 79.5M 0s
    ## 117200K .......... .......... .......... .......... .......... 87%  136M 0s
    ## 117250K .......... .......... .......... .......... .......... 87%  130M 0s
    ## 117300K .......... .......... .......... .......... .......... 87% 71.9M 0s
    ## 117350K .......... .......... .......... .......... .......... 87% 41.2M 0s
    ## 117400K .......... .......... .......... .......... .......... 87% 92.8M 0s
    ## 117450K .......... .......... .......... .......... .......... 87%  102M 0s
    ## 117500K .......... .......... .......... .......... .......... 87%  117M 0s
    ## 117550K .......... .......... .......... .......... .......... 87%  107M 0s
    ## 117600K .......... .......... .......... .......... .......... 87% 95.6M 0s
    ## 117650K .......... .......... .......... .......... .......... 87% 89.4M 0s
    ## 117700K .......... .......... .......... .......... .......... 87%  104M 0s
    ## 117750K .......... .......... .......... .......... .......... 87% 90.2M 0s
    ## 117800K .......... .......... .......... .......... .......... 87%  132M 0s
    ## 117850K .......... .......... .......... .......... .......... 87% 87.6M 0s
    ## 117900K .......... .......... .......... .......... .......... 87% 67.8M 0s
    ## 117950K .......... .......... .......... .......... .......... 88% 46.7M 0s
    ## 118000K .......... .......... .......... .......... .......... 88% 83.6M 0s
    ## 118050K .......... .......... .......... .......... .......... 88% 95.0M 0s
    ## 118100K .......... .......... .......... .......... .......... 88%  123M 0s
    ## 118150K .......... .......... .......... .......... .......... 88%  122M 0s
    ## 118200K .......... .......... .......... .......... .......... 88% 60.3M 0s
    ## 118250K .......... .......... .......... .......... .......... 88% 95.3M 0s
    ## 118300K .......... .......... .......... .......... .......... 88% 90.5M 0s
    ## 118350K .......... .......... .......... .......... .......... 88% 81.8M 0s
    ## 118400K .......... .......... .......... .......... .......... 88%  115M 0s
    ## 118450K .......... .......... .......... .......... .......... 88%  141M 0s
    ## 118500K .......... .......... .......... .......... .......... 88% 56.1M 0s
    ## 118550K .......... .......... .......... .......... .......... 88% 86.7M 0s
    ## 118600K .......... .......... .......... .......... .......... 88% 82.1M 0s
    ## 118650K .......... .......... .......... .......... .......... 88% 88.5M 0s
    ## 118700K .......... .......... .......... .......... .......... 88%  119M 0s
    ## 118750K .......... .......... .......... .......... .......... 88% 44.9M 0s
    ## 118800K .......... .......... .......... .......... .......... 88%  115M 0s
    ## 118850K .......... .......... .......... .......... .......... 88% 91.5M 0s
    ## 118900K .......... .......... .......... .......... .......... 88% 79.5M 0s
    ## 118950K .......... .......... .......... .......... .......... 88% 87.6M 0s
    ## 119000K .......... .......... .......... .......... .......... 88%  115M 0s
    ## 119050K .......... .......... .......... .......... .......... 88%  102M 0s
    ## 119100K .......... .......... .......... .......... .......... 88% 81.8M 0s
    ## 119150K .......... .......... .......... .......... .......... 88% 74.8M 0s
    ## 119200K .......... .......... .......... .......... .......... 88% 84.4M 0s
    ## 119250K .......... .......... .......... .......... .......... 88%  129M 0s
    ## 119300K .......... .......... .......... .......... .......... 89% 53.4M 0s
    ## 119350K .......... .......... .......... .......... .......... 89% 83.4M 0s
    ## 119400K .......... .......... .......... .......... .......... 89% 99.7M 0s
    ## 119450K .......... .......... .......... .......... .......... 89%  126M 0s
    ## 119500K .......... .......... .......... .......... .......... 89% 95.6M 0s
    ## 119550K .......... .......... .......... .......... .......... 89% 90.5M 0s
    ## 119600K .......... .......... .......... .......... .......... 89%  137M 0s
    ## 119650K .......... .......... .......... .......... .......... 89% 51.9M 0s
    ## 119700K .......... .......... .......... .......... .......... 89%  112M 0s
    ## 119750K .......... .......... .......... .......... .......... 89% 86.1M 0s
    ## 119800K .......... .......... .......... .......... .......... 89%  140M 0s
    ## 119850K .......... .......... .......... .......... .......... 89%  131M 0s
    ## 119900K .......... .......... .......... .......... .......... 89% 42.2M 0s
    ## 119950K .......... .......... .......... .......... .......... 89% 60.8M 0s
    ## 120000K .......... .......... .......... .......... .......... 89% 98.8M 0s
    ## 120050K .......... .......... .......... .......... .......... 89%  121M 0s
    ## 120100K .......... .......... .......... .......... .......... 89%  104M 0s
    ## 120150K .......... .......... .......... .......... .......... 89%  134M 0s
    ## 120200K .......... .......... .......... .......... .......... 89% 47.3M 0s
    ## 120250K .......... .......... .......... .......... .......... 89% 74.4M 0s
    ## 120300K .......... .......... .......... .......... .......... 89% 91.3M 0s
    ## 120350K .......... .......... .......... .......... .......... 89% 96.1M 0s
    ## 120400K .......... .......... .......... .......... .......... 89%  155M 0s
    ## 120450K .......... .......... .......... .......... .......... 89%  129M 0s
    ## 120500K .......... .......... .......... .......... .......... 89% 64.0M 0s
    ## 120550K .......... .......... .......... .......... .......... 89% 57.3M 0s
    ## 120600K .......... .......... .......... .......... .......... 89% 82.1M 0s
    ## 120650K .......... .......... .......... .......... .......... 90% 87.6M 0s
    ## 120700K .......... .......... .......... .......... .......... 90%  132M 0s
    ## 120750K .......... .......... .......... .......... .......... 90%  123M 0s
    ## 120800K .......... .......... .......... .......... .......... 90% 55.8M 0s
    ## 120850K .......... .......... .......... .......... .......... 90% 79.4M 0s
    ## 120900K .......... .......... .......... .......... .......... 90% 82.4M 0s
    ## 120950K .......... .......... .......... .......... .......... 90%  114M 0s
    ## 121000K .......... .......... .......... .......... .......... 90%  127M 0s
    ## 121050K .......... .......... .......... .......... .......... 90%  106M 0s
    ## 121100K .......... .......... .......... .......... .......... 90% 48.7M 0s
    ## 121150K .......... .......... .......... .......... .......... 90% 73.6M 0s
    ## 121200K .......... .......... .......... .......... .......... 90% 99.0M 0s
    ## 121250K .......... .......... .......... .......... .......... 90%  118M 0s
    ## 121300K .......... .......... .......... .......... .......... 90%  149M 0s
    ## 121350K .......... .......... .......... .......... .......... 90% 39.4M 0s
    ## 121400K .......... .......... .......... .......... .......... 90% 75.9M 0s
    ## 121450K .......... .......... .......... .......... .......... 90% 91.4M 0s
    ## 121500K .......... .......... .......... .......... .......... 90%  129M 0s
    ## 121550K .......... .......... .......... .......... .......... 90%  125M 0s
    ## 121600K .......... .......... .......... .......... .......... 90%  122M 0s
    ## 121650K .......... .......... .......... .......... .......... 90% 53.7M 0s
    ## 121700K .......... .......... .......... .......... .......... 90% 98.5M 0s
    ## 121750K .......... .......... .......... .......... .......... 90% 92.0M 0s
    ## 121800K .......... .......... .......... .......... .......... 90%  130M 0s
    ## 121850K .......... .......... .......... .......... .......... 90%  127M 0s
    ## 121900K .......... .......... .......... .......... .......... 90%  132M 0s
    ## 121950K .......... .......... .......... .......... .......... 91% 38.8M 0s
    ## 122000K .......... .......... .......... .......... .......... 91%  120M 0s
    ## 122050K .......... .......... .......... .......... .......... 91% 99.3M 0s
    ## 122100K .......... .......... .......... .......... .......... 91%  125M 0s
    ## 122150K .......... .......... .......... .......... .......... 91% 78.1M 0s
    ## 122200K .......... .......... .......... .......... .......... 91%  102M 0s
    ## 122250K .......... .......... .......... .......... .......... 91% 48.7M 0s
    ## 122300K .......... .......... .......... .......... .......... 91% 85.8M 0s
    ## 122350K .......... .......... .......... .......... .......... 91% 78.0M 0s
    ## 122400K .......... .......... .......... .......... .......... 91%  123M 0s
    ## 122450K .......... .......... .......... .......... .......... 91%  134M 0s
    ## 122500K .......... .......... .......... .......... .......... 91% 65.8M 0s
    ## 122550K .......... .......... .......... .......... .......... 91% 88.1M 0s
    ## 122600K .......... .......... .......... .......... .......... 91%  109M 0s
    ## 122650K .......... .......... .......... .......... .......... 91%  107M 0s
    ## 122700K .......... .......... .......... .......... .......... 91%  137M 0s
    ## 122750K .......... .......... .......... .......... .......... 91% 94.1M 0s
    ## 122800K .......... .......... .......... .......... .......... 91% 44.5M 0s
    ## 122850K .......... .......... .......... .......... .......... 91% 72.9M 0s
    ## 122900K .......... .......... .......... .......... .......... 91% 93.6M 0s
    ## 122950K .......... .......... .......... .......... .......... 91%  102M 0s
    ## 123000K .......... .......... .......... .......... .......... 91%  139M 0s
    ## 123050K .......... .......... .......... .......... .......... 91%  106M 0s
    ## 123100K .......... .......... .......... .......... .......... 91% 90.2M 0s
    ## 123150K .......... .......... .......... .......... .......... 91% 70.1M 0s
    ## 123200K .......... .......... .......... .......... .......... 91%  110M 0s
    ## 123250K .......... .......... .......... .......... .......... 91%  105M 0s
    ## 123300K .......... .......... .......... .......... .......... 92%  103M 0s
    ## 123350K .......... .......... .......... .......... .......... 92%  127M 0s
    ## 123400K .......... .......... .......... .......... .......... 92% 32.7M 0s
    ## 123450K .......... .......... .......... .......... .......... 92% 73.0M 0s
    ## 123500K .......... .......... .......... .......... .......... 92% 94.0M 0s
    ## 123550K .......... .......... .......... .......... .......... 92% 99.7M 0s
    ## 123600K .......... .......... .......... .......... .......... 92%  120M 0s
    ## 123650K .......... .......... .......... .......... .......... 92%  132M 0s
    ## 123700K .......... .......... .......... .......... .......... 92% 73.7M 0s
    ## 123750K .......... .......... .......... .......... .......... 92% 83.5M 0s
    ## 123800K .......... .......... .......... .......... .......... 92% 80.1M 0s
    ## 123850K .......... .......... .......... .......... .......... 92% 99.8M 0s
    ## 123900K .......... .......... .......... .......... .......... 92%  111M 0s
    ## 123950K .......... .......... .......... .......... .......... 92% 41.5M 0s
    ## 124000K .......... .......... .......... .......... .......... 92% 78.7M 0s
    ## 124050K .......... .......... .......... .......... .......... 92% 91.4M 0s
    ## 124100K .......... .......... .......... .......... .......... 92% 92.7M 0s
    ## 124150K .......... .......... .......... .......... .......... 92% 76.5M 0s
    ## 124200K .......... .......... .......... .......... .......... 92%  113M 0s
    ## 124250K .......... .......... .......... .......... .......... 92%  122M 0s
    ## 124300K .......... .......... .......... .......... .......... 92% 79.9M 0s
    ## 124350K .......... .......... .......... .......... .......... 92% 89.7M 0s
    ## 124400K .......... .......... .......... .......... .......... 92%  125M 0s
    ## 124450K .......... .......... .......... .......... .......... 92%  148M 0s
    ## 124500K .......... .......... .......... .......... .......... 92%  148M 0s
    ## 124550K .......... .......... .......... .......... .......... 92% 43.4M 0s
    ## 124600K .......... .......... .......... .......... .......... 92%  114M 0s
    ## 124650K .......... .......... .......... .......... .......... 93%  107M 0s
    ## 124700K .......... .......... .......... .......... .......... 93%  105M 0s
    ## 124750K .......... .......... .......... .......... .......... 93% 93.9M 0s
    ## 124800K .......... .......... .......... .......... .......... 93%  124M 0s
    ## 124850K .......... .......... .......... .......... .......... 93% 36.7M 0s
    ## 124900K .......... .......... .......... .......... .......... 93%  102M 0s
    ## 124950K .......... .......... .......... .......... .......... 93%  103M 0s
    ## 125000K .......... .......... .......... .......... .......... 93%  139M 0s
    ## 125050K .......... .......... .......... .......... .......... 93%  149M 0s
    ## 125100K .......... .......... .......... .......... .......... 93% 55.2M 0s
    ## 125150K .......... .......... .......... .......... .......... 93% 70.4M 0s
    ## 125200K .......... .......... .......... .......... .......... 93%  120M 0s
    ## 125250K .......... .......... .......... .......... .......... 93% 98.4M 0s
    ## 125300K .......... .......... .......... .......... .......... 93%  106M 0s
    ## 125350K .......... .......... .......... .......... .......... 93%  121M 0s
    ## 125400K .......... .......... .......... .......... .......... 93% 37.9M 0s
    ## 125450K .......... .......... .......... .......... .......... 93% 77.1M 0s
    ## 125500K .......... .......... .......... .......... .......... 93%  124M 0s
    ## 125550K .......... .......... .......... .......... .......... 93%  113M 0s
    ## 125600K .......... .......... .......... .......... .......... 93%  132M 0s
    ## 125650K .......... .......... .......... .......... .......... 93% 91.0M 0s
    ## 125700K .......... .......... .......... .......... .......... 93% 67.1M 0s
    ## 125750K .......... .......... .......... .......... .......... 93%  102M 0s
    ## 125800K .......... .......... .......... .......... .......... 93%  128M 0s
    ## 125850K .......... .......... .......... .......... .......... 93%  102M 0s
    ## 125900K .......... .......... .......... .......... .......... 93% 92.2M 0s
    ## 125950K .......... .......... .......... .......... .......... 93% 97.2M 0s
    ## 126000K .......... .......... .......... .......... .......... 94% 39.4M 0s
    ## 126050K .......... .......... .......... .......... .......... 94% 87.1M 0s
    ## 126100K .......... .......... .......... .......... .......... 94%  101M 0s
    ## 126150K .......... .......... .......... .......... .......... 94%  121M 0s
    ## 126200K .......... .......... .......... .......... .......... 94%  145M 0s
    ## 126250K .......... .......... .......... .......... .......... 94% 53.2M 0s
    ## 126300K .......... .......... .......... .......... .......... 94%  112M 0s
    ## 126350K .......... .......... .......... .......... .......... 94% 90.9M 0s
    ## 126400K .......... .......... .......... .......... .......... 94%  101M 0s
    ## 126450K .......... .......... .......... .......... .......... 94%  118M 0s
    ## 126500K .......... .......... .......... .......... .......... 94%  104M 0s
    ## 126550K .......... .......... .......... .......... .......... 94% 43.7M 0s
    ## 126600K .......... .......... .......... .......... .......... 94% 84.3M 0s
    ## 126650K .......... .......... .......... .......... .......... 94% 77.6M 0s
    ## 126700K .......... .......... .......... .......... .......... 94%  118M 0s
    ## 126750K .......... .......... .......... .......... .......... 94% 98.8M 0s
    ## 126800K .......... .......... .......... .......... .......... 94% 91.8M 0s
    ## 126850K .......... .......... .......... .......... .......... 94%  101M 0s
    ## 126900K .......... .......... .......... .......... .......... 94% 72.1M 0s
    ## 126950K .......... .......... .......... .......... .......... 94% 89.9M 0s
    ## 127000K .......... .......... .......... .......... .......... 94% 93.1M 0s
    ## 127050K .......... .......... .......... .......... .......... 94%  102M 0s
    ## 127100K .......... .......... .......... .......... .......... 94%  118M 0s
    ## 127150K .......... .......... .......... .......... .......... 94% 47.0M 0s
    ## 127200K .......... .......... .......... .......... .......... 94% 83.4M 0s
    ## 127250K .......... .......... .......... .......... .......... 94% 89.9M 0s
    ## 127300K .......... .......... .......... .......... .......... 94%  137M 0s
    ## 127350K .......... .......... .......... .......... .......... 95%  100M 0s
    ## 127400K .......... .......... .......... .......... .......... 95%  103M 0s
    ## 127450K .......... .......... .......... .......... .......... 95% 74.5M 0s
    ## 127500K .......... .......... .......... .......... .......... 95% 85.7M 0s
    ## 127550K .......... .......... .......... .......... .......... 95% 88.2M 0s
    ## 127600K .......... .......... .......... .......... .......... 95%  107M 0s
    ## 127650K .......... .......... .......... .......... .......... 95%  141M 0s
    ## 127700K .......... .......... .......... .......... .......... 95% 37.0M 0s
    ## 127750K .......... .......... .......... .......... .......... 95% 91.1M 0s
    ## 127800K .......... .......... .......... .......... .......... 95%  104M 0s
    ## 127850K .......... .......... .......... .......... .......... 95%  122M 0s
    ## 127900K .......... .......... .......... .......... .......... 95%  109M 0s
    ## 127950K .......... .......... .......... .......... .......... 95%  102M 0s
    ## 128000K .......... .......... .......... .......... .......... 95% 38.7M 0s
    ## 128050K .......... .......... .......... .......... .......... 95% 74.7M 0s
    ## 128100K .......... .......... .......... .......... .......... 95% 92.1M 0s
    ## 128150K .......... .......... .......... .......... .......... 95% 93.3M 0s
    ## 128200K .......... .......... .......... .......... .......... 95%  118M 0s
    ## 128250K .......... .......... .......... .......... .......... 95%  135M 0s
    ## 128300K .......... .......... .......... .......... .......... 95%  150M 0s
    ## 128350K .......... .......... .......... .......... .......... 95% 72.1M 0s
    ## 128400K .......... .......... .......... .......... .......... 95% 77.4M 0s
    ## 128450K .......... .......... .......... .......... .......... 95% 96.1M 0s
    ## 128500K .......... .......... .......... .......... .......... 95% 94.8M 0s
    ## 128550K .......... .......... .......... .......... .......... 95%  118M 0s
    ## 128600K .......... .......... .......... .......... .......... 95% 55.2M 0s
    ## 128650K .......... .......... .......... .......... .......... 95%  111M 0s
    ## 128700K .......... .......... .......... .......... .......... 96%  126M 0s
    ## 128750K .......... .......... .......... .......... .......... 96% 84.4M 0s
    ## 128800K .......... .......... .......... .......... .......... 96% 97.0M 0s
    ## 128850K .......... .......... .......... .......... .......... 96%  126M 0s
    ## 128900K .......... .......... .......... .......... .......... 96% 37.2M 0s
    ## 128950K .......... .......... .......... .......... .......... 96% 86.3M 0s
    ## 129000K .......... .......... .......... .......... .......... 96%  112M 0s
    ## 129050K .......... .......... .......... .......... .......... 96%  109M 0s
    ## 129100K .......... .......... .......... .......... .......... 96%  152M 0s
    ## 129150K .......... .......... .......... .......... .......... 96% 78.6M 0s
    ## 129200K .......... .......... .......... .......... .......... 96% 71.7M 0s
    ## 129250K .......... .......... .......... .......... .......... 96% 92.8M 0s
    ## 129300K .......... .......... .......... .......... .......... 96% 86.9M 0s
    ## 129350K .......... .......... .......... .......... .......... 96% 81.3M 0s
    ## 129400K .......... .......... .......... .......... .......... 96% 97.7M 0s
    ## 129450K .......... .......... .......... .......... .......... 96%  139M 0s
    ## 129500K .......... .......... .......... .......... .......... 96% 50.4M 0s
    ## 129550K .......... .......... .......... .......... .......... 96% 89.1M 0s
    ## 129600K .......... .......... .......... .......... .......... 96%  110M 0s
    ## 129650K .......... .......... .......... .......... .......... 96%  116M 0s
    ## 129700K .......... .......... .......... .......... .......... 96%  150M 0s
    ## 129750K .......... .......... .......... .......... .......... 96% 53.0M 0s
    ## 129800K .......... .......... .......... .......... .......... 96% 96.9M 0s
    ## 129850K .......... .......... .......... .......... .......... 96% 89.6M 0s
    ## 129900K .......... .......... .......... .......... .......... 96%  106M 0s
    ## 129950K .......... .......... .......... .......... .......... 96% 75.7M 0s
    ## 130000K .......... .......... .......... .......... .......... 97%  106M 0s
    ## 130050K .......... .......... .......... .......... .......... 97% 47.0M 0s
    ## 130100K .......... .......... .......... .......... .......... 97% 78.0M 0s
    ## 130150K .......... .......... .......... .......... .......... 97% 74.2M 0s
    ## 130200K .......... .......... .......... .......... .......... 97% 93.2M 0s
    ## 130250K .......... .......... .......... .......... .......... 97%  148M 0s
    ## 130300K .......... .......... .......... .......... .......... 97%  127M 0s
    ## 130350K .......... .......... .......... .......... .......... 97% 67.2M 0s
    ## 130400K .......... .......... .......... .......... .......... 97%  111M 0s
    ## 130450K .......... .......... .......... .......... .......... 97% 68.8M 0s
    ## 130500K .......... .......... .......... .......... .......... 97% 82.7M 0s
    ## 130550K .......... .......... .......... .......... .......... 97% 89.7M 0s
    ## 130600K .......... .......... .......... .......... .......... 97%  141M 0s
    ## 130650K .......... .......... .......... .......... .......... 97% 62.9M 0s
    ## 130700K .......... .......... .......... .......... .......... 97% 79.2M 0s
    ## 130750K .......... .......... .......... .......... .......... 97% 73.5M 0s
    ## 130800K .......... .......... .......... .......... .......... 97% 97.9M 0s
    ## 130850K .......... .......... .......... .......... .......... 97%  113M 0s
    ## 130900K .......... .......... .......... .......... .......... 97% 81.0M 0s
    ## 130950K .......... .......... .......... .......... .......... 97% 87.9M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 66.9M 0s
    ## 131050K .......... .......... .......... .......... .......... 97% 96.7M 0s
    ## 131100K .......... .......... .......... .......... .......... 97%  111M 0s
    ## 131150K .......... .......... .......... .......... .......... 97%  112M 0s
    ## 131200K .......... .......... .......... .......... .......... 97% 44.5M 0s
    ## 131250K .......... .......... .......... .......... .......... 97%  111M 0s
    ## 131300K .......... .......... .......... .......... .......... 97%  109M 0s
    ## 131350K .......... .......... .......... .......... .......... 98% 87.5M 0s
    ## 131400K .......... .......... .......... .......... .......... 98% 88.4M 0s
    ## 131450K .......... .......... .......... .......... .......... 98% 88.5M 0s
    ## 131500K .......... .......... .......... .......... .......... 98%  119M 0s
    ## 131550K .......... .......... .......... .......... .......... 98% 57.9M 0s
    ## 131600K .......... .......... .......... .......... .......... 98% 88.3M 0s
    ## 131650K .......... .......... .......... .......... .......... 98%  113M 0s
    ## 131700K .......... .......... .......... .......... .......... 98%  119M 0s
    ## 131750K .......... .......... .......... .......... .......... 98% 39.4M 0s
    ## 131800K .......... .......... .......... .......... .......... 98% 94.8M 0s
    ## 131850K .......... .......... .......... .......... .......... 98%  101M 0s
    ## 131900K .......... .......... .......... .......... .......... 98% 87.7M 0s
    ## 131950K .......... .......... .......... .......... .......... 98% 84.2M 0s
    ## 132000K .......... .......... .......... .......... .......... 98%  100M 0s
    ## 132050K .......... .......... .......... .......... .......... 98%  114M 0s
    ## 132100K .......... .......... .......... .......... .......... 98%  101M 0s
    ## 132150K .......... .......... .......... .......... .......... 98% 63.6M 0s
    ## 132200K .......... .......... .......... .......... .......... 98% 94.7M 0s
    ## 132250K .......... .......... .......... .......... .......... 98%  106M 0s
    ## 132300K .......... .......... .......... .......... .......... 98%  130M 0s
    ## 132350K .......... .......... .......... .......... .......... 98% 45.2M 0s
    ## 132400K .......... .......... .......... .......... .......... 98% 93.7M 0s
    ## 132450K .......... .......... .......... .......... .......... 98%  121M 0s
    ## 132500K .......... .......... .......... .......... .......... 98%  102M 0s
    ## 132550K .......... .......... .......... .......... .......... 98% 82.7M 0s
    ## 132600K .......... .......... .......... .......... .......... 98% 96.0M 0s
    ## 132650K .......... .......... .......... .......... .......... 98% 97.0M 0s
    ## 132700K .......... .......... .......... .......... .......... 99%  102M 0s
    ## 132750K .......... .......... .......... .......... .......... 99% 69.8M 0s
    ## 132800K .......... .......... .......... .......... .......... 99%  141M 0s
    ## 132850K .......... .......... .......... .......... .......... 99%  102M 0s
    ## 132900K .......... .......... .......... .......... .......... 99% 42.3M 0s
    ## 132950K .......... .......... .......... .......... .......... 99% 80.4M 0s
    ## 133000K .......... .......... .......... .......... .......... 99% 83.7M 0s
    ## 133050K .......... .......... .......... .......... .......... 99%  125M 0s
    ## 133100K .......... .......... .......... .......... .......... 99%  113M 0s
    ## 133150K .......... .......... .......... .......... .......... 99% 84.4M 0s
    ## 133200K .......... .......... .......... .......... .......... 99%  113M 0s
    ## 133250K .......... .......... .......... .......... .......... 99% 58.3M 0s
    ## 133300K .......... .......... .......... .......... .......... 99% 76.8M 0s
    ## 133350K .......... .......... .......... .......... .......... 99%  102M 0s
    ## 133400K .......... .......... .......... .......... .......... 99%  120M 0s
    ## 133450K .......... .......... .......... .......... .......... 99%  149M 0s
    ## 133500K .......... .......... .......... .......... .......... 99% 50.0M 0s
    ## 133550K .......... .......... .......... .......... .......... 99% 76.4M 0s
    ## 133600K .......... .......... .......... .......... .......... 99%  122M 0s
    ## 133650K .......... .......... .......... .......... .......... 99%  129M 0s
    ## 133700K .......... .......... .......... .......... .......... 99%  110M 0s
    ## 133750K .......... .......... .......... .......... .......... 99% 95.2M 0s
    ## 133800K .......... .......... .......... .......... .......... 99% 40.3M 0s
    ## 133850K .......... .......... .......... .......... .......... 99%  101M 0s
    ## 133900K .......... .......... .......... .......... .......... 99%  111M 0s
    ## 133950K .......... .......... .......... .......... .......... 99%  119M 0s
    ## 134000K .......... .......... .......... .......... .......... 99%  146M 0s
    ## 134050K .......... .....                                      100%  129M=2.3s
    ## 
    ## 2023-01-07 19:11:58 (56.1 MB/s) - ‘silva_nr99_v138.1_train_set.fa.gz?download=1.3’ saved [137283333/137283333]

``` r
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz?download=1", multithread=TRUE)
```

``` r
#Inspecter les affectations taxonomiques
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order               
    ## [1,] "Bacteria" NA             NA            NA                  
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales"     
    ## [3,] "Bacteria" "Firmicutes"   "Bacilli"     "Erysipelotrichales"
    ## [4,] "Bacteria" "Firmicutes"   "Clostridia"  "Clostridiales"     
    ## [5,] "Bacteria" "Firmicutes"   "Clostridia"  "Oscillospirales"   
    ## [6,] "Bacteria" "Firmicutes"   "Clostridia"  "Oscillospirales"   
    ##      Family                                  Genus                        
    ## [1,] NA                                      NA                           
    ## [2,] "F082"                                  NA                           
    ## [3,] "Erysipelatoclostridiaceae"             "Coprobacillus"              
    ## [4,] "Clostridiaceae"                        "Clostridium sensu stricto 1"
    ## [5,] "Oscillospiraceae"                      "UCG-005"                    
    ## [6,] "[Eubacterium] coprostanoligenes group" NA

``` r
library(phyloseq); packageVersion("phyloseq")
```

    ## [1] '1.40.0'

``` r
library(Biostrings); packageVersion("Biostrings")
```

    ## [1] '2.64.1'

``` r
library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '3.3.6'

``` r
theme_set(theme_bw())
```

``` r
#Les donnees de sequencage mise en ligne n'etaient pas separees selon les differents regions ou l'etude a ete effectue, donc on ne peut pas etablir un graphique montrant les differents microorganismes selon la region.
```
