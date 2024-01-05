# AmplificationTimeR
<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/Wedge-lab/AmplificationTimeR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Wedge-lab/AmplificationTimeR?branch=main)
  [![R-CMD-check](https://github.com/Wedge-lab/AmplificationTimeR/actions/workflows/check-standard.yml/badge.svg)](https://github.com/Wedge-lab/AmplificationTimeR/actions/workflows/check-standard.yml)
  <!-- badges: end -->

- [AmplificationTimeR](#amplificationtimer)
  - [Citation](#citation)
  - [Introduction](#introduction)
  - [Installation](#installation)
  - [Input data](#input-data)
    - [Copy number data](#copy-number-data)
    - [Multiplicity data](#multiplicity-data)
    - [Mutation data](#mutation-data)
    - [A Note on Mutation Types](#a-note-on-mutation-types)
    - [Coordinates for a region of interest](#coordinates-for-a-region-of-interest)
    - [WGD status](#wgd-status)
    - [Reference genome](#reference-genome)
    - [Example](#example)
    - [Output flags](#output-flags)


## Citation
If you use `AmplificationTimeR` in your work, please cite our manuscript (currently under review - details TBC).  

## Introduction
`AmplificationTimeR` is an R package for timing individual amplification and whole genome duplication events affecting regions of the tumour genome that have been focally amplified.  

## Installation
`AmplificationTimeR` requires the following packages to be installed before use:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("SomaticSignatures", quietly = TRUE))
    BiocManager::install("SomaticSignatures")
if (!require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
    # OR
if (!require("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")    
```

To install `AmplificationTimeR` run:
```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("Wedge-lab/AmplificationTimeR")
```

## Input data
For demonstration purposes we have included some data with the `AmplificationTimeR` package:
```r
data(demo_cn, demo_mult, demo_muts)
```

### Copy number data
Copy number data should be derived from [`Battenberg`](https://github.com/Wedge-lab/battenberg) and must contain the following columns at the very least: `"chr"`, `"startpos"`, `"endpos"`, `"nMaj1_A"`, `"nMin1_A"`. `"nMaj2_A"` and `"nMin2_A"` should be included if applicable.  Other columns may be present but the 5 columns listed above are a core requirement. 
```r
demo_cn

##   chr startpos endpos nMaj1_A nMin1_A
## 1   1    50000  1e+06       5       2
```

### Multiplicity data
Multiplicity data should be derived from [`dpclust3p`](https://github.com/Wedge-lab/dpclust3p). The multiplicity data frame must have the following columns at a minimum: `"chr"`, `"end"`, and `"no.chrs.bearing.mut"`. Other columns may be present. 

```r
head(demo_mult,20)

##    chr    end no.chrs.bearing.mut
## 1    1  57682                   1
## 2    1  62417                   1
## 3    1  67153                   1
## 4    1  71888                   1
## 5    1  76623                   1
## 6    1  81358                   1
## 7    1  86094                   1
## 8    1  90829                   1
## 9    1  95564                   1
## 10   1 100299                   1
## 11   1 105035                   1
## 12   1 109770                   1
## 13   1 114505                   1
## 14   1 119240                   1
## 15   1 123976                   1
## 16   1 128711                   1
## 17   1 133446                   1
## 18   1 138182                   1
## 19   1 142917                   1
## 20   1 147652                   1
```
```r
str(demo_mult)
 
## 'data.frame':	200 obs. of  3 variables:
##  $ chr                : num  1 1 1 1 1 1 1 1 1 1 ...
##  $ end                : num  57682 62417 67153 71888 76623 ...
##  $ no.chrs.bearing.mut: num  1 1 1 1 1 1 1 1 1 1 ...
```

### Mutation data
The mutation data frame should contain either all mutations present in the sample or region, or only mutations attributed to mutational signatures SBS1 and SBS5. The presence of the following columns is required:`"chr"`, `"start"`, `"end"`, `"ref"`, `"alt"`. Other columns may be present. 

```r
head(demo_muts, 20)

##    chr  start    end ref alt
## 1    1  57682  57682   C   T
## 2    1  62417  62417   C   T
## 3    1  67153  67153   C   T
## 4    1  71888  71888   C   T
## 5    1  76623  76623   A   G
## 6    1  81358  81358   C   T
## 7    1  86094  86094   C   T
## 8    1  90829  90829   C   T
## 9    1  95564  95564   C   T
## 10   1 100299 100299   A   G
## 11   1 105035 105035   C   T
## 12   1 109770 109770   C   T
## 13   1 114505 114505   C   T
## 14   1 119240 119240   C   T
## 15   1 123976 123976   A   G
## 16   1 128711 128711   C   T
## 17   1 133446 133446   C   T
## 18   1 138182 138182   C   T
## 19   1 142917 142917   C   T
## 20   1 147652 147652   A   G
```
```r
str(demo_muts)

## 'data.frame':	200 obs. of  5 variables:
##  $ chr  : num  1 1 1 1 1 1 1 1 1 1 ...
##  $ start: num  57682 62417 67153 71888 76623 ...
##  $ end  : num  57682 62417 67153 71888 76623 ...
##  $ ref  : chr  "C" "C" "C" "C" ...
##  $ alt  : chr  "T" "T" "T" "T" ...
```

### A Note on Mutation Types
For best results we recommend restricting analysis to clock-like mutations, such as mutations attributed to mutational signatures SBS1 and SBS5, or C>T mutations at CpG sites.  If the `muts_type` argument is set to `muts_type="All"`, as is the default,`AmplificationTimeR` will subset the provided list of mutations for C>T mutations occurring at CpG sites using the `mutationContext` function from Julian Gehring's [`SomaticSignatures`](https://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html) Bioconductor package. If `muts_type="SBS1 and SBS5"`, a data frame containing only mutations attributed to mutational signatures SBS1 and SBS5 can be provided.  If `muts_type="SBS1 and SBS5"` and a data frame that has not been filtered is provided, `AmplificationTimeR` will use all mutations in that data frame.  This is not recommended. 

### Coordinates for a region of interest
The numeric coordinates of a region of interest.  

```r
amp_chrom <- 1
amp_start <- 50100
amp_stop <- 1500000
```

### WGD status
A `TRUE` or `FALSE` logical value encoding whether a sample has been identified as whole genome duplicated or not. 

### Reference genome
Either `"hg19"` or `"hg38"` depending on which is appropriate. This is required so that `AmplificationTimeR` can filter mutations for C>T variants occurring specifically at CpG sites, to ensure that only clock-like mutations are used. 

### Example

```r
library(AmplificationTimeR)
library(BSgenome.Hsapiens.UCSC.hg19)

data(demo_cn, demo_mult, demo_muts)

amp_chrom <- 1
amp_start <- 50100
amp_stop <- 1500000

segment_time <- time_amplification(
  cn_data = demo_cn,
  multiplicity_data = demo_mult,
  mutation_data = demo_muts,
  muts_type = "All",
  sample_id = "test_sample",
  amplification_chrom = amp_chrom,
  amplification_start = amp_start,
  amplification_stop = amp_stop,
  is_WGD = TRUE,
  genome = "hg19")
```

```r
segment_time

## Data frame with 1 row and 4 columns
##        sample        region highest_copy_number event_order num_mutations_used clonality_status flags      t_1
## 1 test_sample 1:50000-1e+06                 5+2        WGGG                 34           clonal    NA 0.893617
##   t_1_mean_bootstrap t_1_lower_ci t_1_upper_ci      t_2 t_2_mean_bootstrap t_2_lower_ci t_2_upper_ci      t_3
## 1          0.8827177    0.8723866    0.8930488 0.893617          0.8827177    0.8723866    0.8930488 0.893617
##   t_3_mean_bootstrap t_3_lower_ci t_3_upper_ci      t_4 t_4_mean_bootstrap t_4_lower_ci t_4_upper_ci t_5
## 1          0.8827177    0.8723866    0.8930488 0.893617          0.9165478    0.8942302    0.9388653  NA
##   t_5_mean_bootstrap t_5_lower_ci t_5_upper_ci t_6 t_6_mean_bootstrap t_6_lower_ci t_6_upper_ci t_7 t_7_mean_bootstrap
## 1                 NA           NA           NA  NA                 NA           NA           NA  NA                 NA
##   t_7_lower_ci t_7_upper_ci t_8 t_8_mean_bootstrap t_8_lower_ci t_8_upper_ci t_9 t_9_mean_bootstrap t_9_lower_ci
## 1           NA           NA  NA                 NA           NA           NA  NA                 NA           NA
##   t_9_upper_ci t_10 t_10_mean_bootstrap t_10_lower_ci t_10_upper_ci
## 1           NA   NA                  NA            NA            NA
```
### Output flags
We have provided a number of flags to assist users in identifying `AmplificationTimeR` results that may warrant further investigation. 

| Flag                              | Description |
| --------------------------------- | ----------- |
| "Points not in order"             | This flag indicates that one or more time points have been timed in the wrong order (e.g. `t_2_median_bootsrap` occurs after `t_3_median_bootstrap`). This may indicate that an incorrect order of events has been inferred based on the available data, potentially due to copy number losses occurring after gains. |
| "Time > 1"                        | Occasionally, `AmplificationTimeR` will provide timing estimates that exceed 1. Such time estimates may represent events that happened very late in the lifetime of the tumour. This may indicate that an incorrect order of events has been inferred based on the available data. |
| "Time < 0"                        | Occasionally, `AmplificationTimeR` will provide negative timing estimates. Such time estimates may represent a violation of the assumptions on which `AmplificationTimeR` is based. Equally, this may indicate that an incorrect order of events has been inferred based on the available data. |
| "Missing multiplicity states"     | This flag will be raised if one or more multiplicity states that are expected to be present are not identified in the segment. This situation can arise when gains happen in close succession without any mutations occurring between gain events. Alternatively, this may indicate that an incorrect order of events has been inferred from the data available. |
