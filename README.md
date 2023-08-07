# AmplificationTimeR
<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/Wedge-lab/AmplificationTimeR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Wedge-lab/AmplificationTimeR?branch=main)
  [![R-CMD-check](https://github.com/Wedge-lab/AmplificationTimeR/actions/workflows/check-standard.yml/badge.svg)](https://github.com/Wedge-lab/AmplificationTimeR/actions/workflows/check-standard.yml)
  <!-- badges: end -->

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
##    chr startpos endpos nMaj1_A nMin1_A
## 1   1    50000  1e+06       2       1
```

### Multiplicity data
Multiplicity data should be derived from [`dpclust3p`](https://github.com/Wedge-lab/dpclust3p). The multiplicity data frame must have the following columns at a minimum: `"chr"`, `"end"`, and `"no.chrs.bearing.mut"`. Other columns may be present. 

```r
deom_mult

##    chr   end no.chrs.bearing.mut
## 1    1 57682                   2
## 2    1 57683                   2
## 3    1 57684                   2
## 4    1 57685                   2
## 5    1 57686                   2
## 6    1 57687                   1
## 7    1 57688                   1
## 8    1 57689                   1
## 9    1 57690                   1
## 10   1 57691                   1
## 11   1 57692                   1
## 12   1 57693                   1
## 13   1 57694                   1
## 14   1 57695                   1
## 15   1 57696                   1
## 16   1 57697                   1
## 17   1 57698                   1
## 18   1 57699                   1
## 19   1 57700                   1
## 20   1 57701                   1
```

### Mutation data
The mutation data frame should contain either all mutations present in the sample or region, or only mutations attributed to mutational signatures SBS1 and SBS5. The presence of the following columns is required:`"chr"`, `"start"`, `"end"`, `"ref"`, `"alt"`. Other columns may be present. 

```r
demo_muts

##    chr start   end ref alt
## 1    1 57682 57682   C   T
## 2    1 57683 57683   C   T
## 3    1 57684 57684   C   T
## 4    1 57685 57685   C   T
## 5    1 57686 57686   C   T
## 6    1 57687 57687   C   T
## 7    1 57688 57688   C   T
## 8    1 57689 57689   C   T
## 9    1 57690 57690   C   T
## 10   1 57691 57691   C   T
## 11   1 57692 57692   C   T
## 12   1 57693 57693   C   T
## 13   1 57694 57694   C   T
## 14   1 57695 57695   C   T
## 15   1 57696 57696   C   T
## 16   1 57697 57697   A   G
## 17   1 57698 57698   A   G
## 18   1 57699 57699   A   G
## 19   1 57700 57700   A   G
## 20   1 57701 57701   A   G
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
Either `"hg19"` or `"hg38"` depending on which is appropriate. This is required so that `AmplificationTimeR` can filter mutations for C>T variants occurring specifically at CpG sites, to ensre that only clock-like mutations are used. 

### Example

```r
library(AmplificationTimeR)
library(BSgenome.Hsapiens.UCSC.hg19)

data(demo_cn, demo_mult, demo_muts)

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

## Data frame with 1 row and 46 columns
##        sample        region highest_copy_number event_order num_mutations_used
## 1 test_sample 1:50000-1e+06                 2+1           W                  4
##   clonality_status t_1 t_1_median_bootstrap t_1_lower_ci t_1_upper_ci t_2
## 1           clonal   0                    0            0    0.5555556  NA
##   t_2_median_bootstrap t_2_lower_ci t_2_upper_ci t_3 t_3_median_bootstrap
## 1                   NA           NA           NA  NA                   NA
##   t_3_lower_ci t_3_upper_ci t_4 t_4_median_bootstrap t_4_lower_ci t_4_upper_ci t_5
## 1           NA           NA  NA                   NA           NA           NA  NA
##   t_5_median_bootstrap t_5_lower_ci t_5_upper_ci t_6 t_6_median_bootstrap
## 1                   NA           NA           NA  NA                   NA
##   t_6_lower_ci t_6_upper_ci t_7 t_7_median_bootstrap t_7_lower_ci t_7_upper_ci t_8
## 1           NA           NA  NA                   NA           NA           NA  NA
##   t_8_median_bootstrap t_8_lower_ci t_8_upper_ci t_9 t_9_median_bootstrap
## 1                   NA           NA           NA  NA                   NA
##   t_9_lower_ci t_9_upper_ci t_10 t_10_median_bootstrap t_10_lower_ci t_10_upper_ci
## 1           NA           NA   NA                    NA            NA            NA
```
