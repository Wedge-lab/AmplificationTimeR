# AmplificationTimeR
<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/mjakobs/AmplificationTimeR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mjakobs/AmplificationTimeR?branch=main)
  [![R](https://github.com/mjakobs/AmplificationTimeR/actions/workflows/r.yml/badge.svg)](https://github.com/mjakobs/AmplificationTimeR/actions/workflows/r.yml)
  <!-- badges: end -->

## Citation
If you use `AmplificationTimeR` in your work, please cite xyz.  

## Introduction
`AmplificationTimeR` is an R package for timing individual amplification and whole genome duplication events affecting regions of the tumour genome that have been focally amplified.  

## Installation
`AmplificationTimeR` requires the following packages to be installed before use:
```
if (!require("SomaticSignatures", quietly = TRUE))
    install.packages("SomaticSignatures")
if (!require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
    install.packages("BSgenome.Hsapiens.UCSC.hg19")
if (!require("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
    install.packages("BSgenome.Hsapiens.UCSC.hg38")    
```

To install `AmplificationTimeR` run:
```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("mjakobs/AmplificationTimeR")
```

## Usage

For best results we recommend restricting analysis to clock-like mutations, such as mutations attributed to mutational signatures SBS1 and SBS5, or C>T mutations at CpG sites.  If the `muts_type` is set to `muts_type="All`, as is the default,`AmplificationTimeR` will subset the provided list of mutations for C>T mutations occurring at CpG sites using the `mutationContext` function from Julian Gehring's [`SomaticSignatures`](https://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html) Bioconductor package. 