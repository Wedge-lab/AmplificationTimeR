# AmplificationTimeR
<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/mjakobs/AmplificationTimeR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mjakobs/AmplificationTimeR?branch=main)
  [![R](https://github.com/mjakobs/AmplificationTimeR/actions/workflows/r.yml/badge.svg)](https://github.com/mjakobs/AmplificationTimeR/actions/workflows/r.yml)
  <!-- badges: end -->

`AmplificationTimeR` is an R package for timing individual amplification and whole genome duplication events affecting regions regions of the tumour genome that have been focally amplified.  `AmplificationTimeR` can be thought of as an extension to `gerstung-lab/MutationTimeR` with slightly different intentions and functionality.  

To install AmplificationTimeR run:
```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("mjakobs/AmplificationTimeR")
```
