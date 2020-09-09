tidySE - part of tidytranscriptomics
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidySE/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/stemangiola/tidySE/actions)
<!-- badges: end -->

**Brings SummarizedExperiment to the tidyverse\!**

website:
[stemangiola.github.io/tidySE/](https://stemangiola.github.io/tidySE/)

Please have a look also to

  - [tidySCE](https://stemangiola.github.io/tidySCE/) for tidy
    manipulation of Seurat objects
  - [tidyseurat](https://stemangiola.github.io/tidyseurat/) for tidy
    manipulation of Seurat objects
  - [tidybulk](https://stemangiola.github.io/tidybulk/) for tidy
    high-level data analysis and manipulation
  - [nanny](https://github.com/stemangiola/nanny) for tidy high-level
    data analysis and manipulation
  - [tidygate](https://github.com/stemangiola/tidygate) for adding
    custom gate information to your tibble
  - [tidyHeatmap](https://stemangiola.github.io/tidyHeatmap/) for
    heatmaps produced with tidy principles

<!---

[![Build Status](https://travis-ci.org/stemangiola/tidySE.svg?branch=master)](https://travis-ci.org/stemangiola/tidySE) [![Coverage Status](https://coveralls.io/repos/github/stemangiola/tidySE/badge.svg?branch=master)](https://coveralls.io/github/stemangiola/tidySE?branch=master)

-->

## Functions/utilities available

| SummarizedExperiment-compatible Functions | Description                                                      |
| ----------------------------------------- | ---------------------------------------------------------------- |
| `all`                                     | After all `tidySE` is a SummarizedExperiment object, just better |

| tidyverse Packages | Description                          |
| ------------------ | ------------------------------------ |
| `dplyr`            | All `dplyr` APIs like for any tibble |
| `tidyr`            | All `tidyr` APIs like for any tibble |
| `ggplot2`          | `ggplot` like for any tibble         |
| `plotly`           | `plot_ly` like for any tibble        |

| Utilities   | Description                                                     |
| ----------- | --------------------------------------------------------------- |
| `tidy`      | Add `tidySE` invisible layer over a SummarizedExperiment object |
| `as_tibble` | Convert cell-wise information to a `tbl_df`                     |

## Installation

From Github

``` r
devtools::install_github("stemangiola/tidySE")
```

## Create `tidySE`, the best of both worlds\!

This is a SummarizedExperiment object but it is evaluated as tibble. So
it is fully compatible both with SummarizedExperiment and tidyverse
APIs.

``` r
pasilla_tidy = tidySE::pasilla %>% tidySE::tidy()
```

**It looks like a tibble**

``` r
pasilla_tidy
```

    ## # A tibble: 102,193 x 5
    ##    sample condition type       transcript  counts
    ##    <chr>  <chr>     <chr>      <chr>        <int>
    ##  1 untrt1 untreated single_end FBgn0000003      0
    ##  2 untrt1 untreated single_end FBgn0000008     92
    ##  3 untrt1 untreated single_end FBgn0000014      5
    ##  4 untrt1 untreated single_end FBgn0000015      0
    ##  5 untrt1 untreated single_end FBgn0000017   4664
    ##  6 untrt1 untreated single_end FBgn0000018    583
    ##  7 untrt1 untreated single_end FBgn0000022      0
    ##  8 untrt1 untreated single_end FBgn0000024     10
    ##  9 untrt1 untreated single_end FBgn0000028      0
    ## 10 untrt1 untreated single_end FBgn0000032   1446
    ## # … with 102,183 more rows

**But it is a SummarizedExperiment object after all**

``` r
pasilla_tidy@assays
```

    ## An object of class "SimpleAssays"
    ## Slot "data":
    ## List of length 1
    ## names(1): counts

## Annotation polishing using tidyverse

We may have a column that contains the directory each run was taken
from. We may want to extract the run/sample name out of it.

``` r
pasilla_polished =
  pasilla_tidy %>%
  mutate(type = gsub("_end", "", type)) 

pasilla_polished 
```

    ## # A tibble: 102,193 x 5
    ##    sample condition type   transcript  counts
    ##    <chr>  <chr>     <chr>  <chr>        <int>
    ##  1 untrt1 untreated single FBgn0000003      0
    ##  2 untrt1 untreated single FBgn0000008     92
    ##  3 untrt1 untreated single FBgn0000014      5
    ##  4 untrt1 untreated single FBgn0000015      0
    ##  5 untrt1 untreated single FBgn0000017   4664
    ##  6 untrt1 untreated single FBgn0000018    583
    ##  7 untrt1 untreated single FBgn0000022      0
    ##  8 untrt1 untreated single FBgn0000024     10
    ##  9 untrt1 untreated single FBgn0000028      0
    ## 10 untrt1 untreated single FBgn0000032   1446
    ## # … with 102,183 more rows

## Preliminary plots

We can treat `pasilla_polished` effectively as a normal tibble for
plotting.

Here we plot the distribution of counts per sample

``` r
pasilla_polished %>%
    ggplot(aes(counts + 1, group=sample, color=`type`)) +
    geom_density() +
    scale_x_log10() +
    my_theme
```

![](man/figures/plot1-1.png)<!-- -->

## Nested analyses

A powerful tool we can use with tidySE is `nest`. We can easily perform
independent analyses on subsets of the dataset. First we classify cell
types in lymphoid and myeloid; then, nest based on the new
classification

``` r
pasilla_nested = 
  pasilla_polished %>%
  nest(data = -type)

pasilla_nested
```

    ## # A tibble: 2 x 2
    ##   type   data    
    ##   <chr>  <list>  
    ## 1 single <tidySE>
    ## 2 paired <tidySE>
