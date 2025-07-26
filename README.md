
# lit

<!-- badges: start -->

[![R-CMD-check](https://github.com/ajbass/lit/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ajbass/lit/actions/workflows/R-CMD-check.yml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/lit)](https://cran.r-project.org/package=lit)

<!-- badges: end -->

## Overview

The `lit` package implements a kernel-based multivariate testing
procedure, called Latent Interaction Testing (LIT), to test for latent
genetic interactions in genome-wide association studies. See our
manuscript for additional details:

> Bass AJ, Bian S, Wingo AP, Wingo TS, Culter DJ, Epstein MP.
> Identifying latent genetic interactions in genome-wide association
> studies using multiple traits. Genome Medicine; 2024.

## Installation

This software is implemented in the `R` statistical programming
language. To install the release version, type the following in the `R`
terminal:

``` r
# release version
install.packages("lit")
```

The development version of `lit` can be installed using the following
code:

``` r
# install devtools
install.packages("devtools")
devtools::install_github("ajbass/lit")
```

The vignette can be viewed by typing:

``` r
browseVignettes(package = "lit")
```

If you run into issues with `gfortran` on Mac, see the answer
[here](https://stackoverflow.com/questions/69639782/installing-gfortran-on-macbook-with-apple-m1-chip-for-use-in-r)
for additional details.

## Quick start

We provide two ways to use the `lit` package. When the genotypes can be
loaded in R (small GWAS datasets), the `lit()` function can be used:

``` r
library(lit)
# set seed
set.seed(123)

# generate 10 SNPs for 10 individuals
X <- matrix(rbinom(10 * 10, size = 2, prob = 0.25), ncol = 10)

# generate 4 phenotypes for 10 individuals
Y <- matrix(rnorm(10 * 4), ncol = 4) 

# test for latent genetic interactions
out <- lit(Y, X)
head(out)
#>        wlit      ulit      alit
#> 1 0.2681410 0.3504852 0.3056363
#> 2 0.7773637 0.3504852 0.6044655
#> 3 0.4034423 0.3504852 0.3760632
#> 4 0.7874949 0.3504852 0.6157108
#> 5 0.8701189 0.3504852 0.7337565
#> 6 0.2352616 0.3504852 0.2847600
```

The output is a data frame of p-values where the rows are SNPs and the
columns are different implementations of LIT to test for latent genetic
interactions:

- `wlit` uses a linear kernel to measure pairwise similarity for the
  genotype and trait matrices
- `ulit` uses a projection kernel to measure pairwise similarity for the
  genotype and trait matrices
- `alit` combines the p-values of `wlit` and `ulit` using a Cauchy
  combination test to maximize the number of discoveries

For large GWAS datasets (e.g., biobank-sized), the `lit()` function is
not computationally feasible because the genotypes cannot be loaded in
`R`. Instead, the `lit_plink()` function can be applied directly to
plink files. To demonstrate how to use the function, we use the example
plink files from the `genio` package:

``` r
# load genio package
library(genio)

# path to plink files
file <- system.file("extdata", 'sample.bed', package = "genio", mustWork = TRUE)

# generate trait expression
Y <- matrix(rnorm(10 * 4), ncol = 4)

# apply lit to plink file
out <- lit_plink(Y, file = file, verbose = FALSE)
head(out)
#>   chr         id     pos alt ref       maf      wlit      ulit      alit
#> 1   1  rs3094315  752566   G   A 0.3888889 0.7908763 0.3422960 0.6150572
#> 2   1  rs7419119  842013   T   G 0.3888889 0.1552580 0.3422960 0.2194972
#> 3   1 rs13302957  891021   G   A 0.2500000 0.4088937 0.3325939 0.3687589
#> 4   1  rs6696609  903426   C   T 0.3125000 0.5857829 0.3325939 0.4519475
#> 5   1     rs8997  949654   A   G 0.4375000 0.6628300 0.3325939 0.4969663
#> 6   1  rs9442372 1018704   A   G 0.2500000 0.3192430 0.3325939 0.3258332
```

See `?lit` and `?lit_plink` for additional details and input arguments.

Note that a marginal testing procedure for latent genetic interactions
based on the squared residuals and cross products (Marginal (SQ/CP)) can
also be implemented using the `marginal` and `marginal_plink` functions:

``` r
# apply Marginal (SQ/CP) to loaded genotypes
out <- marginal(Y, X)

# apply Marginal (SQ/CP) to plink file
out <- marginal_plink(Y, file = file, verbose = FALSE)
```

<!---
The columns represent the cross products, e.g., `col1_col2` is the cross product of the first trait multiplied by the second trait (adjusting for the additive genetic effect)). See our preprint for additional details on this approach.

## To do

List of vario
NA's in phenotypes.

-->
