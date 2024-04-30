# Adaptive Catalyst Discovery Using Multicriteria Bayesian Optimization with Representation Learning

# Overview

High-performance catalysts are crucial for sustainable energy conversion and human health. However, the discovery of catalysts faces challenges due to the absence of efficient approaches to navigating vast and high-dimensional structure and composition spaces. In this package, we provide a high-throughput computational catalyst screening approach integrating density functional theory (DFT) and Bayesian Optimization (BO). Within the BO framework, we propose an uncertainty-aware atomistic machine learning model, UPNet, which enables automated representation learning directly from high-dimensional catalyst structures and achieves principled uncertainty quantification. Utilizing a constrained expected improvement acquisition function, our BO framework simultaneously considers multiple evaluation criteria. Using the proposed methods, we explore catalyst discovery for the CO2 reduction reaction. 


# System Requirements

## Hardware Requirements

The package development version is tested on a computer with the following specs:

CPU: Intel(R) Xeon(R) W-2295  @ 3.00GHz
GPU: NVIDIA TU104GL

## Software Requirements

### OS Requirements

The package development version is tested on *Linux* operating systems. The developmental version of the package has been tested on the following systems:

Linux:  5.4.0-176-generic #196-Ubuntu SMP  

### Software Dependencies
Python 3.9.7
pandas 1.4.3
numpy 1.23.1
scipy 1.4.1
matplotlib 3.5.2
tensorflow 2.2.0



# Installation Guide

## Stable Release

`lolR` is available in a stable release on CRAN:

```
install.packages('lolR')
```

## Development Version

### Package dependencies

Users should install the following packages prior to installing `lolR`, from an `R` terminal:

```
install.packages(c('ggplot2', 'abind', 'irlba', 'knitr', 'rmarkdown', 'latex2exp', 'MASS', 'randomForest'))
```

which will install in about 30 seconds on a machine with the recommended specs.

The `lolR` package functions with all packages in their latest versions as they appear on `CRAN` on December 13, 2017. Users can check [CRAN snapshot](https://mran.microsoft.com/timemachine/) for details. The versions of software are, specifically:
```
abind_1.4-5
latex2exp_0.4.0
ggplot2_2.2.1
irlba_2.3.1
Matrix_1.2-3
MASS_7.3-47
randomForest_4.6-12
```

If you are having an issue that you believe to be tied to software versioning issues, please drop us an [Issue](https://github.com/neurodata/lol/issues). 

### Package Installation

From an `R` session, type:

```
require(devtools)
install_github('neurodata/lol', build_vignettes=TRUE, force=TRUE)  # install lol with the vignettes
require(lolR)
vignette("lol", package="lolR")  # view one of the basic vignettes
```

The package should take approximately 40 seconds to install with vignettes on a recommended computer. 

# Demo

## Functions

For interactive demos of the functions, please check out the vignettes built into the package. They can be accessed as follows:

```
require(lolR)
vignette('lol')
vignette('pca')
vignette('cpca')
vignette('lrcca')
vignette('mdp')
vignette('xval')
vignette('qoq')
vignette('simulations')
vignette('nearestCentroid')
```

## Extending the lolR Package

The lolR package makes many useful resources available (such as embedding and cross-validation) for simple extension. 

To extend the lolR package, check out the vignettes:

```
require(lolR)
vignette('extend_embedding')
vignette('extend_classification')
```

# Results

In this [benchmark comparison](http://docs.neurodata.io/lol/lol-paper/figures/real_data.html), we show that LOL does better than all linear embedding techniques in supervised HDLSS settings when dimensionality is high (d > 100, ntrain <= d) on 20 benchmark problems from the [UCI](https://archive.ics.uci.edu/ml/index.php) and [PMLB](https://github.com/EpistasisLab/penn-ml-benchmarks) datasets. LOL provides a good tradeoff between maintaining the class conditional difference (good misclassification rate) in a small number of dimensions (low number of embedding dimensions).

# Citation

For usage of the package and associated manuscript, please cite according to the enclosed [citation.bib](./citation.bib).
