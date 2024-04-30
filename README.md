# Adaptive Catalyst Discovery Using Multicriteria Bayesian Optimization with Representation Learning

# Overview

In this package, we provide a high-throughput computational catalyst screening approach integrating density functional theory (DFT) and Bayesian Optimization (BO). Within the BO framework, we propose an uncertainty-aware atomistic machine learning model, UPNet, which enables automated representation learning directly from high-dimensional catalyst structures and achieves principled uncertainty quantification. Utilizing a constrained expected improvement acquisition function, our BO framework simultaneously considers multiple evaluation criteria. Using the proposed methods, we explore catalyst discovery for the CO2 reduction reaction. 


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
```
Python 3.9.7
pandas 1.4.3
numpy 1.23.1
scipy 1.4.1
matplotlib 3.5.2
tensorflow 2.2.0
```


# Installation Guide and Instructions for Use

Please download all files and run "bayesian optimization CO2R.py".

The code takes approximately 60 seconds before running BO iterations. 

# Demo

The total number of BO iterations, "sampling_itr," is set to 3 for the sake of demonstration brevity. In the paper, the number of iterations is 80.

# Results

The output is the initial samples and samples selected by BO. The results are plotted on the volcano scaling relationship for activity.

