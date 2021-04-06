# GBLM-SMOOTH  
GBLM-SMOOTH is a model, implemented in MATLAB, that can track long-term and short-term synaptic plasticity simultaneously. This repository contains **four** folders.  
  
## 1. "core"  
**Three** sub-folders are included:  
  
**(i) "diagnose"**: tools to evaluate fitting.  
**(ii) "modelFit"**: implementation of the model. The basic model (without tuning **_Q_** in adaptive smoothing) starts from **smooth_gblm.m**. The **_Q_**-tuned versions are generally called **tune_smooth_gblm_xd_xxxx.m**. "xd"" means whether it is 2D optimization or 1D approximation. The last four letters denote optimization methods, i.e. "_grad" for gradient descent and "_grid" for grid search. (all "v2" turn on the estimation of history effects)   
**(iii) "simulation"**: code for generating neural spikes.

## 2. "demo"
This folder provides code for generating all figures in the paper. Generally, "_caseDemo" denotes code for simulations and model fitting, while "_caseAnalysis" denotes code for plotting.

## 3. "helper"
It contains helper functions for simulations and model fitting.

## 4. "plot"
This folder contains all the plots used in the paper.







