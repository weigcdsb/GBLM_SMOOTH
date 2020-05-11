# GBLM_SMOOTH
GBLM-smoothing model for estimations of LTP and STP

This repository contains 6 folders: helper, core, demo, absResQ, convergence and references

## helper folder
'helper' folder includes functions for alpha function and basis function

## core folder
'core' folder includes functions for the core algorithm.

## demo folder
Code in 'demo' folder is to show some examples for estimations. 

## absResQ folder
Code in 'absResQ' folder is to analysis the relationship between absolute residuals for modification (ISI-dependent) function and Q for long-term effect estimations. 
 
## convergence folder
Code in 'convergence' folder is to show the convergence of the algorithm, by 1) switching estimations of long-term effect (adaptive smoothing) and short-term effect (basis spline regression) and 2) setting random start.

## references folder
This folder contains some references

### chooseQ
It contains 3 references: 2017-- adaptive adjustment for <img src="https://render.githubusercontent.com/render/math?math=Q_k">, i.e. Q is changing at each point; 2016-1 and 2016-2-- estimations for time-constant Q.

