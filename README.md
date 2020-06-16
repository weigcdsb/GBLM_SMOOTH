# GBLM_SMOOTH
GBLM-smoothing model for estimations of LTP and STP

## convergence
check convergence when considering 1) random start points and 2) order of estimations

## core
This includes 3 subfolders:

 1) diagnose: functions to help plot

 2) modelFit: 2 main functions are called smooth_gblm.m and tune_smooth_gblm.m. tune_smooth_gblm is the Q-tuned version for smooth_gblm. 

 3) simulation: functions to help simulations, such as data generation, beta0 matching and wt_long matching.

## demo
Show some examples of model fitting

## draft
drafts for paper

## fireRate
Investigate influence of pre- and post-synaptic firing rate

## helper

## interaction
Temporarily not included in the paper. Maybe we can investigate interactions by fixing 1 (true, under-estimate, over-estimate) of these 3, and then see the performance of the other 2?

## LTPUncertainty
See the influence of pre-synaptic firing rate on variance of wt_long

## Qchoose
A demo to show function tune_smooth_gblm

## refereces

## STP_type
Investigate influence of different STP (depression vs. facilitation)

## synType
Investigate influence of different synapses (excitatory vs. inhibitory)





