addpath(genpath('D:/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('D:/GitHub/GBLM_SMOOTH/core'));

%% pre-checking

isiDisCheck(T, dt, pPreSpike, seed);
isiCheck(trueParam, Nq_true, Nm, Nst);
lamCehck(T, dt, pPreSpike, beta0, wt_long,t_alpha, tau_alpha, trueParam, Nq_true, Nm, Nst, seed);
llhdCheck(T, dt, pPreSpike, beta0, wt_long, t_alpha, tau_alpha, trueParam, Nq_true, Nm, Nst, seed);


%% post-checking
clc;clear all;close all;
saveDir = 'C:/Users/gaw19004/Desktop/ganchaoResearch/results/new/cases/jump/';
load(fullfile(saveDir, '6_noPlasticity.mat'));


elpsT/60
static.tAlpha
static.tauAlpha

plot(dev(2:k));
plot(llhd(2:k));

e = llhdCheck2(llhd, k, T, dt, pPreSpike, beta0,...
    wt_long, t_alpha, tau_alpha, trueParam, Nq_true, Nm, Nst, seed);

simPlot(beta0_final, W_beta0_final, beta0,...
    wt_long_final, W_wt_long_final,...
    wt_long, wt_short_final, wt_short, covB_final, e);

isiCheck2(Nq_true, Nq_fit, Nm, Nst, trueParam, wt_short_param_final,covB_final)