addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%% set up simulation
clc; clear all; close all;
rng(3)

sim.seed = randi(10);
sim.T = 10*60;
sim.dt = 1e-3;
sim.vecN = round(sim.T/sim.dt);
sim.pPreSpike = 5*sim.dt;
sim.alpha_dt = 0.004;
sim.alpha_tau = 0.001;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1;
sim.hist_tau = .01;
sim.hist_beta = -1;
sim.stp_B = [0 1 2 3 4]'*-.04;
sim.wt_long = wtLongAdjust(1.5, sim);
sim.beta0 = beta0Adjust(20, sim);
data.dt = sim.dt;
[data,sim] = sim_model(data,sim);

%% fix beta0

% over
[fitBeta0_over,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_beta0', sim.beta0 + 1,...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

% under
[fitBeta0_under,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_beta0', sim.beta0 - 1,...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

% true
[fitBeta0_true,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_beta0', sim.beta0,...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

%% fix wt_long

% over
[fitwtLong_over,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_wt_long', sim.wt_long + 1,...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

% under
[fitwtLong_under,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_wt_long', sim.wt_long - 1,...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

% true
[fitwtLong_true,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_wt_long', sim.wt_long - 1,...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

%% fix wt_short

% over

% under

% true











