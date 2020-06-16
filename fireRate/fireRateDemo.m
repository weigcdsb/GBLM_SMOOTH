addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%% set up simulation
clc; clear all; close all;
rng(3)

sim.seed = randi(10);
sim.T = 10*60;
sim.dt = 1e-3;
sim.vecN = round(sim.T/sim.dt);

sim.alpha_dt = 0.004;
sim.alpha_tau = 0.001;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1;
sim.hist_tau = .01;
sim.hist_beta = -1;
sim.stp_B = [0 1 2 3 4]'*-.03;

%% low pre - low post
simLL = sim;
simLL.pPreSpike = 5*simLL.dt;
simLL.wt_long = wtLongAdjust(1.5, simLL);
simLL.beta0 = beta0Adjust(15, simLL);
dataLL.dt = simLL.dt;
[dataLL,simLL] = sim_model(dataLL,simLL);

[fitLL,~] = smooth_gblm(dataLL.pre_spk_vec, dataLL.post_spk_vec,...
    'iter',30,...
    'hist_tau', simLL.hist_tau, 'hist_beta', simLL.hist_beta);

%% low pre - high post
simLH = sim;
simLH.pPreSpike = 5*simLH.dt;
simLH.wt_long = wtLongAdjust(1.5, simLH);
simLH.beta0 = beta0Adjust(30, simLH);
dataLH.dt = simLH.dt;
[dataLH,simLH] = sim_model(dataLH,simLH);

[fitLH,~] = smooth_gblm(dataLH.pre_spk_vec, dataLH.post_spk_vec,...
    'iter',30,...
    'hist_tau', simLH.hist_tau, 'hist_beta', simLH.hist_beta);


%% high pre - low post
simHL = sim;
simHL.pPreSpike = 10*simHL.dt;
simHL.wt_long = wtLongAdjust(1.5, simHL);
simHL.beta0 = beta0Adjust(15, simHL);
dataHL.dt = simHL.dt;
[dataHL,simHL] = sim_model(dataHL,simHL);

[fitHL,~] = smooth_gblm(dataHL.pre_spk_vec, dataHL.post_spk_vec,...
    'iter',30,...
    'hist_tau', simHL.hist_tau, 'hist_beta', simHL.hist_beta);


%% high pre - high post
simHH = sim;
simHH.pPreSpike = 10*simHH.dt;
simHH.wt_long = wtLongAdjust(1.5, simHH);
simHH.beta0 = beta0Adjust(30, simHH);
dataHH.dt = simHH.dt;
[dataHH,simHH] = sim_model(dataHH,simHH);

[fitHH,~] = smooth_gblm(dataHH.pre_spk_vec, dataHH.post_spk_vec,...
    'iter',30,...
    'hist_tau', simHH.hist_tau, 'hist_beta', simHH.hist_beta);

save('fireRate.mat');

%% plots
fireRatePlots(simLL, dataLL, fitLL);
fireRatePlots(simLH, dataLH, fitLH);
fireRatePlots(simHL, dataHL, fitHL);
fireRatePlots(simHH, dataHH, fitHH);


