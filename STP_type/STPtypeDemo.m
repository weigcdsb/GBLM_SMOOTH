addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%% set up simulation
clc; clear all; close all;
rng(2)

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

%% facilitation
simFac = sim;
simFac.stp_B = [0 1 2 3 4]'*.04;
simFac.wt_long = wtLongAdjust(1.5, simFac);
simFac.beta0 = beta0Adjust(20, simFac);
dataFac.dt = simFac.dt;
[dataFac,simFac] = sim_model(dataFac,simFac);

[fitFac,~] = smooth_gblm(dataFac.pre_spk_vec, dataFac.post_spk_vec,...
    'iter',30,...
    'hist_tau', simFac.hist_tau, 'hist_beta', simFac.hist_beta);

%% depression
simDep = sim;
simDep.stp_B = [0 1 2 3 4]'*-.04;
simDep.wt_long = wtLongAdjust(1.5, simDep);
simDep.beta0 = beta0Adjust(20, simDep);
dataDep.dt = simDep.dt;
[dataDep,simDep] = sim_model(dataDep,simDep);


[fitDep,~] = smooth_gblm(dataDep.pre_spk_vec, dataDep.post_spk_vec,...
    'iter',30,...
    'hist_tau', simDep.hist_tau, 'hist_beta', simDep.hist_beta);

save('STPtype.mat');

%% plot

paramPlot(fitFac.beta0, squeeze(fitFac.W(1, 1, :)), simFac.beta0,...
    fitFac.wt_long, squeeze(fitFac.W(2, 2, :)),...
    simFac.wt_long, fitFac.wt_short, 1 + simFac.stp_X*simFac.stp_B,...
    fitFac.covB, fitFac.stp_X)

modPlot(simFac, fitFac)

paramPlot(fitDep.beta0, squeeze(fitDep.W(1, 1, :)), simDep.beta0,...
    fitDep.wt_long, squeeze(fitDep.W(2, 2, :)),...
    simDep.wt_long, fitDep.wt_short, 1 + simDep.stp_X*simDep.stp_B,...
    fitDep.covB, fitDep.stp_X)

modPlot(simDep, fitDep)


