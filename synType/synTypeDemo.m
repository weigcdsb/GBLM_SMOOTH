addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%% set up simulation
clc; clear all; close all;

sim.seed = 3;
sim.T = 60*60;
sim.dt = 1e-3;
sim.vecN = round(sim.T/sim.dt);
sim.pPreSpike = 5*sim.dt;
sim.alpha_dt = 0.004;
sim.alpha_tau = 0.001;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1;
sim.stp_B = [0 1 2 3 4]'*.04;
sim.hist_tau = .01;
sim.hist_beta = -3;

% ecitatory
simEx = sim;
simEx.wt_long = ones(simEx.vecN, 1)*1.5;
simEx.beta0 = beta0Adjust(10, simEx);
dataEx.dt = simEx.dt;
[dataEx,simEx] = sim_model(dataEx,simEx);

% inhibitory
simIn = sim;
simIn.wt_long = -ones(simIn.vecN, 1)*2;
simIn.beta0 = beta0Adjust(20, simIn);
dataIn.dt = simIn.dt;
[dataIn,simIn] = sim_model(dataIn,simIn);


%% fit the model

[fitEx,~] = smooth_gblm(dataEx.pre_spk_vec, dataEx.post_spk_vec,...
    'iter',30,...
    'hist_tau', simEx.hist_tau, 'hist_beta', simEx.hist_beta);

[fitIn,~] = smooth_gblm(dataIn.pre_spk_vec, dataIn.post_spk_vec,...
    'iter',30,...
    'hist_tau', simIn.hist_tau, 'hist_beta', simIn.hist_beta);

save('synType.mat')

%% show the results
paramPlot(fitEx.beta0, squeeze(fitEx.W(1, 1, :)), simEx.beta0,...
    fitEx.wt_long, squeeze(fitEx.W(2, 2, :)),...
    simEx.wt_long, fitEx.wt_short, 1 + simEx.stp_X*simEx.stp_B,...
    fitEx.covB, fitEx.stp_X)

modPlot(simEx, fitEx)


paramPlot(fitIn.beta0, squeeze(fitIn.W(1, 1, :)), simIn.beta0,...
    fitIn.wt_long, squeeze(fitIn.W(2, 2, :)),...
    simIn.wt_long, fitIn.wt_short, 1 + simIn.stp_X*simIn.stp_B,...
    fitIn.covB, fitIn.stp_X)

modPlot(simIn, fitIn)






















