addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%% set up simulation
clc; clear all; close all;
rng(2)
Q_true = diag([1e-5 1e-5]);


sim.seed = randi(10);
sim.T = 20*60;
sim.dt = 1e-3;
sim.vecN = round(sim.T/sim.dt);
sim.pPreSpike = 5*sim.dt;
sim.alpha_dt = 0.004;
sim.alpha_tau = 0.001;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1;
sim.stp_B = [0 0 1 3 1]'*(-0.05);
sim.hist_tau = .01;
sim.hist_beta = -2;
sim.beta0 = ones(1, sim.T/sim.dt)'*3 + ...
    detrend(cumsum(randn(sim.vecN,1)*sqrt(Q_true(1, 1))));
sim.wt_long = ones(1, sim.T/sim.dt)'*3 + ...
    detrend(cumsum(randn(sim.vecN,1)*sqrt(Q_true(2, 2))));

data.dt = sim.dt;
[data,sim] = sim_model(data,sim);

%% Q tune
[fit, fit_trace, Qvec, qbllhd_pred, qwllhd_pred, Qopt] = ...
    tune_smooth_gblm_1d_grid(data.pre_spk_vec, data.post_spk_vec,...
    'iter',15, 'llhdPlot', true, 'nq', 20,...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

%% no Q tune
% set Q too small
[fitSmall,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',15, 'Q', eye(2)*Qvec(4),...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

% sett Q too large
[fitLarge,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',15, 'Q', eye(2)*Qvec(18),...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

save('QTune.mat')
%% plots

subplot(1,2,1)
semilogx(Qvec,qbllhd_pred)
ylabel('log likelihood')
xlabel('Q')
title('beta0')
subplot(1,2,2)
semilogx(Qvec,qwllhd_pred)
xlabel('Q')
title('wtlong')


paramPlot(fitLarge.beta0, squeeze(fitLarge.W(1, 1, :)), sim.beta0,...
    fitLarge.wt_long, squeeze(fitLarge.W(2, 2, :)),...
    sim.wt_long, fitLarge.wt_short, 1 + sim.stp_X*sim.stp_B,...
    fitLarge.covB, fitLarge.stp_X)

paramPlot(fitSmall.beta0, squeeze(fitSmall.W(1, 1, :)), sim.beta0,...
    fitSmall.wt_long, squeeze(fitSmall.W(2, 2, :)),...
    sim.wt_long, fitSmall.wt_short, 1 + sim.stp_X*sim.stp_B,...
    fitSmall.covB, fitSmall.stp_X)

paramPlot(fit.beta0, squeeze(fit.W(1, 1, :)), sim.beta0,...
    fit.wt_long, squeeze(fit.W(2, 2, :)),...
    sim.wt_long, fit.wt_short, 1 + sim.stp_X*sim.stp_B,...
    fit.covB, fit.stp_X)


