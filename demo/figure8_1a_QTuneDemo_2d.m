addpath(genpath('~/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('~/GitHub/GBLM_SMOOTH/core'));

%%
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
tic
[fit, ~, Qvec, llhdmesh, Qopt_grid] = ...
    tune_smooth_gblm_2d_grid(data.pre_spk_vec, data.post_spk_vec,...
    'iter',15, 'nq', 20,'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta,...
    'doFit', false);
t_grid = toc;

tic
[~, ~, Qopt_grad] = ...
    tune_smooth_gblm_2d_grad(data.pre_spk_vec, data.post_spk_vec,...
    'iter',15, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta,...
    'doFit', false, 'synParams', fit.synParams);
t_grad = toc;


data.vecN = length(data.pre_spk_vec);
nQ = length(Qvec);
llhd_mle_QLTP = zeros(1, nQ);
llhd_mle_Qbase = zeros(1, nQ);

for j = 1:nQ
    
    fit.Q = diag([Qopt_grad(1) Qvec(j)]);
    [~, fit, ~] = evalc('loopCore(data, fit)');
    llhd_mle_QLTP(j) = fit.llhd_pred;
    
    fit.Q = diag([Qvec(j) Qopt_grad(2)]);
    [~, fit, ~] = evalc('loopCore(data, fit)');
    llhd_mle_Qbase(j) = fit.llhd_pred;

end

save('QoptDemo.mat')

%%

% opt_grad
[fitMLE,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',15, 'Q', diag(Qopt_grad),...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta, 'synParams', fit.synParams);

% set Q too small
[fitSmall,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',15, 'Q', eye(2)*Qvec(4),...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta, 'synParams', fit.synParams);

% sett Q too large
[fitLarge,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',15, 'Q', eye(2)*Qvec(18),...
    'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta, 'synParams', fit.synParams);

save('QoptDemo2.mat');



















