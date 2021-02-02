addpath(genpath('~/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('~/GitHub/GBLM_SMOOTH/core'));

%% set up parameters
clc;clear all;close all;
rng(8)
T = 20*60;
dt = 0.001;
Q = diag([1e-5 1e-5]);
trueParam = [0 0 1 3 1]'*(-0.05);

n = 30;
nb0 = normrnd(3, 1, n, 1);
beta0 = (interp1(linspace(0,T,length(nb0)),nb0,linspace(0,T,T/dt),'spline'))';
% plot(beta0)

wt_long = ones(1, T/dt)'*3;
sim.seed = randperm(50, 1);
sim.T = T;
sim.dt = dt;
sim.vecN = round(sim.T/sim.dt);

np0 = normrnd(8*sim.dt, 4*sim.dt, n, 1);
sim.pPreSpike = interp1(linspace(0,T,length(np0)),np0,linspace(0,T,T/dt),'spline');
sim.pPreSpike(sim.pPreSpike < 0) = 0;
plot(sim.pPreSpike)


sim.alpha_dt = 0.004;
sim.alpha_tau = 0.002;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1;
sim.stp_B = trueParam;
sim.hist_tau = .01;
sim.hist_beta = -2;

sim.beta0 = beta0;
sim.wt_long = wt_long;
data.dt = sim.dt;
[data,sim] = sim_model(data,sim);

[fit,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',10, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta, 'Q', Q);

save('complete2.mat')

disp(fit.synParams.syn_params(1))
disp(fit.synParams.syn_params(2))
%% beta0 constant
Q = diag([0 1e-5]);

synParams = fit.synParams;
clear fit;
[fit,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',20, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta,...
    'Q', Q, 'synParams', synParams);

save('b0Constant.mat')

%% no estimation on wt_short
data.vecN = length(data.pre_spk_vec);
clear fit;
fit.hist_tau = sim.hist_tau;
fit.hist_beta = sim.hist_beta;

% parameters for GBLM for STP
fit.stp_Nq = 5;
fit.stp_Nm = 450;
fit.stp_Ns = 50;
fit.stp_tau= 1;

% parameters for synaptic connection
fit.syn_hyper_params.coupling_timescale = 0.05;
fit.syn_hyper_params.bin_width = (0.05)*data.dt;
fit.syn_hyper_params.baseline_nsplines = 4;

% parameters for smoothing
fit.iter = 10;
fit.toleranceValue= 1e-6;
fit.F = eye(2);
fit.Q = diag([1e-5 1e-5]);
fit.doFiltOnly = false;
fit.doOracle = false;

% Synaptic Connection
fit.synParams = synParams;
fit = synConEst(data,fit);

fit.W = zeros(2, 2, data.vecN);
offset = log(data.dt) + fit.hist*fit.hist_beta;
alph = glmfit([fit.Xc],data.post_spk_vec,'poisson','Offset',offset);
fit.wt_long = ones(data.vecN, 1)*alph(2);
fit.beta0 = ones(data.vecN,1)*alph(1);

W0 = eye(2)*0.1;
b0 = [fit.beta0(1, :) fit.wt_long(1, :)]';
X = [ones(data.vecN,1) fit.Xc];
offset = fit.hist*fit.hist_beta;
[b, W, ~, ~] = ppasmoo_poissexp(data.post_spk_vec, X, b0, W0, fit.F, fit.Q, offset);

fit.beta0 = b(1,:)';
fit.wt_long = b(2,:)';
fit.W = W;


save('noWt_short.mat')


% diagnose
idx = 1:size(fit.beta0);
subplot(2, 1, 1)
plot(idx, fit.beta0,...
    idx, fit.beta0 + sqrt(squeeze(fit.W(1, 1, :))), 'b:',...
    idx, fit.beta0 - sqrt(squeeze(fit.W(1, 1, :))), 'b:',...
    idx, sim.beta0, 'r');
title('beta_0')

subplot(2, 1, 2)
hold on
plot(idx, sim.wt_long, 'r',...
    idx, fit.wt_long, 'b',...
    idx, fit.wt_long + sqrt(squeeze(fit.W(2, 2, :))), 'b:',...
    idx, fit.wt_long - sqrt(squeeze(fit.W(2, 2, :))), 'b:');
plot(idx, sim.wt_long*mean(1 + sim.stp_X*sim.stp_B))

title('wt_{long}')
hold off

