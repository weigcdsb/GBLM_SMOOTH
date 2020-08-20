addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%% generate simulation data
clc;clear all;close all;
sim.seed = 123;
sim.dt = 0.001;
sim.pPreSpike = 5*sim.dt;
sim.T = 30*60;

% parameters for history basis
sim.postBaseRate = 20; %Hz
sim.mprops.nfilt = 5;
sim.mprops.delay = 400/(sim.dt*1000);
sim.b = [0 0.1 0.9 0.4 0.1]'/1.75; % coupling parameters
sim.hist_tau = 0.005; % history filter
sim.hist_beta = -10;

% stdp parameters
sim.stdp_params.noise = 0; % in s
sim.stdp_params.tau_forgetting = 15; % in s

sim.stdp_params.type = 'dexp';
sim.stdp_params.tau_plus = 20/1000; % in s
sim.stdp_params.tau_minus = 20/1000; % in s
sim.stdp_params.A_plus = 0.003; 
sim.stdp_params.A_minus = 0.5*sim.stdp_params.A_plus;

% sim.stdp_params.type = 'ahebb';
% sim.stdp_params.tau_plus = 20/1000; % in s
% sim.stdp_params.tau_minus = 50/1000; % in s
% sim.stdp_params.A_plus = 0.005;
% sim.stdp_params.A_minus = 0.475*sim.stdp_params.A_plus;

sim.stdp_params.g_max = 50;
sim.stdp_params.g_init = 1;

[data,sim] = sim_model_stdp(sim);

%% fit model
[fit,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',10, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

save('stdpSim.mat');

