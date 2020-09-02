addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%% generate simulation data
clc;clear all;close all;
sim.seed = 123;
sim.dt = 0.001;
sim.pPreSpike = 5*sim.dt;
sim.T = 20*60;

% parameters for history basis
sim.postBaseRate = 15; %Hz
sim.hist_tau = 0.0005; % history filter
sim.hist_beta = -10;

sim.alpha_dt = 0.004;
sim.alpha_tau = 0.001;

% parameters for (short-term) plasticity
sim.wt_long = 1;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1;
sim.stp_B = [0 0 1 3 1]'*(-0.05);

% stdp parameters
sim.stdp_params.noise = 0; % in s
sim.stdp_params.tau_forgetting = 20; % in s

sim.stdp_params.type = 'dexp';
sim.stdp_params.tau_plus = 20/1000; % in s
sim.stdp_params.tau_minus = 20/1000; % in s
sim.stdp_params.A_plus = 0.006;
sim.stdp_params.A_minus = 0.002;

sim.stdp_params.g_max = 50;
sim.stdp_params.g_init = 1;

[data,sim] = sim_model_stdp(sim);

%% fit model
[fit,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',10, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

save('stdpSim.mat');
%%

figure(1)
subplot(1,3,1:2)
plot(sim.g)
subplot(1,3,3)
d = corr_fast_v3(data.pre_spk_times,data.post_spk_times(:,1),-.025,.025,64);
t = linspace(-.025,.025,64);
t = t+mean(diff(t))/2;
bar(t,d,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')
