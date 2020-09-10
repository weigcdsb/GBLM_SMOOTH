addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/helper'));
addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/core'));

%% set up true parameters
clc;clear all;close all;
T = 20*60;
dt = 0.001;
Q = diag([1e-5 1e-5]);

sim.seed = 123;
sim.T = T;
sim.dt = dt;
sim.vecN = round(sim.T/sim.dt);
% sim.pPreSpike = ones(sim.vecN,1)*3*sim.dt;
% sim.pPreSpike(round(sim.vecN/3):round(sim.vecN*2/3)) = 15*sim.dt;

sim.pPreSpike = ones(sim.vecN,1)*5*sim.dt;
sim.pPreSpike(round(sim.vecN/3):round(sim.vecN*2/3)) = 0*sim.dt;

sim.alpha_dt = 0.004;
sim.alpha_tau = 0.001;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1;
sim.stp_B = [0 0 1 3 1]'*.05;
sim.hist_tau = .005;
sim.hist_beta = -10;

sim.beta0 = ones(1, T/dt)'*4;
sim.wt_long = ones(1, T/dt)';
data.dt = sim.dt;
[data,sim] = sim_model(data,sim);

%%
[fit,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
    'iter',30, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);

save('uncertainty2.mat')



