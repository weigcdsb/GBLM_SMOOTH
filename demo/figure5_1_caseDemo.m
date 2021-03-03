addpath(genpath('~/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('~/GitHub/GBLM_SMOOTH/core'));

%% set up true parameters
clc;clear all;close all;
T = 20*60;
dt = 0.001;
Q = diag([1e-5 1e-5]);

fileName_seq =...
    {'1_supp_depression_linearJump.mat',...
    '2_supp_depression_sinJump.mat',...
    '3_supp_depression_linearLinear.mat',...
    '4_supp_depression_sinLinear.mat'};

period = T/(2*dt);

beta0_seq = ...
    [linspace(2, 4, T/dt)' ...
    3 + sin((2*pi/period)*(1:T/dt))' ...
    linspace(2, 4, T/dt)' ...
    3 + sin((2*pi/period)*(1:T/dt))'];


wt_long_seq = ...
    [[repmat(1.5, 1, round(T/(dt*2))) repmat(3.5, 1, T/dt - round(T/(dt*2)))]' ...
    [repmat(1.5, 1, round(T/(dt*2))) repmat(3.5, 1, T/dt - round(T/(dt*2)))]' ...
    linspace(1.5, 3.5, T/dt)' ...
    linspace(1.5, 3.5, T/dt)'];


%% run...
rng(123);
seed_seq = randperm(50, 4);

for k = 1:4
    fprintf('Case %02i...',k)
    
    % simulation set up
    sim.seed = seed_seq(k);
    sim.T = T;
    sim.dt = dt;
    sim.vecN = round(sim.T/sim.dt);
    sim.pPreSpike = 5*dt;
    sim.alpha_dt = 0.004;
    sim.alpha_tau = 0.001;
    sim.stp_Nq = 5;
    sim.stp_Nm = 450;
    sim.stp_Ns = 50;
    sim.stp_tau= 1;
    sim.stp_B = [0 0 1 3 1]'*(-0.1);
    
    sim.hist_tau = .01;
    sim.hist_beta = -2;
    
    sim.beta0 = beta0_seq(:, k);
    sim.wt_long = wt_long_seq(:, k);
    data.dt = sim.dt;
    [data,sim] = sim_model(data,sim);
    
    [fit,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
        'iter',10, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);
    
    save(fileName_seq{k});
    fprintf('\n')
end