addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%% set up true parameters
clc;clear all;close all;
T = 20*60;
dt = 0.001;
Q = diag([1e-5 1e-5]);

pPreSpike = [3 15]*dt;

fileName_seq =...
    {'1_supp_depression_conJump_3hzPre.mat',...
    '2_supp_depression_conJump_10hzPre.mat'};

beta0_seq = ...
    [ones(1, T/dt)'*5 ...
    ones(1, T/dt)'*5];

wt_long_seq = ...
    [[repmat(1, 1, round(T/(dt*2))) repmat(3, 1, T/dt - round(T/(dt*2)))]' ...
    [repmat(1, 1, round(T/(dt*2))) repmat(3, 1, T/dt - round(T/(dt*2)))]'];


%% run...
rng(123);
seed_seq = randperm(50, 2);

for k = 1:2
    fprintf('Case %02i...',k)
    
    % simulation set up
    sim.seed = seed_seq(k);
    sim.T = T;
    sim.dt = dt;
    sim.vecN = round(sim.T/sim.dt);
    sim.pPreSpike = pPreSpike(k);
    sim.alpha_dt = 0.004;
    sim.alpha_tau = 0.001;
    sim.stp_Nq = 5;
    sim.stp_Nm = 450;
    sim.stp_Ns = 50;
    sim.stp_tau= 1;
    sim.stp_B = [0 0 1 3 1]'*(-0.1);
    
    sim.hist_tau = .005;
    sim.hist_beta = -10;
    
    sim.beta0 = beta0_seq(:, k);
    sim.wt_long = wt_long_seq(:, k);
    data.dt = sim.dt;
    [data,sim] = sim_model(data,sim);
    
    [fit,~] = smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
        'iter',10, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);
    
    save(fileName_seq{k});
    fprintf('\n')
end