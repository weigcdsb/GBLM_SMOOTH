addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/helper'));
addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/core'));

%% set up true parameters
clc;clear all;close all;
T = 20*60;
dt = 0.001;
Q = diag([1e-5 1e-5]);

trueParam_seq =...
    [[0 0 1 3 1]'*(0.1)...
    [0 0 1 3 1]'*(-0.1)...
    [0 0 0 0 0]'...
    [0 0 1 3 1]'*(-0.1)...
    [0 0 1 3 1]'*(-0.1)];

fileName_seq =...
    {'1_facilitation_conJump.mat',...
    '2_depression_conJump.mat',...
    '3_noPlasticity_conJump.mat',...
    '4_depression_conLinear.mat',...
    '5_depression_conSin.mat'};

period = T/(2*dt);

beta0_seq = ...
    ... constant1 ...
    [ones(1, T/dt)'*2.5 ...
    ones(1, T/dt)'*3 ...
    ones(1, T/dt)'*3 ...
    ones(1, T/dt)'*3 ...
    ones(1, T/dt)'*3];


% plot([repmat(2, 1, round(T/(dt*2))) repmat(4, 1, T/dt - round(T/(dt*2)))]')

wt_long_seq = ...
    [[repmat(0.5, 1, round(T/(dt*2))) repmat(1.5, 1, T/dt - round(T/(dt*2)))]' ...
    [repmat(1.5, 1, round(T/(dt*2))) repmat(3.5, 1, T/dt - round(T/(dt*2)))]' ...
    [repmat(1, 1, round(T/(dt*2))) repmat(3, 1, T/dt - round(T/(dt*2)))]' ...
    linspace(1.5, 3.5, T/dt)'...
    2.5 + sin((2*pi/period)*(1:T/dt))'];

%% run...
rng(123);
seed_seq = randperm(50, 5);

for k = 3
    fprintf('Case %02i...',k)
    
    % simulation set up
    sim.seed = seed_seq(k);
    sim.T = T;
    sim.dt = dt;
    sim.vecN = round(sim.T/sim.dt);
    sim.pPreSpike = 5*sim.dt;
    sim.alpha_dt = 0.004;
    sim.alpha_tau = 0.001;
    sim.stp_Nq = 5;
    sim.stp_Nm = 450;
    sim.stp_Ns = 50;
    sim.stp_tau= 1;
    sim.stp_B = trueParam_seq(:, k);
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