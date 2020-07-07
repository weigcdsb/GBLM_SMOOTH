addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%% set up true parameters
clc;clear all;close all;
T = 20*60;
dt = 0.001;
Q = diag([1e-5 1e-5]);

trueParam_seq =...
    [[5 4 3 2 1]'*(0.04)...
    [5 4 3 2 1]'*(-0.04)...
    [0 0 0 0 0]'...
    [5 4 3 2 1]'*(-0.04)...
    [5 4 3 2 1]'*(0.04)];

fileName_seq =...
    {'1_facilitation_constant.mat',...
    '2_depression_constant.mat',...
    '3_noPlasticity_constant.mat',...
    '4_depression_linearJump.mat',...
    '5_facilitation_sinSin.mat'};

period = T/(2*dt);

beta0_seq = ...
    ... constant1 ...
    [ones(1, T/dt)'*3 ...
    ones(1, T/dt)'*3 ...
    ones(1, T/dt)'*4 ...
    ... constant2 ...
    linspace(4, 2, T/dt)'...
    3 + 1*sin(pi + (2*pi/period)*(1:T/dt))'];

wt_long_seq = ...
    ... constant ...
    [ones(1, T/dt)'*1.5 ...
    ones(1, T/dt)'*1.5 ...
    ones(1, T/dt)' ...
    ... change1 ...
    [repmat(0.5, 1, round(T/(dt*2))) repmat(2.5, 1, T/dt - round(T/(dt*2)))]'...
    1.5 + 1*sin((2*pi/period)*(1:T/dt))'];


%% run...
rng(123);
seed_seq = randperm(50, 5);

for k = 1:5
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
        'iter',5, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta);
    
    save(fileName_seq{k});
    fprintf('\n')
end










