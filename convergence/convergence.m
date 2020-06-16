addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%% set up simulation
clc; clear all; close all;
simParam.T = 10*60;
simParam.dt = 1e-3;
simParam.vecN = round(simParam.T/simParam.dt);
simParam.pPreSpike = 5*simParam.dt;
simParam.alpha_dt = 0.004;
simParam.alpha_tau = 0.001;
simParam.stp_Nq = 5;
simParam.stp_Nm = 450;
simParam.stp_Ns = 50;
simParam.stp_tau= 1;
simParam.hist_tau = .01;
simParam.hist_beta = -1;
data.dt = simParam.dt;

%% facilitation
dataFac = data;
simFac = simParam;

simFac.seed = 1;
simFac.stp_B = [0 1 2 3 4]'*.04;
simFac.wt_long = wtLongAdjust(1.5, simFac);
simFac.beta0 = beta0Adjust(20, simFac); % mean lambda is approximately 10Hz
[dataFac,simFac] = sim_model(dataFac,simFac);


[fitL_fac, fitL_fac_trace] = smooth_gblm_rand(dataFac.pre_spk_vec, dataFac.post_spk_vec,...
    'iter', 30, 'longTermFirst', true, 'hist_tau', simFac.hist_tau, 'hist_beta', simFac.hist_beta);
[fitS_fac, fitS_fac_trace] = smooth_gblm_rand(dataFac.pre_spk_vec, dataFac.post_spk_vec,...
    'iter', 30, 'longTermFirst', false, 'hist_tau', simFac.hist_tau, 'hist_beta', simFac.hist_beta);

%% depression
dataDep = data;
simDep = simParam;

simDep.seed = 2;
simDep.stp_B = [0 1 2 3 4]'*(-.04);
simDep.wt_long = wtLongAdjust(1.5, simDep);
simDep.beta0 = beta0Adjust(20, simDep); % mean lambda is approximately 10Hz
[dataDep,simDep] = sim_model(dataDep,simDep);


[fitL_dep, fitL_dep_trace] = smooth_gblm_rand(dataDep.pre_spk_vec, dataDep.post_spk_vec,...
    'iter', 30, 'longTermFirst', true, 'hist_tau', simDep.hist_tau, 'hist_beta', simDep.hist_beta);
[fitS_dep, fitS_dep_trace] = smooth_gblm_rand(dataDep.pre_spk_vec, dataDep.post_spk_vec,...
    'iter', 30, 'longTermFirst', false, 'hist_tau', simDep.hist_tau, 'hist_beta', simDep.hist_beta);

save('convergence.mat')




