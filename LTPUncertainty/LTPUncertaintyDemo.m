addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%% set up simulation
clc; clear all; close all;
rng(123)

sim.seed = randi(10);
sim.T = 20*60;
sim.dt = 1e-3;
sim.vecN = round(sim.T/sim.dt);
sim.pPreSpike = ones(sim.vecN,1)*5*sim.dt;
sim.pPreSpike(round(sim.vecN/3):round(sim.vecN*2/3)) = 0;
sim.alpha_dt = 0.004;
sim.alpha_tau = 0.001;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1;
sim.hist_tau = .01;
sim.hist_beta = -1;


%% depression
simDep = sim;
simDep.stp_B = [0 1 2 3 4]'*-.04;
simDep.wt_long = wtLongAdjust(1.5, simDep);
simDep.beta0 = beta0Adjust(20, simDep);
dataDep.dt = simDep.dt;
[dataDep,simDep] = sim_model(dataDep,simDep);

[fitDepAll, ~] = smooth_gblm(dataDep.pre_spk_vec, dataDep.post_spk_vec,...
    'iter',30,...
    'hist_tau', simDep.hist_tau, 'hist_beta', simDep.hist_beta);

[fitDep,~] = smooth_gblm(dataDep.pre_spk_vec, dataDep.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_stp_B', simDep.stp_B,...
    'hist_tau', simDep.hist_tau, 'hist_beta', simDep.hist_beta);

[fitDepMore,~] = smooth_gblm(dataDep.pre_spk_vec, dataDep.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_stp_B', simDep.stp_B,...
    'oracle_beta0', simDep.beta0,...
    'hist_tau', simDep.hist_tau, 'hist_beta', simDep.hist_beta);


%% facilitation
simFac = sim;
simFac.stp_B = [0 1 2 3 4]'*.04;
simFac.wt_long = wtLongAdjust(1.5, simFac);
simFac.beta0 = beta0Adjust(20, simFac);
dataFac.dt = simFac.dt;
[dataFac,simFac] = sim_model(dataFac,simFac);

[fitFacAll,~] = smooth_gblm(dataFac.pre_spk_vec, dataFac.post_spk_vec,...
    'iter',30,...
    'hist_tau', simFac.hist_tau, 'hist_beta', simFac.hist_beta);

[fitFac,~] = smooth_gblm(dataFac.pre_spk_vec, dataFac.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_stp_B', simFac.stp_B,...
    'hist_tau', simFac.hist_tau, 'hist_beta', simFac.hist_beta);

[fitFacMore,~] = smooth_gblm(dataFac.pre_spk_vec, dataFac.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_stp_B', simFac.stp_B,...
    'oracle_beta0', simFac.beta0,...
    'hist_tau', simFac.hist_tau, 'hist_beta', simFac.hist_beta);

%% no plasticity
simNull = sim;
simNull.stp_B = zeros(5, 1);
simNull.wt_long = wtLongAdjust(1.5, simNull);
simNull.beta0 = beta0Adjust(20, simNull);
dataNull.dt = simNull.dt;
[dataNull,simNull] = sim_model(dataNull,simNull);

[fitNullAll,~] = smooth_gblm(dataNull.pre_spk_vec, dataNull.post_spk_vec,...
    'iter',30,...
    'hist_tau', simNull.hist_tau, 'hist_beta', simNull.hist_beta);

[fitNull,~] = smooth_gblm(dataNull.pre_spk_vec, dataNull.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_stp_B', simNull.stp_B,...
    'hist_tau', simNull.hist_tau, 'hist_beta', simNull.hist_beta);

[fitNullMore,~] = smooth_gblm(dataNull.pre_spk_vec, dataNull.post_spk_vec,...
    'iter',30, 'doOracle', true, 'oracle_stp_B', simNull.stp_B,...
    'oracle_beta0', simNull.beta0,...
    'hist_tau', simNull.hist_tau, 'hist_beta', simNull.hist_beta);

save('uncertainty.mat')

%% plot
LTPUPlot(simDep, dataDep, fitDepAll)
LTPUPlot(simDep, dataDep, fitDep)
LTPUPlot(simDep, dataDep, fitDepMore)

LTPUPlot(simFac, dataFac, fitFacAll)
LTPUPlot(simFac, dataFac, fitFac)
LTPUPlot(simFac, dataFac, fitFacMore)

LTPUPlot(simNull, dataNull, fitNullAll)
LTPUPlot(simNull, dataNull, fitNull)
LTPUPlot(simNull, dataNull, fitNullMore)
