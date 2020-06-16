addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%% facilitation
clc; clear all; close all;
load('F:\COVID-19\updateResults\convergence.mat')

rstart = 6;
sim = simDep;
fitL = fitL_dep;
fitS = fitS_dep;
fitL_trace = fitL_dep_trace;
fitS_trace = fitS_dep_trace;

%%
nTraceL = length(fitL_trace);
devTrace_fitL = ones(nTraceL, rstart)*NaN;

for k = 2:nTraceL-2
    devTrace_fitL(k, :) = fitL_trace(k).dev;
end

devTrace_fitL(devTrace_fitL == 0) = NaN;
plot(devTrace_fitL)

plot(fitL_trace(1).beta0)
plot(fitL_trace(3).beta0)
plot(fitL_trace(5).beta0)

plot(fitL_trace(1).wt_long)
plot(fitL_trace(3).wt_long)
plot(fitL_trace(5).wt_long)

plot(sim.stp_basis'*fitL_trace(1).wt_short_param)
plot(sim.stp_basis'*fitL_trace(3).wt_short_param)
plot(sim.stp_basis'*fitL_trace(5).wt_short_param)

hold on
plot(fitL.beta0, 'r')
plot(fitS.beta0, 'b')
plot(sim.beta0, 'k')
hold off

hold on
plot(fitL.wt_long, 'r')
plot(fitS.wt_long, 'b')
plot(sim.wt_long, 'k')
hold off

hold on
plot(sim.stp_basis'*fitL.wt_short_param, 'r')
plot(sim.stp_basis'*fitS.wt_short_param, 'b')
plot(sim.stp_basis'*sim.stp_B, 'k')
hold off



%% depression













%%
hold on
plot(fitL.beta0, 'r')
plot(fitS.beta0, 'b')
hold off

hold on
plot(fitL.wt_long, 'r')
plot(fitS.wt_long, 'b')
hold off

hold on
plot(fitL.wt_short, 'r')
plot(fitS.wt_short, 'b')
hold off

hold on
plot(fitL_trace(5).beta0, 'r')
plot(fitS_trace(5).beta0, 'b')
hold off


