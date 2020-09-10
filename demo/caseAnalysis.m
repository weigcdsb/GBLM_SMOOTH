addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%%
clc;clear all;close all;

%% pre fit
figure(1)
subplot(1,3,1:2)
plot(sim.stp_X*sim.stp_B+1)
subplot(1,3,3)
d = corr_fast_v3(data.pre_spk_times,data.post_spk_times(:,1),-.025,.025,64);
t = linspace(-.025,.025,64);
t = t+mean(diff(t))/2;
bar(t,d,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')

hold on
plot(filter(ones(2000,1),1,data.pre_spk_vec)/2, 'b');
plot(filter(ones(2000,1),1,data.post_spk_vec)/2, 'k');
hold off

plot(1 + sim.stp_basis'*sim.stp_B)

plot(sim.lam*dt)

%% post fit
figure(2)
paramPlot(fit.beta0, squeeze(fit.W(1, 1, :)), sim.beta0,...
    fit.wt_long, squeeze(fit.W(2, 2, :)),...
    sim.wt_long, fit.wt_short, 1 + sim.stp_X*sim.stp_B,...
    fit.covB, fit.stp_X)

modPlot(sim, fit, 3, 2)
