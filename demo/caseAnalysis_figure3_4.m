addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%%
clc;clear all;close all;

% load('E:\suppResults_0803\1_supp_depression_linearJump.mat')
% load('E:\suppResults_0803\2_supp_depression_sinJump.mat')
% load('E:\suppResults_0803\3_supp_depression_linearLinear.mat')
% load('E:\suppResults_0803\4_supp_depression_linearSin.mat')
% load('E:\suppResults_0803\5_supp_depression_sinLinear.mat')
% load('E:\suppResults_0803\6_supp_depression_sinSin.mat')
% load('E:\suppResults_0803\7_supp_depression_conJump_5hzPre.mat')
% load('E:\suppResults_0803\8_supp_depression_conJump_10hzPre.mat')

% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig3_plot\1_supp_depression_linearJump';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig3_plot\2_supp_depression_sinJump';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig3_plot\3_supp_depression_linearLinear';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig3_plot\4_supp_depression_linearSin';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig3_plot\5_supp_depression_sinLinear';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig3_plot\6_supp_depression_sinSin';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig4_plot\7_supp_depression_conJump_5hzPre';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig4_plot\8_supp_depression_conJump_10hzPre';

cd(plotFolder)

%%
Tpre = data.pre_spk_times;
Tpost = data.post_spk_times(:, 1);
data.vecN = length(data.pre_spk_vec);

lam = exp(fit.beta0 + fit.wt_long.*fit.wt_short.*fit.Xc +...
    fit.hist*fit.hist_beta)*data.dt;

%% baseline
baseLine = figure;
idx = 1:size(fit.beta0);
hold on
plot(idx, fit.beta0, 'r', 'LineWidth', 3)
plot(idx, sim.beta0, 'k', 'LineWidth', 3)
plot(idx, fit.beta0 + sqrt(squeeze(fit.W(1, 1, :))), 'r:', 'LineWidth', 2)
plot(idx, fit.beta0 - sqrt(squeeze(fit.W(1, 1, :))), 'r:', 'LineWidth', 2)
ylim([min(sim.beta0)-1 max(sim.beta0)+1])
xlim([0 sim.T/sim.dt])
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(baseLine,'PaperUnits','inches','PaperPosition',[0 0 5 3])
saveas(baseLine, '1_baseLine.svg')
saveas(baseLine, '1_baseLine.png')

%% LTP
ltp = figure;
idx = 1:size(fit.beta0);
hold on
plot(idx, fit.wt_long, 'r', 'LineWidth', 3)
plot(idx, sim.wt_long, 'k', 'LineWidth', 3)
plot(idx, fit.wt_long + sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
plot(idx, fit.wt_long - sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
ylim([min(sim.wt_long)-1 max(sim.wt_long)+1])
xlim([0 sim.T/sim.dt])
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(ltp,'PaperUnits','inches','PaperPosition',[0 0 5 3])
saveas(ltp, '2_ltp.svg')
saveas(ltp, '2_ltp.png')

%% STP (modification function)

modFun = figure;
hold on
modPlot(sim, fit, 3, 2)
ylim([0.5 1.5])
if k >= 7; ylim([0.7 1.3]); end;
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(modFun,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(modFun, '3_modFun.svg')
saveas(modFun, '3_modFun.png')

%% firing rate (special for case 7 and 8)
if k >= 7
    firRate = figure;
    hold on
    plot(linspace(0,sim.T,sim.vecN), filter(ones(2000,1),1,data.pre_spk_vec)/2, 'b', 'LineWidth', 1.5);
    plot(linspace(0,sim.T,sim.vecN), filter(ones(2000,1),1,data.post_spk_vec)/2, 'k', 'LineWidth', 1.5);
    plot(linspace(0,sim.T,sim.vecN), filter(ones(2000,1),1,lam)/2, 'r', 'LineWidth', 1.5);
    hold off
    ylim([0 50])
    xlim([0 sim.T])
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    
    set(firRate,'PaperUnits','inches','PaperPosition',[0 0 5 3])
    saveas(firRate, '4_firRate.svg')
    saveas(firRate, '4_firRate.png')
end


%% split cross-correlogram (show LTP and baseline, seperate)
cd(strcat(plotFolder, '\corrT'))
maxd = -Inf;

for q=1:4
    qspk = Tpre(Tpre >= (q-1)*0.25*sim.T & Tpre < q*0.25*sim.T);
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    if maxd<max(d),maxd=max(d); end
end

for q=1:4
    qspk = Tpre(Tpre >= (q-1)*0.25*sim.T & Tpre < q*0.25*sim.T);
    qvec = zeros(1,data.vecN);
    qvec(round(qspk/data.dt))=1;
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    tvec = linspace(-0.025,0.025,64);
    tvec = tvec+mean(diff(tvec))/2;
    [dfit,lag_fit] = xcorr(qvec, lam, 25);
    
    corrT = figure;
    bar(tvec(1:end-1),d(1:end-1),1,'k','EdgeColor','none');
    hold on
    plot(-lag_fit*data.dt, dfit*mean(diff(tvec))/data.dt, 'r', 'LineWidth',3)
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    hold off
    ylim([0 ceil(maxd/50)*50])
    xlim([-0.01 0.02])
    
    set(corrT,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    saveas(corrT, strcat('corrT_Q', string(q), '.svg'))
    saveas(corrT, strcat('corrT_Q', string(q), '.png'))
end