addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/helper'));
addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/core'));

%%
clc;clear all;close all;

% load('E:\casedDemoResults_supp_0801\1_facilitation_conJump.mat');
% load('E:\casedDemoResults_supp_0801\2_depression_conJump.mat');
% load('E:\casedDemoResults_supp_0801\3_noPlasticity_conJump.mat');
% load('E:\casedDemoResults_supp_0801\4_depression_conLinear.mat');
% load('E:\casedDemoResults_supp_0801\5_depression_conSin.mat');
% load('E:\caseDemoResults_0715\4_depression_linearJump.mat');
load('E:\caseDemoResults_0715\5_depression_sinSin.mat');

% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case_supp2\1_facilitation_conJump';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case_supp2\2_depression_conJump';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case_supp2\3_noPlasticity_conJump';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case_supp2\4_depression_conLinear';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case_supp2\5_depression_conSin';
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case_supp2\6_depression_linearJump';
plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case_supp2\7_depression_sinSin';

cd(plotFolder)
%% baseline
baseLine = figure;
idx = 1:size(fit.beta0);
hold on
plot(idx, fit.beta0, 'r', 'LineWidth', 3)
plot(idx, sim.beta0, 'k', 'LineWidth', 3)
plot(idx, fit.beta0 + sqrt(squeeze(fit.W(1, 1, :))), 'r:', 'LineWidth', 2)
plot(idx, fit.beta0 - sqrt(squeeze(fit.W(1, 1, :))), 'r:', 'LineWidth', 2)
ylim([min(sim.beta0)-1 max(sim.beta0)+1])
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off

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
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off

saveas(ltp, '2_ltp.svg')
saveas(ltp, '2_ltp.png')

%% STP (modification function)
modFun = figure;
hold on
modPlot(sim, fit, 3, 2)
ylim([0.5 1.5])
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off

saveas(modFun, '3_modFun.svg')
saveas(modFun, '3_modFun.png')


