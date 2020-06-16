addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%%
clc;clear all;close all;
load('F:\COVID-19\updateResults\caseResults\2_depression_constant.mat');

%% diagnose (informal)
paramPlot(fit.beta0, squeeze(fit.W(1, 1, :)), sim.beta0,...
    fit.wt_long, squeeze(fit.W(2, 2, :)),...
    sim.wt_long, fit.wt_short, 1 + sim.stp_X*sim.stp_B,...
    fit.covB, fit.stp_X)

modPlot(sim, fit)
%% give the fitted post spike
data_fit.dt = data.dt;
fit.seed = sim.seed;
fit.pPreSpike = sim.pPreSpike;
fit.T = sim.T;
fit.vecN = sim.vecN;
fit.alpha_dt = fit.synParams.syn_params(1);
fit.alpha_tau = fit.synParams.syn_params(2);
fit.stp_B = fit.wt_short_param;
[data_fit,~] = sim_model(data_fit,fit);

Tpre = data.pre_spk_times;
Tpost = data.post_spk_times(:, 1);
Tpre_fit = data_fit.pre_spk_times;
Tpost_fit = data_fit.post_spk_times(:, 1);

%% pre and post firing rate
hold on
lpre = plot(filter(ones(2000,1),1,data.pre_spk_vec)/2, 'b');
lpst = plot(filter(ones(2000,1),1,data.post_spk_vec)/2, 'k');
lpst_fit = plot(filter(ones(2000,1),1,data_fit.post_spk_vec)/2, 'r');
lpre.Color(4) = 0.8;
lpst.Color(4) = 0.8;
lpst_fit.Color(4) = 0.2;
hold off
legend('pre', 'post', 'post-fit')
ylabel('Firing Rates')

%% Overall Cross Correlogram
[d,~] = corr_fast_v3(Tpre, Tpost,-.01,.1,110);
[d_fit, ~] = corr_fast_v3(Tpre_fit, Tpost_fit,-.01,.1,110);

hold on
bar(linspace(-.01,.1,110),d,'k');
plot(linspace(-.01,.1,110),d_fit, 'r', 'LineWidth',2);
xlim([-.01 .1]);
title('Overall Cross Correlogram')
xlabel('s')
ylabel('Count')
legend('data', 'fitted results')
hold off

%% show STP
isi = diff(find(data.pre_spk_vec>0)*data.dt);
isiQ1 = prctile(isi,25);
isiQ2 = prctile(isi,50); 
isiQ3 = prctile(isi,75);

TpreISI1 = Tpre(find(isi<isiQ1)+1);
TpreISI2 = Tpre(find(isi >= isiQ1 & isi<isiQ2)+1);
TpreISI3 = Tpre(find(isi >= isiQ2 & isi<isiQ3)+1);
TpreISI4 = Tpre(find(isi >= isiQ3)+1);

[dISI1,~] = corr_fast_v3(TpreISI1, Tpost,-.01,.1,110);
[dISI2,~] = corr_fast_v3(TpreISI2, Tpost,-.01,.1,110);
[dISI3,~] = corr_fast_v3(TpreISI3, Tpost,-.01,.1,110);
[dISI4,~] = corr_fast_v3(TpreISI4, Tpost,-.01,.1,110);

[dISI1_fit,~] = corr_fast_v3(TpreISI1, Tpost_fit,-.01,.1,110);
[dISI2_fit,~] = corr_fast_v3(TpreISI2, Tpost_fit,-.01,.1,110);
[dISI3_fit,~] = corr_fast_v3(TpreISI3, Tpost_fit,-.01,.1,110);
[dISI4_fit,~] = corr_fast_v3(TpreISI4, Tpost_fit,-.01,.1,110);


yULim = ceil(max([dISI1; dISI2; dISI3; dISI4])/50)*50;
% subplot(4,1,1)
subplot(2,2,1)
hold on
bar(linspace(-.01,.1,110),dISI1,'k');
plot(linspace(-.01,.1,110),dISI1_fit,'r', 'LineWidth',2)
title('Correlogram for < Q_1 of ISI')
xlim([-.01,.1]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off

% subplot(4,1,2)
subplot(2,2,2)
hold on
bar(linspace(-.01,.1,110),dISI2,'k');
plot(linspace(-.01,.1,110),dISI2_fit,'r', 'LineWidth',2)
title('Correlogram for Q_1 to Q_2 of ISI')
xlim([-.01,.1]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off

% subplot(4,1,3)
subplot(2,2,3)
hold on
bar(linspace(-.01,.1,110),dISI3,'k');
plot(linspace(-.01,.1,110),dISI3_fit,'r', 'LineWidth',2)
title('Correlogram for Q_2 to Q_3 of ISI')
xlim([-.01,.1]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off

% subplot(4,1,4)
subplot(2,2,4)
hold on
bar(linspace(-.01,.1,110),dISI4,'k');
plot(linspace(-.01,.1,110),dISI4_fit,'r', 'LineWidth',2)
title('Correlogram for > Q_3 of ISI')
xlim([-.01,.1]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off


%% show LTP

Tpre1 = Tpre(Tpre<(0.25*sim.T));
Tpre2 = Tpre(Tpre >= (0.25*sim.T) & Tpre<(0.5*sim.T));
Tpre3 = Tpre(Tpre >= (0.5*sim.T) & Tpre<(0.75*sim.T));
Tpre4 = Tpre(Tpre >= (0.75*sim.T));

[d1,~] = corr_fast_v3(Tpre1, Tpost,-.01,.1,110);
[d2,~] = corr_fast_v3(Tpre2, Tpost,-.01,.1,110);
[d3,~] = corr_fast_v3(Tpre3, Tpost,-.01,.1,110);
[d4,~] = corr_fast_v3(Tpre4, Tpost,-.01,.1,110);

[d1_fit,~] = corr_fast_v3(Tpre1, Tpost_fit,-.01,.1,110);
[d2_fit,~] = corr_fast_v3(Tpre2, Tpost_fit,-.01,.1,110);
[d3_fit,~] = corr_fast_v3(Tpre3, Tpost_fit,-.01,.1,110);
[d4_fit,~] = corr_fast_v3(Tpre4, Tpost_fit,-.01,.1,110);

yULim = ceil(max([d1; d2; d3; d4])/50)*50;
% subplot(4,1,1)
subplot(2,2,1)
hold on
bar(linspace(-.01,.1,110),d1,'k');
plot(linspace(-.01,.1,110),d1_fit,'r', 'LineWidth',2)
title('Correlogram for <Q_1 of T')
xlim([-.01,.1]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off

% subplot(4,1,2)
subplot(2,2,2)
hold on
bar(linspace(-.01,.1,110),d2,'k');
plot(linspace(-.01,.1,110),d2_fit,'r', 'LineWidth',2)
title('Correlogram for Q_1 to Q_2 of T')
xlim([-.01,.1]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off


% subplot(4,1,3)
subplot(2,2,3)
hold on
bar(linspace(-.01,.1,110),d3,'k');
plot(linspace(-.01,.1,110),d3_fit,'r', 'LineWidth',2)
title('Correlogram for Q_2 to Q_3 of T')
xlim([-.01,.1]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off


% subplot(4,1,4)
subplot(2,2,4)
hold on
bar(linspace(-.01,.1,110),d4,'k');
plot(linspace(-.01,.1,110),d4_fit,'r', 'LineWidth',2)
title('Correlogram for > Q_3 of T')
xlim([-.01,.1]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off
