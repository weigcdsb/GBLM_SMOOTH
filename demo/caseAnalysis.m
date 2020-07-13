addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%%
clc;clear all;close all;
load('F:\COVID-19\updateResults\caseResults\2_depression_constant.mat');

% cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\1_facilitation_constant'
% cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\2_depression_constant'
% cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\3_no_plasticity_constant'
% cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\4_depression_linear_jump'
cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\5_depression_sin_sin'

%% diagnose (informal)
paramPlot(fit.beta0, squeeze(fit.W(1, 1, :)), sim.beta0,...
    fit.wt_long, squeeze(fit.W(2, 2, :)),...
    sim.wt_long, fit.wt_short, 1 + sim.stp_X*sim.stp_B,...
    fit.covB, fit.stp_X)

modFun = figure;
modPlot(sim, fit)

saveas(modFun, '1_modificationFuntion.svg')
saveas(modFun, '1_modificationFuntion.png')
%% calculate lambda

Tpre = data.pre_spk_times;
Tpost = data.post_spk_times(:, 1);

lam = exp(fit.beta0 + fit.wt_long.*fit.wt_short.*fit.Xc +...
    fit.hist*fit.hist_beta)*data.dt;

%% pre and post firing rate
firRate = figure;
hold on
plot(filter(ones(2000,1),1,data.pre_spk_vec)/2, 'b');
plot(filter(ones(2000,1),1,data.post_spk_vec)/2, 'k');
plot(filter(ones(2000,1),1,lam)/2, 'r');
hold off
legend('pre', 'post', 'post-fit')
ylabel('Firing Rates')

saveas(firRate, '2_firRate.svg')
saveas(firRate, '2_firRate.png')

%% Overall Cross Correlogram
[d,~] = corr_fast_v3(Tpre, Tpost,-.02,.02,102);
[d_fit, lag_fit] = xcorr(data.pre_spk_vec, lam, 20);

tvec = linspace(-0.02,0.02,102);
tvec = tvec+mean(diff(tvec))/2;

corrOverall = figure;
hold on
bar(tvec(1:end-1),d(1:end-1),1,'k','EdgeColor','none');
plot(-lag_fit*data.dt, d_fit*mean(diff(tvec))/data.dt, 'r', 'LineWidth',2)
xlim([-.01 .02]);
title('Overall Cross Correlogram')
xlabel('s')
ylabel('Count')
legend('data', 'fitted results')
hold off

saveas(corrOverall, '3_overallCorr.svg')
saveas(corrOverall, '3_overallCorr.png')

%% show STP

isi = [Inf; diff(data.pre_spk_times)];
quantiles = prctile(isi,linspace(0,100,5));
maxd = -Inf;
for q=1:length(quantiles)-1
    qspk = data.pre_spk_times(isi>=quantiles(q) & isi<quantiles(q+1));
    qvec = zeros(1,data.vecN);
    qvec(round(qspk/data.dt))=1;
    d = corr_fast_v3(qspk,data.post_spk_times,-.025,.025,64);
    tvec = linspace(-0.025,0.025,64);
    tvec = tvec+mean(diff(tvec))/2;
    [dfit,lag_fit] = xcorr(qvec, lam, 25);
    
    subplot(1,length(quantiles)-1,q)
    bar(tvec(1:end-1),d(1:end-1),1,'k','EdgeColor','none');
    hold on
    plot(-lag_fit*data.dt, dfit*mean(diff(tvec))/data.dt, 'r', 'LineWidth',2)
    title(sprintf('%0.2f ms',quantiles(q)*1000))
    xlabel('ms')
    ylabel('Count')
    hold off
    
    if maxd<max(d),maxd=max(d); end
end

for q=1:length(quantiles)-1
    subplot(1,length(quantiles)-1,q)
    ylim([0 maxd])
end

% saveas(corrISI, '4_corrISI.svg')
% saveas(corrISI, '4_corrISI.png')

%% show LTP

Tpre1 = Tpre(Tpre<(0.25*sim.T));
Tpre2 = Tpre(Tpre >= (0.25*sim.T) & Tpre<(0.5*sim.T));
Tpre3 = Tpre(Tpre >= (0.5*sim.T) & Tpre<(0.75*sim.T));
Tpre4 = Tpre(Tpre >= (0.75*sim.T));

Tpre1_spk_vec = zeros(1,sim.vecN);Tpre1_spk_vec(round(Tpre1/data.dt))=1;
Tpre2_spk_vec = zeros(1,sim.vecN);Tpre2_spk_vec(round(Tpre2/data.dt))=1;
Tpre3_spk_vec = zeros(1,sim.vecN);Tpre3_spk_vec(round(Tpre3/data.dt))=1;
Tpre4_spk_vec = zeros(1,sim.vecN);Tpre4_spk_vec(round(Tpre4/data.dt))=1;


[d1,~] = corr_fast_v3(Tpre1, Tpost,-.02,.02,40);
[d2,~] = corr_fast_v3(Tpre2, Tpost,-.02,.02,40);
[d3,~] = corr_fast_v3(Tpre3, Tpost,-.02,.02,40);
[d4,~] = corr_fast_v3(Tpre4, Tpost,-.02,.02,40);

[d1_fit, lag1_fit] = xcorr(Tpre1_spk_vec, lam, 20);
[d2_fit, lag2_fit] = xcorr(Tpre2_spk_vec, lam, 20);
[d3_fit, lag3_fit] = xcorr(Tpre3_spk_vec, lam, 20);
[d4_fit, lag4_fit] = xcorr(Tpre4_spk_vec, lam, 20);


yULim = ceil(max([d1; d2; d3; d4])/25)*25;

corrT = figure;
subplot(1,4,1)
hold on
bar(linspace(-.02,.02,40),d1,'k');
plot(-lag1_fit*data.dt, d1_fit, 'r', 'LineWidth',2)
title('<Q_1 of T')
xlim([-.01,.02]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off

subplot(1,4,2)
hold on
bar(linspace(-.02,.02,40),d2,'k');
plot(-lag2_fit*data.dt, d2_fit, 'r', 'LineWidth',2)
title('Q_1 to Q_2 of T')
xlim([-.01,.02]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off


subplot(1,4,3)
hold on
bar(linspace(-.02,.02,40),d3,'k');
plot(-lag3_fit*data.dt, d3_fit, 'r', 'LineWidth',2)
title('Q_2 to Q_3 of T')
xlim([-.01,.02]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off

subplot(1,4,4)
hold on
bar(linspace(-.02,.02,40),d4,'k');
plot(-lag4_fit*data.dt, d4_fit, 'r', 'LineWidth',2)
title('> Q_3 of T')
xlim([-.01,.02]);ylim([0, yULim]);
xlabel('s')
ylabel('Count')
hold off

saveas(corrT, '5_corrT.svg')
saveas(corrT, '5_corrT.png')

