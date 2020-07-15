addpath(genpath('D:/cleanTemp/helper'));
addpath(genpath('D:/cleanTemp/core'));

%%
clc;clear all;close all;
% load('F:\COVID-19\updateResults\caseResults\2_depression_constant.mat');

% cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\1_facilitation_constant'
% cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\2_depression_constant'
% cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\3_no_plasticity_constant'
cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\4_depression_linear_jump'
% cd 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case\5_depression_sin_sin'

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
data.vecN = length(data.pre_spk_vec);

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

isi = [Inf; diff(Tpre)];
quantiles = prctile(isi,linspace(0,100,5));
maxd = -Inf;

corrISI = figure;
for q=1:length(quantiles)-1
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    qvec = zeros(1,data.vecN);
    qvec(round(qspk/data.dt))=1;
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
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
    ylim([0 ceil(maxd/50)*50])
end

saveas(corrISI, '4_corrISI.svg')
saveas(corrISI, '4_corrISI.png')

%% show LTP
maxd = -Inf;

corrT = figure;
for q=1:4
    qspk = Tpre(Tpre >= (q-1)*0.25*sim.T & Tpre < q*0.25*sim.T);
    qvec = zeros(1,data.vecN);
    qvec(round(qspk/data.dt))=1;
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    tvec = linspace(-0.025,0.025,64);
    tvec = tvec+mean(diff(tvec))/2;
    [dfit,lag_fit] = xcorr(qvec, lam, 25);
    
    subplot(1,4,q)
    bar(tvec(1:end-1),d(1:end-1),1,'k','EdgeColor','none');
    hold on
    plot(-lag_fit*data.dt, dfit*mean(diff(tvec))/data.dt, 'r', 'LineWidth',2)
    title(sprintf('%0.2f s', q*0.25*sim.T))
    xlabel('ms')
    ylabel('Count')
    hold off
    
    if maxd<max(d),maxd=max(d); end
end

for q=1:4
    subplot(1,4,q)
    ylim([0 ceil(maxd/50)*50])
end

saveas(corrT, '5_corrT.svg')
saveas(corrT, '5_corrT.png')




