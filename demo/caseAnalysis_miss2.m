addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%%
clc;clear all;close all;

plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\miss\noWt_short';
cd(plotFolder)

%%
Tpre = data.pre_spk_times;
Tpost = data.post_spk_times(:, 1);
data.vecN = length(data.pre_spk_vec);

lam = exp(fit.beta0 + fit.wt_long.*fit.Xc +...
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
plot(idx, sim.wt_long, 'k', 'LineWidth', 3)
plot(idx, fit.wt_long, 'r', 'LineWidth', 3)
plot(idx, fit.wt_long + sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
plot(idx, fit.wt_long - sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
ylim([min(sim.wt_long)-3 max(sim.wt_long)+3])
xlim([0 sim.T/sim.dt])
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(ltp,'PaperUnits','inches','PaperPosition',[0 0 5 3])
saveas(ltp, '2_ltp.svg')
saveas(ltp, '2_ltp.png')

ltp = figure;
idx = 1:size(fit.beta0);
hold on
plot(idx, sim.wt_long.*(1 + sim.stp_X*sim.stp_B), 'k', 'LineWidth', 3)
plot(idx, fit.wt_long, 'r', 'LineWidth', 3)
plot(idx, fit.wt_long + sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
plot(idx, fit.wt_long - sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
ylim([min(sim.wt_long)-3 max(sim.wt_long)+3])
xlim([0 sim.T/sim.dt])
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(ltp,'PaperUnits','inches','PaperPosition',[0 0 5 3])
saveas(ltp, '2_ltpv2.svg')
saveas(ltp, '2_ltpv2.png')

%% overall cross-correlogram
[d,~] = corr_fast_v3(Tpre, Tpost,-.02,.02,102);
[d_fit, lag_fit] = xcorr(data.pre_spk_vec, lam, 20);

tvec = linspace(-0.02,0.02,102);
tvec = tvec+mean(diff(tvec))/2;

corrOverall = figure;
hold on
bar(tvec(1:end-1)*1e3,d(1:end-1),1,'k','EdgeColor','none');
plot(-lag_fit*data.dt*1e3, d_fit*mean(diff(tvec))/data.dt, 'r', 'LineWidth',3)
xlim([-.01 .02]*1e3);
ylim([0 450])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(corrOverall,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(corrOverall, '3_corrOverall.svg')
saveas(corrOverall, '3_corrOverall.png')


%% split cross-correlogram (STP)
cd(strcat(plotFolder, '\corrISI'))
maxd = -Inf;

isi = [Inf; diff(Tpre)];
quantiles = prctile(isi,linspace(0,100,5));

for q=1:length(quantiles)-1
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    if maxd<max(d),maxd=max(d); end
end

for q=1:length(quantiles)-1
    
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    qvec = zeros(1,data.vecN);
    qvec(round(qspk/data.dt))=1;
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    tvec = linspace(-0.025,0.025,64);
    tvec = tvec+mean(diff(tvec))/2;
    [dfit,lag_fit] = xcorr(qvec, lam, 25);
    
    corrISI = figure;
    bar(tvec(1:end-1),d(1:end-1),1,'k','EdgeColor','none');
    hold on
    plot(-lag_fit*data.dt, dfit*mean(diff(tvec))/data.dt, 'r', 'LineWidth',3)
    %     title(sprintf('%0.2f ms',quantiles(q)*1000))
    %     xlabel('Time (s)')
    %     ylabel('Count')
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    hold off
    %     ylim([0 ceil(maxd/50)*50])
    ylim([0 300])
    xlim([-0.01 0.02])
    
    set(corrISI,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    saveas(corrISI, strcat('corrISI_Q', string(q), '.svg'))
    saveas(corrISI, strcat('corrISI_Q', string(q), '.png'))
end

%% split cross-correlogram (LTP)
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
    %     title(sprintf('%0.2f s', q*0.25*sim.T))
    %     xlabel('ms')
    %     ylabel('Count')
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    hold off
    %     ylim([0 ceil(maxd/50)*50])
    ylim([0 300])
    xlim([-0.01 0.02])
    
    set(corrT,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    saveas(corrT, strcat('corrT_Q', string(q), '.svg'))
    saveas(corrT, strcat('corrT_Q', string(q), '.png'))
end
