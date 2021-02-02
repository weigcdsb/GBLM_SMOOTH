addpath(genpath('~/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('~/GitHub/GBLM_SMOOTH/core'));

%% set up true parameters
clc;clear all;close all;
T = 20*60;
dt = 0.001;

trueParam = [0 0 1 3 1]'*(-0.05);
beta0 = ones(1, T/dt)'*3;
wt_long = [repmat(2.5, 1, round(T/(dt*2))) repmat(3.5, 1, T/dt - round(T/(dt*2)))]';

%% run...
sim.seed = 11;
sim.T = T;
sim.dt = dt;
sim.vecN = round(sim.T/sim.dt);
sim.pPreSpike = 5*sim.dt;
sim.alpha_dt = 0.004;
sim.alpha_tau = 0.001;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1;
sim.stp_B = trueParam;
sim.hist_tau = .01;
sim.hist_beta = -2;

sim.beta0 = beta0;
sim.wt_long = wt_long;
data.dt = sim.dt;
data.T = sim.T;
data.vecN = sim.vecN;
[data,sim] = sim_model(data,sim);


%% plot
Tpre = data.pre_spk_times;
Tpost = data.post_spk_times(:, 1);
lam = sim.lam*sim.dt;

plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig1_plot';
cd(plotFolder)

%% wt_long (simulation)
wtLong = figure;
hold on
plot(1 : sim.T/(sim.dt*2), sim.wt_long(1:end/2),...
    'Color', [0, 0.4470, 0.7410], 'LineWidth', 3)
plot(sim.T/(sim.dt*2)+1 : sim.T/sim.dt, sim.wt_long(end/2 + 1:end),...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
plot(sim.T/(sim.dt*2): sim.T/(sim.dt*2)+1, sim.wt_long(end/2: end/2+1),...
    'k', 'LineWidth', 3)
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
ylim([min(sim.wt_long)-1 max(sim.wt_long)+1]);
xlim([0 sim.T/sim.dt])
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(wtLong,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(wtLong, '1_LTP.svg')
saveas(wtLong, '1_LTP.png')

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
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(corrOverall,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(corrOverall, '2_corrOverall.svg')
saveas(corrOverall, '2_corrOverall.png')

%% cross-correlogram (before & after)
qspk1 = Tpre(Tpre < 0.5*sim.T);
qspk2 = Tpre(Tpre >= 0.5*sim.T);
tvec = linspace(-0.02,0.02,102);
tvec = tvec+mean(diff(tvec))/2;

qvec1 = zeros(1,data.vecN);
qvec1(round(qspk1/data.dt))=1;
qvec2 = zeros(1,data.vecN);
qvec2(round(qspk2/data.dt))=1;

[d1,~] = corr_fast_v3(qspk1, Tpost,-.02,.02,102);
[d2,~] = corr_fast_v3(qspk2, Tpost,-.02,.02,102);
[d1_fit, lag1_fit] = xcorr(qvec1, lam, 20);
[d2_fit, lag2_fit] = xcorr(qvec2, lam, 20);


corrBefore = figure;
hold on
bar(tvec(1:end-1)*1e3,d1(1:end-1),1,'k','EdgeColor','none');
plot(-lag1_fit*data.dt*1e3, d1_fit*mean(diff(tvec))/data.dt,...
    'Color', [0, 0.4470, 0.7410], 'LineWidth',3)
xlim([-.01 .02]*1e3);
ylim([0 300]);
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(corrBefore,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(corrBefore, '3_corrBefore.svg')
saveas(corrBefore, '3_corrBefore.png')

corrAfter = figure;
hold on
bar(tvec(1:end-1)*1e3,d2(1:end-1),1,'k','EdgeColor','none');
plot(-lag2_fit*data.dt*1e3, d2_fit*mean(diff(tvec))/data.dt,...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth',3)
xlim([-.01 .02]*1e3);
ylim([0 300]);
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(corrAfter,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(corrAfter, '4_corrAfter.svg')
saveas(corrAfter, '4_corrAfter.png')

%% efficacy (before & after)
Tpre1 = Tpre(Tpre < 0.5*sim.T);
Tpre2 = Tpre(Tpre >= 0.5*sim.T);

isi1 = [Inf; diff(Tpre1)];
isi2 = [Inf; diff(Tpre2)];
quantiles1 = prctile(isi1,linspace(0,100,25));
quantiles2 = prctile(isi2,linspace(0,100,25));
qx1 = quantiles1(1:end-1)+diff(quantiles1)/2;
qx2 = quantiles2(1:end-1)+diff(quantiles2)/2;


efficacy1 = zeros(1, length(quantiles1)-1);
efficacy2 = zeros(1, length(quantiles2)-1);
model_efficacy1 = zeros(1, length(quantiles1)-1);
model_efficacy2 = zeros(1, length(quantiles2)-1);
syn_kern = sim.syn(linspace(0, 1, 1/data.dt));
sk1 = exp(syn_kern);
sk1=sk1(sk1>1);

wt = (1+sim.stp_X*sim.stp_B).*sim.wt_long;
for q=1:length(quantiles1)-1
    qspk1 = Tpre1(isi1>=quantiles1(q) & isi1<quantiles1(q+1));
    qspk2 = Tpre2(isi2>=quantiles2(q) & isi2<quantiles2(q+1));
    
    d1base = corr_fast_v3(qspk1,Tpost,-.008,-0.004,20);
    d1 = corr_fast_v3(qspk1,Tpost,0.004,.008,20);
    d2base = corr_fast_v3(qspk1,Tpost,-.008,-0.004,20);
    d2 = corr_fast_v3(qspk2,Tpost,0.004,.008,20);
    
    n1 = zeros(1,data.T/data.dt);
    n1(ceil(qspk1/data.dt))=1;
    model_efficacy1(q) = mean(exp(wt(n1>0)));
    
    n2 = zeros(1,data.T/data.dt);
    n2(ceil(qspk2/data.dt))=1;
    model_efficacy2(q) = mean(exp(wt(n2>0)));
    
    efficacy1(q) = (sum(d1) - sum(d1base))/length(qspk1);
    efficacy2(q) = (sum(d2) - sum(d2base))/length(qspk2);
end

effPlot = figure;
clf
hold on
plot(qx1*1000, efficacy1, 'o', 'Color', [0, 0.4470, 0.7410],...
    'LineWidth', 2, 'markerfacecolor', [0, 0.4470, 0.7410])
plot(qx1*1000, (model_efficacy1)/nanmean(model_efficacy1)*nanmean(efficacy1),...
    'Color', [0, 0.4470, 0.7410], 'LineWidth', 3)
plot(qx2*1000, efficacy2, 'o', 'Color', [0.8500, 0.3250, 0.0980],...
    'LineWidth', 2, 'markerfacecolor', [0.8500, 0.3250, 0.0980])
plot(qx2*1000, (model_efficacy2)/nanmean(model_efficacy2)*nanmean(efficacy2),...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(effPlot,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(effPlot, '5_effPlot.svg')
saveas(effPlot, '5_effPlot.png')

%% cross-correlogram split by ISI (before & after)
cd(strcat(plotFolder, '\corrISI'))

Tpre1 = Tpre(Tpre < 0.5*sim.T);
Tpre2 = Tpre(Tpre >= 0.5*sim.T);

isi1 = [Inf; diff(Tpre1)];
isi2 = [Inf; diff(Tpre2)];
quantiles1 = prctile(isi1,linspace(0,100,5));
quantiles2 = prctile(isi2,linspace(0,100,5));

for q=1:length(quantiles1)-1
    qspk1 = Tpre1(isi1>=quantiles1(q) & isi1<quantiles1(q+1));
    qspk2 = Tpre2(isi2>=quantiles2(q) & isi2<quantiles2(q+1));
    
    qvec1 = zeros(1,data.vecN);
    qvec1(round(qspk1/data.dt))=1;
    qvec2 = zeros(1,data.vecN);
    qvec2(round(qspk2/data.dt))=1;
    
    d1 = corr_fast_v3(qspk1,Tpost,-.025,.025,64);
    d2 = corr_fast_v3(qspk2,Tpost,-.025,.025,64);
    tvec = linspace(-0.025,0.025,64);
    tvec = tvec+mean(diff(tvec))/2;
    [dfit1,lag1_fit] = xcorr(qvec1, lam, 25);
    [dfit2,lag2_fit] = xcorr(qvec2, lam, 25);
    
    corrISIBefore = figure;
    hold on
    bar(tvec(1:end-1)*1e3,d1(1:end-1),1,'k','EdgeColor','none');
    plot(-lag1_fit*data.dt*1e3, dfit1*mean(diff(tvec))/data.dt,...
        'Color', [0, 0.4470, 0.7410], 'LineWidth',3)
    ylim([0 300])
    hold off
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    
    set(corrISIBefore,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    saveas(corrISIBefore, strcat('corrISIBefore_Q', string(q), '.svg'))
    saveas(corrISIBefore, strcat('corrISIBefore_Q', string(q), '.png'))
    
    corrISIAfter = figure;
    hold on
    bar(tvec(1:end-1)*1e3,d2(1:end-1),1,'k','EdgeColor','none');
    plot(-lag2_fit*data.dt*1e3, dfit2*mean(diff(tvec))/data.dt,...
        'Color', [0.8500, 0.3250, 0.0980], 'LineWidth',3)
    ylim([0 300])
    hold off
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    
    set(corrISIAfter,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    saveas(corrISIAfter, strcat('corrISIAfter_Q', string(q), '.svg'))
    saveas(corrISIAfter, strcat('corrISIAfter_Q', string(q), '.png'))
end




