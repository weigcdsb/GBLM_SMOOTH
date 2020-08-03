addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/helper'));
addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/core'));

%% set up true parameters
clc;clear all;close all;
T = 20*60;
dt = 0.001;

trueParam = [0 0 1 3 1]'*(-0.05);
beta0 = ones(1, T/dt)'*3;
wt_long = [repmat(2.5, 1, round(T/(dt*2))) repmat(4.5, 1, T/dt - round(T/(dt*2)))]';

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
[data,sim] = sim_model(data,sim);
%% plot
Tpre = data.pre_spk_times;
Tpost = data.post_spk_times(:, 1);

plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig1_plot';
cd(plotFolder)

%% wt_long (simulation)
wtLong = figure;
plot(sim.wt_long, 'k', 'LineWidth', 3)
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
ylim([min(sim.wt_long)-2 max(sim.wt_long)+2]);
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off


saveas(wtLong, '1_LTP.svg')
saveas(wtLong, '1_LTP.png')

%% overall cross-correlogram
[d,~] = corr_fast_v3(Tpre, Tpost,-.02,.02,102);
tvec = linspace(-0.02,0.02,102);
tvec = tvec+mean(diff(tvec))/2;

corrOverall = figure;
hold on
bar(tvec(1:end-1)*1e3,d(1:end-1),1,'k','EdgeColor','none');
xlim([-.01 .02]*1e3);
ylim([0 500]);
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off
hold off

saveas(corrOverall, '2_corrOverall.svg')
saveas(corrOverall, '2_corrOverall.png')

%% cross-correlogram (before & after)
qspk1 = Tpre(Tpre < 0.5*sim.T);
qspk2 = Tpre(Tpre >= 0.5*sim.T);
tvec = linspace(-0.02,0.02,102);
tvec = tvec+mean(diff(tvec))/2;

[d1,~] = corr_fast_v3(qspk1, Tpost,-.02,.02,102);
[d2,~] = corr_fast_v3(qspk2, Tpost,-.02,.02,102);

corrBefore = figure;
hold on
bar(tvec(1:end-1)*1e3,d1(1:end-1),1,'k','EdgeColor','none');
xlim([-.01 .02]*1e3);
ylim([0 400]);
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off
hold off

saveas(corrBefore, '3_corrBefore.svg')
saveas(corrBefore, '3_corrBefore.png')

corrAfter = figure;
hold on
bar(tvec(1:end-1)*1e3,d2(1:end-1),1,'k','EdgeColor','none');
xlim([-.01 .02]*1e3);
ylim([0 400]);
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off
hold off

saveas(corrAfter, '4_corrAfter.svg')
saveas(corrAfter, '4_corrAfter.png')

%% effect modification function (before & after)

modBefore = figure;
plot((1 + sim.stp_basis'*sim.stp_B)*sim.wt_long(1), 'k', 'LineWidth', 3)
ylim([0 5]);
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off

saveas(modBefore, '5_modBefore.svg')
saveas(modBefore, '5_modBefore.png')

modAfter = figure;
plot((1 + sim.stp_basis'*sim.stp_B)*sim.wt_long(end), 'k', 'LineWidth', 3)
ylim([0 5]);
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off

saveas(modAfter, '6_modAfter.svg')
saveas(modAfter, '6_modAfter.png')

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
    
    d1 = corr_fast_v3(qspk1,Tpost,-.025,.025,64);
    d2 = corr_fast_v3(qspk2,Tpost,-.025,.025,64);
    tvec = linspace(-0.025,0.025,64);
    tvec = tvec+mean(diff(tvec))/2;
    
    corrISIBefore = figure;
    bar(tvec(1:end-1)*1e3,d1(1:end-1),1,'k','EdgeColor','none');
    ylim([0 300])
    set(gca,'FontSize',15, 'LineWidth', 1.5)
    box off
    
    saveas(corrISIBefore, strcat('corrISIBefore_Q', string(q), '.svg'))
    saveas(corrISIBefore, strcat('corrISIBefore_Q', string(q), '.png'))
    
    corrISIAfter = figure;
    bar(tvec(1:end-1)*1e3,d2(1:end-1),1,'k','EdgeColor','none');
    ylim([0 300])
    set(gca,'FontSize',15, 'LineWidth', 1.5)
    box off
    
    saveas(corrISIAfter, strcat('corrISIAfter_Q', string(q), '.svg'))
    saveas(corrISIAfter, strcat('corrISIAfter_Q', string(q), '.png'))
end










