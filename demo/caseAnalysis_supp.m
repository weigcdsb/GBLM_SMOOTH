addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/helper'));
addpath(genpath('C:/Users/gaw19004/Desktop/ganchaoResearch/GBLM_SMOOTH-master/core'));

%%
clc;clear all;close all;
load('E:\casedDemoResults_supp_0801\2_depression_conJump.mat');
plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\case_supp\2_depression_conJump\supp';
cd(plotFolder)

Tpre = data.pre_spk_times;
Tpost = data.post_spk_times(:, 1);
data.vecN = length(data.pre_spk_vec);
%% wt_long (simulation)

wtLong = figure;
plot(sim.wt_long, 'k', 'LineWidth', 3)
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
ylim([0 5]);
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
ylim([0 150]);
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
ylim([0 150]);
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off
hold off

saveas(corrBefore, '3_corrBefore.svg')
saveas(corrBefore, '3_corrBefore.png')

corrAfter = figure;
hold on
bar(tvec(1:end-1)*1e3,d2(1:end-1),1,'k','EdgeColor','none');
xlim([-.01 .02]*1e3);
ylim([0 150]);
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off
hold off

saveas(corrAfter, '4_corrAfter.svg')
saveas(corrAfter, '4_corrAfter.png')

%% effect modification function (before & after)

modBefore = figure;
plot((1 + sim.stp_basis'*sim.stp_B)*sim.wt_long(1), 'k', 'LineWidth', 3)
ylim([0 4]);
set(gca,'FontSize',15, 'LineWidth', 1.5)
box off

saveas(modBefore, '5_modBefore.svg')
saveas(modBefore, '5_modBefore.png')

modAfter = figure;
plot((1 + sim.stp_basis'*sim.stp_B)*sim.wt_long(end), 'k', 'LineWidth', 3)
ylim([0 4]);
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
    ylim([0 100])
    set(gca,'FontSize',15, 'LineWidth', 1.5)
    box off
    
    saveas(corrISIBefore, strcat('corrISIBefore_Q', string(q), '.svg'))
    saveas(corrISIBefore, strcat('corrISIBefore_Q', string(q), '.png'))
    
    corrISIAfter = figure;
    bar(tvec(1:end-1)*1e3,d2(1:end-1),1,'k','EdgeColor','none');
    ylim([0 100])
    set(gca,'FontSize',15, 'LineWidth', 1.5)
    box off
    
    saveas(corrISIAfter, strcat('corrISIAfter_Q', string(q), '.svg'))
    saveas(corrISIAfter, strcat('corrISIAfter_Q', string(q), '.png'))
end




