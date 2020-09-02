plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\uncertainty';
cd(plotFolder)

Tpre = data.pre_spk_times;
Tpost = data.post_spk_times(:, 1);
data.vecN = length(data.pre_spk_vec);

lam = exp(fit.beta0 + fit.wt_long.*fit.wt_short.*fit.Xc +...
    fit.hist*fit.hist_beta)*data.dt;


%% firing rate
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

set(firRate,'PaperUnits','inches','PaperPosition',[0 0 6 3])
saveas(firRate, '1_firRate.svg')
saveas(firRate, '1_firRate.png')

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

set(ltp,'PaperUnits','inches','PaperPosition',[0 0 6 3])
saveas(ltp, '2_ltp.svg')
saveas(ltp, '2_ltp.png')







