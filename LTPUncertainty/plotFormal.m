plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\fig4_plot\uncertainty';
cd(plotFolder)

close all;
% % baseline
% baselinePlot(simDep, fitDepAll, 'depAll_baseline')
% baselinePlot(simFac, fitFacAll, 'facAll_baseline')
% baselinePlot(simNull, fitNullAll, 'nullAll_baseline')
% 
% %% LTP
% ltpPlot(simDep, fitDepAll, 'depAll_ltp')
% ltpPlot(simDep, fitDepMore, 'depMore_ltp')
% 
% ltpPlot(simFac, fitFacAll, 'facAll_ltp')
% ltpPlot(simFac, fitFacMore, 'facMore_ltp')
% 
% ltpPlot(simNull, fitNullAll, 'nullAll_ltp')
% ltpPlot(simNull, fitNullMore, 'nullMore_ltp')
% 
% %% STP
% stpPlot(simDep, fitDepAll, 'depAll_stp')
% stpPlot(simFac, fitFacAll, 'facAll_stp')
% stpPlot(simNull, fitNullAll, 'nullAll_stp')

%% firing rate
fireRatePlot(simDep, dataDep, 'dep_fireRate')
fireRatePlot(simFac, dataFac, 'fac_fireRate')
fireRatePlot(simNull, dataNull, 'null_fireRate')










%% functions
function baselinePlot(sim, fit, fileName)

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
saveas(baseLine, strcat(fileName, '.svg'))
saveas(baseLine, strcat(fileName, '.png'))

end


function ltpPlot(sim, fit, fileName)

ltp = figure;
idx = 1:size(fit.beta0);
hold on
plot(idx, fit.wt_long, 'r', 'LineWidth', 3)
plot(idx, sim.wt_long, 'k', 'LineWidth', 3)
plot(idx, fit.wt_long + sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
plot(idx, fit.wt_long - sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
ylim([min(sim.wt_long)-2 max(sim.wt_long)+2])
xlim([0 sim.T/sim.dt])
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(ltp,'PaperUnits','inches','PaperPosition',[0 0 5 3])
saveas(ltp, strcat(fileName, '.svg'))
saveas(ltp, strcat(fileName, '.png'))

end


function stpPlot(sim, fit, fileName)

stp = figure;
idx = 1:size(fit.beta0);

se_wt_short = zeros(size(fit.stp_X, 1), 1);
for k = 1:size(fit.stp_X, 1)
    se_wt_short(k) = sqrt(fit.stp_X(k,:)*fit.covB*(fit.stp_X(k,:))');
end
sim.wt_short = 1 + sim.stp_X*sim.stp_B;

hold on
plot(idx, fit.wt_short, 'r', 'LineWidth', 3)
plot(idx, sim.wt_short, 'k', 'LineWidth', 3)
plot(idx, fit.wt_short +  se_wt_short, 'r:', 'LineWidth', 2)
plot(idx, fit.wt_short -  se_wt_short, 'r:', 'LineWidth', 2)
ylim([min(fit.wt_short)-2 max(fit.wt_short)+2])
xlim([0 sim.T/sim.dt])
xticks([0 2 4 6 8 10 12]*1e5)
xticklabels({'0','200','400','600','800','1000','1200'})
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(stp,'PaperUnits','inches','PaperPosition',[0 0 5 3])
saveas(stp, strcat(fileName, '.svg'))
saveas(stp, strcat(fileName, '.png'))

end

function fireRatePlot(sim, data, fileName)


firRate = figure;
hold on
plot(linspace(0,sim.T,sim.vecN), filter(ones(2000,1),1,data.pre_spk_vec)/2, 'b', 'LineWidth', 1.5);
plot(linspace(0,sim.T,sim.vecN), filter(ones(2000,1),1,data.post_spk_vec)/2, 'k', 'LineWidth', 1.5);
hold off
ylim([0 30])
xlim([0 sim.T])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(firRate,'PaperUnits','inches','PaperPosition',[0 0 5 3])
saveas(firRate, strcat(fileName, '.svg'))
saveas(firRate, strcat(fileName, '.png'))

end



