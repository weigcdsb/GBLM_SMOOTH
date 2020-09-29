addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%%
plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\Q';
cd(plotFolder)

%% log-likelihood
% baseline
baseLlhd = figure;
semilogx(Qvec,qbllhd_pred, 'k', 'LineWidth', 2)
xline(Q_true(1, 1), 'Color', 'r', 'LineWidth', 2);
hold on
semilogx(Qvec(4), qbllhd_pred(4), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b')
semilogx(Qvec(18), qbllhd_pred(18), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b')
semilogx(Qopt(1, 1), qbllhd_pred(Qvec == Qopt(1, 1)), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b')
xlim([5*1e-10 5*1e-3])
ylim([min(qbllhd_pred) - 1e2 max(qbllhd_pred) + 1e2])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(baseLlhd,'PaperUnits','inches','PaperPosition',[0 0 6 10])
saveas(baseLlhd, '1_baseLlhd.svg')
saveas(baseLlhd, '1_baseLlhd.png')

% wt_long
ltpLlhd = figure;
semilogx(Qvec,qwllhd_pred, 'k', 'LineWidth', 2)
xline(Q_true(2, 2), 'Color', 'r', 'LineWidth', 2);
hold on
semilogx(Qvec(4), qwllhd_pred(4), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b')
semilogx(Qvec(18), qwllhd_pred(18), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b')
semilogx(Qopt(1, 1), qwllhd_pred(Qvec == Qopt(2, 2)), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b')
xlim([5*1e-10 5*1e-3])
ylim([min(qwllhd_pred) - 1e1 max(qwllhd_pred) + 1e1])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(ltpLlhd,'PaperUnits','inches','PaperPosition',[0 0 6 10])
saveas(ltpLlhd, '2_ltpLlhd.svg')
saveas(ltpLlhd, '2_ltpLlhd.png')

%% baseline
basePlot(sim, fit, '3_baseLine_mle')
basePlot(sim, fitSmall, '4_baseLine_small')
basePlot(sim, fitLarge, '5_baseLine_large')

%% wt_long
ltpPlot(sim, fit, '6_ltp_mle')
ltpPlot(sim, fitSmall, '7_ltp_small')
ltpPlot(sim, fitLarge, '8_ltp_large')

%% modification function
modPlotLoc(sim, fit, '9_mod_mle')
modPlotLoc(sim, fitSmall, '10_mod_small')
modPlotLoc(sim, fitLarge, '11_mod_large')


%%
function basePlot(sim, fit, fileName)
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


function modPlotLoc(sim, fit, fileName)

modFun = figure;
hold on
modPlot(sim, fit, 3, 2)
ylim([0.7 1.1])
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(modFun,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(modFun, strcat(fileName, '.svg'))
saveas(modFun, strcat(fileName, '.png'))
end




