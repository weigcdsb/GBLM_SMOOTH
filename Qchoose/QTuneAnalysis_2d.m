% addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
% addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));
addpath(genpath('D:\GitHub\GBLM_SMOOTH\helper'));
addpath(genpath('D:\GitHub\GBLM_SMOOTH\core'));


%%
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\Q\2d';
plotFolder = 'D:\GitHub\GBLM_SMOOTH\plot\Q\2d';
cd(plotFolder)

%% llhd plot
% heatmap version
llhdPlot_heat = figure; 
set(llhdPlot_heat,'color','w');
hold on
xlabel('log_{10}(Q_{LTP})'); ylabel('log_{10}(Q_{baseline})');
colormap(gray(256));
colorbar;
imagesc(log10(Qvec), log10(Qvec), llhdmesh);
xlim([min(log10(Qvec))-.5, max(log10(Qvec))+.5]);
ylim([min(log10(Qvec))-.5, max(log10(Qvec))+.5]);

plot(log10(Qopt_grad(2)), log10(Qopt_grad(1)), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b', 'MarkerSize',5)
plot(log10(Qvec(4)), log10(Qvec(4)), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b', 'MarkerSize',5)
plot(log10(Qvec(18)), log10(Qvec(18)), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b', 'MarkerSize',5)
plot(log10(Q_true(2, 2)), log10(Q_true(1, 1)), '^', 'Color', 'r',...
    'LineWidth', 2, 'markerfacecolor', 'r', 'MarkerSize',6)
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(llhdPlot_heat,'PaperUnits','inches','PaperPosition',[0 0 7 6])
saveas(llhdPlot_heat, '1_llhdPlot_heat.svg')
saveas(llhdPlot_heat, '1_llhdPlot_heat.png')


% surface plot version
data.vecN = length(data.pre_spk_vec);
fit.Q = Q_true;[~, fit, ~] = evalc('loopCore(data, fit)');llhd_true = fit.llhd_pred;
fit.Q = diag([Qopt_grad(1) Qopt_grad(2)]);[~, fit, ~] = evalc('loopCore(data, fit)');llhd_mle = fit.llhd_pred;

llhdPlot_surf = figure; 
set(llhdPlot_surf,'color','w');
colormap(gray(256));
surf(log10(Qvec), log10(Qvec), llhdmesh);
% colorbar;
hold on
xlabel('log_{10}(Q_{LTP})'); ylabel('log_{10}(Q_{baseline})');
zlabel('Log-likelihood')
xlim([min(log10(Qvec))-.5, max(log10(Qvec))+.5]);
ylim([min(log10(Qvec))-.5, max(log10(Qvec))+.5]);
plot3(log10(Qopt_grad(2)), log10(Qopt_grad(1)), llhd_mle, 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b', 'MarkerSize',5)
plot3(log10(Qvec(4)), log10(Qvec(4)), llhdmesh(4, 4), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b', 'MarkerSize',5)
plot3(log10(Qvec(18)), log10(Qvec(18)), llhdmesh(18, 18), 'o', 'Color', 'b',...
    'LineWidth', 2, 'markerfacecolor', 'b', 'MarkerSize',5)
plot3(log10(Q_true(2, 2)), log10(Q_true(1, 1)), llhd_true, '^', 'Color', 'r',...
    'LineWidth', 2, 'markerfacecolor', 'r', 'MarkerSize',6)
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(llhdPlot_surf,'PaperUnits','inches','PaperPosition',[0 0 7 6])
% plot2svg("2_llhdPlot_surf.svg", llhdPlot_surf);
% saveas(llhdPlot_surf, '2_llhdPlot_surf.svg')
exportgraphics(llhdPlot_surf,'2_llhdPlot_surf.eps','ContentType','vector')
saveas(llhdPlot_surf, '2_llhdPlot_surf.png')


%% baseline
basePlot(sim, fitMLE, '3_baseLine_mle')
basePlot(sim, fitSmall, '4_baseLine_small')
basePlot(sim, fitLarge, '5_baseLine_large')

%% wt_long
ltpPlot(sim, fitMLE, '6_ltp_mle')
ltpPlot(sim, fitSmall, '7_ltp_small')
ltpPlot(sim, fitLarge, '8_ltp_large')

%% modification function
modPlotLoc(sim, fitMLE, '9_mod_mle')
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


