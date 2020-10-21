addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));
%% pre-firing rate
clc; clear all; close all;
rng(4)
Q_true = diag([1e-6 1e-6]);
dt = 1e-3;

nRate = 6;
nSeed = 50;
preRate_seq = linspace(2, 22, nRate)*dt;
QbEst = zeros(nRate, nSeed);
QwEst = zeros(nRate, nSeed);
seed_seq = randperm(100, nSeed);


for j = 1:nRate
    %     seed_seq = randperm(100, nSeed);
    for k = 1:nSeed
        fprintf('j= %i...', j)
        fprintf('k= %i...', k)
        fprintf('\n')
        
        sim.seed = seed_seq(k);
        sim.T = 5*60;
        sim.dt = 1e-3;
        sim.vecN = round(sim.T/sim.dt);
        sim.pPreSpike = preRate_seq(j);
        sim.alpha_dt = 0.004;
        sim.alpha_tau = 0.001;
        sim.stp_Nq = 5;
        sim.stp_Nm = 450;
        sim.stp_Ns = 50;
        sim.stp_tau= 1;
        sim.stp_B = [0 0 1 3 1]'*(-0.05);
        sim.hist_tau = .01;
        sim.hist_beta = -2;
        
        sim.beta0 = ones(1, sim.T/sim.dt)'*3 + ...
            detrend(cumsum(randn(sim.vecN,1)*sqrt(Q_true(1, 1))));
        sim.wt_long = ones(1, sim.T/sim.dt)'*3 + ...
            detrend(cumsum(randn(sim.vecN,1)*sqrt(Q_true(2, 2))));
        
        data.dt = sim.dt;
        [data,sim] = sim_model(data,sim);
        
        % Q tune
        %         [~, ~, Qopt] = ...
        %             tune_smooth_gblm_1d_golden(data.pre_spk_vec, data.post_spk_vec,...
        %             'iter',15, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta,...
        %             'doFit', false);
        %         [~, ~, ~, ~, ~, Qopt] = ...
        %             tune_smooth_gblm_1d_grid(data.pre_spk_vec, data.post_spk_vec,...
        %             'nq', 10,'hist_tau', sim.hist_tau,...
        %             'hist_beta', sim.hist_beta, 'doFit', false);
        %
        [~, ~, Qopt] = ...
            tune_smooth_gblm_2d_grad(data.pre_spk_vec, data.post_spk_vec,...
            'iter',15, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta,...
            'doFit', false, 'QUB', [1e-4 1e-4]);
        
        %         [~, ~, ~, ~, Qopt] = ...
        %             tune_smooth_gblm_2d_grid(data.pre_spk_vec, data.post_spk_vec,...
        %             'iter',15, 'nq', 10,'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta,...
        %             'doFit', false);
        
        
        QbEst(j, k) = Qopt(1);
        QwEst(j, k) = Qopt(2);
    end
end

clear data sim
save('QpreRate_2d_grad_5min.mat')

plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\Q\2d';
cd(plotFolder)

preRatet_Qb = figure; 
hold on
for m = 1:nRate
    plot(preRate_seq(m)*1e3,log10(QbEst(m,:)), 'o', 'Color', 'k',...
        'LineWidth', 2, 'markerfacecolor', 'k')
end
yline(log10(Q_true(1, 1)), 'r--', 'LineWidth', 2);
xlim([min(preRate_seq*1e3)-3 max(preRate_seq*1e3)+3])
ylim([log10(Q_true(1, 1)) - 4 log10(Q_true(1, 1)) + 4])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(preRatet_Qb,'PaperUnits','inches','PaperPosition',[0 0 6 5])
saveas(preRatet_Qb, 'A3_preRatet_Qb_1min.svg')
saveas(preRatet_Qb, 'A3_preRatet_Qb_1min.png')

preRatet_Qw = figure;
hold on
for m = 1:nRate
    plot(preRate_seq(m)*1e3,log10(QwEst(m,:)), 'o', 'Color', 'k',...
        'LineWidth', 2, 'markerfacecolor', 'k')
end
yline(log10(Q_true(2, 2)), 'r--', 'LineWidth', 2);
xlim([min(preRate_seq*1e3)-3 max(preRate_seq*1e3)+3])
ylim([log10(Q_true(2, 2)) - 4 log10(Q_true(2, 2)) + 4])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(preRatet_Qw,'PaperUnits','inches','PaperPosition',[0 0 6 5])
saveas(preRatet_Qw, 'A3_preRatet_Qw_1min.svg')
saveas(preRatet_Qw, 'A3_preRatet_Qw_1min.png')

%% post-firing rate
clc; clear all; close all;
rng(5)
Q_true = diag([1e-6 1e-6]);
dt = 1e-3;

nBase = 6;
nSeed = 50;

base_seq = log(linspace(15, 40, nBase));
QbEst = zeros(nBase, nSeed);
QwEst = zeros(nBase, nSeed);
seed_seq = randperm(100, nSeed);

for j = 1:nBase
    %     seed_seq = randperm(100, nSeed);
    for k = 1:nSeed
        fprintf('j= %i...', j)
        fprintf('k= %i...', k)
        fprintf('\n')
        
        sim.seed = seed_seq(k);
        sim.T = 5*60;
        sim.dt = 1e-3;
        sim.vecN = round(sim.T/sim.dt);
        sim.pPreSpike = 5*sim.dt;
        sim.alpha_dt = 0.004;
        sim.alpha_tau = 0.001;
        sim.stp_Nq = 5;
        sim.stp_Nm = 450;
        sim.stp_Ns = 50;
        sim.stp_tau= 1;
        sim.stp_B = [0 0 1 3 1]'*(-0.05);
        sim.hist_tau = .01;
        sim.hist_beta = -2;
        
        sim.beta0 = ones(1, sim.T/sim.dt)'*base_seq(j) + ...
            detrend(cumsum(randn(sim.vecN,1)*sqrt(Q_true(1, 1))));
        sim.wt_long = ones(1, sim.T/sim.dt)'*3 + ...
            detrend(cumsum(randn(sim.vecN,1)*sqrt(Q_true(2, 2))));
        
        data.dt = sim.dt;
        [data,sim] = sim_model(data,sim);
        
        %         Q tune
        %         [~, ~, Qopt] = ...
        %             tune_smooth_gblm_1d_golden(data.pre_spk_vec, data.post_spk_vec,...
        %             'iter',15, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta,...
        %             'doFit', false);
        %         [~, ~, ~, ~, ~, Qopt] = ...
        %             tune_smooth_gblm_1d_grid(data.pre_spk_vec, data.post_spk_vec,...
        %             'nq', 10,'hist_tau', sim.hist_tau,...
        %             'hist_beta', sim.hist_beta, 'doFit', false);
        %
        
        [~, ~, Qopt] = ...
            tune_smooth_gblm_2d_grad(data.pre_spk_vec, data.post_spk_vec,...
            'iter',15, 'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta,...
            'doFit', false, 'QUB', [1e-4 1e-4]);
        %
        %         [~, ~, ~, ~, Qopt] = ...
        %             tune_smooth_gblm_2d_grid(data.pre_spk_vec, data.post_spk_vec,...
        %             'iter',15, 'nq', 10,'hist_tau', sim.hist_tau, 'hist_beta', sim.hist_beta,...
        %             'doFit', false);
        
        QbEst(j, k) = Qopt(1);
        QwEst(j, k) = Qopt(2);
    end
end

clear data sim
save('QpostRate_2d_grad_5min.mat')

postRatet_Qb = figure; 
hold on
for m = 1:nBase
   plot(exp(base_seq(m)),log10(QbEst(m,:)), 'o', 'Color', 'k',...
    'LineWidth', 2, 'markerfacecolor', 'k')
end
xlim([exp(min(base_seq))-3 exp(max(base_seq))+3])
ylim([log10(Q_true(1, 1)) - 4 log10(Q_true(1, 1)) + 4])
yline(log10(Q_true(1, 1)), 'r--', 'LineWidth', 2);
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(postRatet_Qb,'PaperUnits','inches','PaperPosition',[0 0 6 5])
saveas(postRatet_Qb, 'A4_postRatet_Qb_1min.svg')
saveas(postRatet_Qb, 'A4_postRatet_Qb_1min.png')

postRatet_Qw = figure; 
hold on
for m = 1:nBase
   plot(exp(base_seq(m)),log10(QwEst(m,:)), 'o', 'Color', 'k',...
    'LineWidth', 2, 'markerfacecolor', 'k')
end
xlim([exp(min(base_seq))-3 exp(max(base_seq))+3])
ylim([log10(Q_true(2, 2)) - 4 log10(Q_true(2, 2)) + 4])
yline(log10(Q_true(2, 2)), 'r--', 'LineWidth', 2);
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(postRatet_Qw,'PaperUnits','inches','PaperPosition',[0 0 6 5])
saveas(postRatet_Qw, 'A4_postRatet_Qw_1min.svg')
saveas(postRatet_Qw, 'A4_postRatet_Qw_1min.png')

