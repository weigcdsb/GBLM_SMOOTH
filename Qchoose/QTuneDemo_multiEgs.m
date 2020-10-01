addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('C:/Users/gaw19004/Documents/GitHub/GBLM_SMOOTH/core'));

%% set up simulation
clc; clear all; close all;
rng(3)

nQ = 4;
nSeed = 5;
Q_true_seq = logspace(-7, -4, nQ);
QbEst = zeros(nQ, nQ, nSeed);
QwEst = zeros(nQ, nQ, nSeed);

for n = 1:nQ
    for m = 1:nQ
        seed_seq = randperm(100, nSeed);
        for k = 1:nSeed
            fprintf('n= %i...', n)
            fprintf('m= %i...', m)
            fprintf('k= %i...', k)
            fprintf('\n')
            
            Q_true = diag([Q_true_seq(n) Q_true_seq(m)]);
            sim.seed = seed_seq(k);
            sim.T = 20*60;
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
            
            sim.beta0 = ones(1, sim.T/sim.dt)'*3 + ...
                detrend(cumsum(randn(sim.vecN,1)*sqrt(Q_true(1, 1))));
            sim.wt_long = ones(1, sim.T/sim.dt)'*3 + ...
                detrend(cumsum(randn(sim.vecN,1)*sqrt(Q_true(2, 2))));
            
            data.dt = sim.dt;
            [data,sim] = sim_model(data,sim);
            
            % Q tune
            [~, ~, ~, ~, ~, Qopt] = ...
                tune_smooth_gblm(data.pre_spk_vec, data.post_spk_vec,...
                'nq', 10,'hist_tau', sim.hist_tau,...
                'hist_beta', sim.hist_beta, 'doFit', false);
            
            QbEst(n, m, k) = Qopt(1, 1);
            QwEst(n, m, k) = Qopt(2, 2);
            
        end
    end
end

save('Q_multi2.mat')

%%
cd('C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\Q')

LCol = [0.9290, 0.6940, 0.1250]; % yellow
UCol = [0.25, 0.25, 0.25];
colGrad = [linspace(LCol(1),UCol(1),nQ)',...
    linspace(LCol(2),UCol(2),nQ)', linspace(LCol(3),UCol(3),nQ)'];

baseQ = figure;
hold on
for i = 1:nQ
    for j = 1:nSeed
        jitters = normrnd(0, 0.1, 1, nQ);
        plot(log10(Q_true_seq) + jitters, log10(QbEst(:, i, j)), 'o',...
            'Color', colGrad(i,:),'LineWidth', 2, 'markerfacecolor',colGrad(i,:))
    end
end
plot([-10, -1], [-10, -1], 'r--', 'LineWidth', 2)
xlim([-10 -1])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(baseQ,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(baseQ, 'A1_baseQ.svg')
saveas(baseQ, 'A1_baseQ.png')

%%

for i = 1:nSeed
   QwEst(:, :, i) =  QwEst(:, :, i)';
end

LCol = [0.7490, 0.5608, 0.8902]; % purple
UCol = [0.25, 0.25, 0.25];
colGrad = [linspace(LCol(1),UCol(1),nQ)',...
    linspace(LCol(2),UCol(2),nQ)', linspace(LCol(3),UCol(3),nQ)'];

wtQ = figure;
hold on
for i = 1:nQ
    for j = 1:nSeed
        jitters = normrnd(0, 0.1, 1, nQ);
        plot(log10(Q_true_seq) + jitters, log10(QwEst(:, i, j)), 'o',...
            'Color', colGrad(i,:),'LineWidth', 2, 'markerfacecolor',colGrad(i,:))
    end
end
plot([-10, -1], [-10, -1], 'r--', 'LineWidth', 2)
xlim([-10 -1])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(wtQ,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(wtQ, 'A2_wtQ.svg')
saveas(wtQ, 'A2_wtQ.png')
