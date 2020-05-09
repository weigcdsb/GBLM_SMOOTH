addpath(genpath('D:/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('D:/GitHub/GBLM_SMOOTH/core'));

%%
clc;clear all;close all;
d = load('mseResults_raw.mat');
d = deleteSing(d);

Nq_true = 5; Nq_fit = 5;
Nm = 450; Nst = 50;
Bm0_true = getBasis('rcos', Nq_true, Nm,Nst,0);
Bm0_fit = getBasis('rcos', Nq_fit, Nm,Nst,0);

trueParam_seq = ...
    [[0 1 2 3 4]'*(0.06)...
    [0 -1 -2 -3 -4]'*(0.06)...
    [0 0 0 0 0]'];

numSeed = length(d.seed_seq);

mod_fac = 1 +Bm0_true'*trueParam_seq(:, 1);
mod_dep = 1 +Bm0_true'*trueParam_seq(:, 2);
mod_null = 1 +Bm0_true'*trueParam_seq(:, 3);

% In this case, Bm0_true == Bm0_fit, so we can do in a more efficient way
modRes_fac = trueParam_seq(:, 1) - d.short_param_track_fac;
modRes_dep = trueParam_seq(:, 2) - d.short_param_track_dep;
modRes_null = trueParam_seq(:, 3) - d.short_param_track_null;

nMod = length(Bm0_true);

res_fac = zeros(nMod*numSeed, length(d.Q_wt_long_seq));
res_dep = zeros(nMod*numSeed, length(d.Q_wt_long_seq));
res_null = zeros(nMod*numSeed, length(d.Q_wt_long_seq));

for k = 1:numSeed
    
    res_fac((nMod*(k-1) + 1: nMod*k), :) = Bm0_true'*modRes_fac(:, :, k);
    res_dep((nMod*(k-1) + 1: nMod*k), :) = Bm0_true'*modRes_dep(:, :, k);
    res_null((nMod*(k-1) + 1: nMod*k), :) = Bm0_true'*modRes_null(:, :, k);
    
end

absRes_fac = abs(res_fac);
absRes_dep = abs(res_dep);
absRes_null = abs(res_null);

boxplot(absRes_fac);
boxplot(absRes_dep);
boxplot(absRes_null);

%
org = [0.8500, 0.3250, 0.0980];
blu = [0.3010, 0.7450, 0.9330];
grn = [0.4660, 0.6740, 0.1880];

log10Q = log10(d.Q_wt_long_seq);
log10Q(1) = -10;

%% Mean + STD
% include Q = 1e-3;
hold on
plot(log10Q, mean(absRes_fac, 1), 'Color', org,'LineWidth',2)
plot(log10Q, mean(absRes_dep, 1), 'Color', blu,'LineWidth',2)
plot(log10Q, mean(absRes_null, 1), 'Color', grn,'LineWidth',2)
xticks(log10Q)
xticklabels({'Inf','-7','-6','-5','-4', '-3'})
xlim([-0.5 + min(log10Q)  0.5 + max(log10Q)])
legend('Facilitation', 'Depression', 'No Plasticity')
hold off

hold on
errorbar(log10Q, mean(absRes_fac, 1), std(absRes_fac, 1), 'Color', org,'LineWidth',2)
errorbar(log10Q+0.05, mean(absRes_dep, 1), std(absRes_dep, 1), 'Color', blu,'LineWidth',2)
errorbar(log10Q+0.1, mean(absRes_null, 1), std(absRes_null, 1), 'Color', grn,'LineWidth',2)
xticks(log10Q)
xticklabels({'Inf','-7','-6','-5','-4', '-3'})
xlim([-0.5 + min(log10Q)  0.5 + max(log10Q)])
legend('Facilitation', 'Depression', 'No Plasticity')
hold off


% exclude Q = 1e-3
absRes2_fac = absRes_fac(:, 1:end - 1);
absRes2_dep = absRes_dep(:, 1:end - 1);
absRes2_null = absRes_null(:, 1:end - 1);
log10Q2 = log10Q(1:end-1);

hold on
plot(log10Q2, mean(absRes2_fac, 1), 'Color', org,'LineWidth',2)
plot(log10Q2, mean(absRes2_dep, 1), 'Color', blu,'LineWidth',2)
plot(log10Q2, mean(absRes2_null, 1), 'Color', grn,'LineWidth',2)
xticks(log10Q2)
xticklabels({'Inf','-7','-6','-5','-4'})
xlim([-0.5 + min(log10Q2)  0.5 + max(log10Q2)])
legend('Facilitation', 'Depression', 'No Plasticity')
hold off

hold on
errorbar(log10Q2, mean(absRes2_fac, 1), std(absRes2_fac, 1), 'Color', org,'LineWidth',2)
errorbar(log10Q2+0.05, mean(absRes2_dep, 1), std(absRes2_dep, 1), 'Color', blu,'LineWidth',2)
errorbar(log10Q2+0.1, mean(absRes2_null, 1), std(absRes2_null, 1), 'Color', grn,'LineWidth',2)
xticks(log10Q2)
xticklabels({'Inf','-7','-6','-5','-4'})
xlim([-0.5 + min(log10Q2)  0.5 + max(log10Q2)])
legend('Facilitation', 'Depression', 'No Plasticity')
hold off

%% Median + IQR

% include Q = 1e-3;
hold on
plot(log10Q, median(absRes_fac, 1), 'Color', org,'LineWidth',2)
plot(log10Q, median(absRes_dep, 1), 'Color', blu,'LineWidth',2)
plot(log10Q, median(absRes_null, 1), 'Color', grn,'LineWidth',2)
xticks(log10Q)
xticklabels({'Inf','-7','-6','-5','-4', '-3'})
xlim([-0.5 + min(log10Q)  0.5 + max(log10Q)])
legend('Facilitation', 'Depression', 'No Plasticity')
hold off

neg_fac = median(absRes_fac, 1) - quantile(absRes_fac, 0.25, 1);
neg_dep = median(absRes_dep, 1) - quantile(absRes_dep, 0.25, 1);
neg_null = median(absRes_null, 1) - quantile(absRes_null, 0.25, 1);
pos_fac = quantile(absRes_fac, 0.75, 1) - median(absRes_fac, 1);
pos_dep = quantile(absRes_dep, 0.75, 1) - median(absRes_dep, 1);
pos_null = quantile(absRes_null, 0.75, 1) - median(absRes_null, 1);

hold on
errorbar(log10Q, median(absRes_fac, 1), neg_fac, pos_fac, 'Color', org,'LineWidth',2)
errorbar(log10Q+0.05, median(absRes_dep, 1), neg_dep, pos_dep, 'Color', blu,'LineWidth',2)
errorbar(log10Q+0.1, median(absRes_null, 1), neg_null, pos_null, 'Color', grn,'LineWidth',2)
xticks(log10Q)
xticklabels({'Inf','-7','-6','-5','-4', '-3'})
xlim([-0.5 + min(log10Q)  0.5 + max(log10Q)])
legend('Facilitation', 'Depression', 'No Plasticity')
hold off


% exclude Q = 1e-3
hold on
plot(log10Q2, median(absRes2_fac, 1), 'Color', org,'LineWidth',2)
plot(log10Q2, median(absRes2_dep, 1), 'Color', blu,'LineWidth',2)
plot(log10Q2, median(absRes2_null, 1), 'Color', grn,'LineWidth',2)
xticks(log10Q2)
xticklabels({'Inf','-7','-6','-5','-4'})
xlim([-0.5 + min(log10Q2)  0.5 + max(log10Q2)])
legend('Facilitation', 'Depression', 'No Plasticity')
hold off

neg2_fac = median(absRes2_fac, 1) - quantile(absRes2_fac, 0.25, 1);
neg2_dep = median(absRes2_dep, 1) - quantile(absRes2_dep, 0.25, 1);
neg2_null = median(absRes2_null, 1) - quantile(absRes2_null, 0.25, 1);
pos2_fac = quantile(absRes2_fac, 0.75, 1) - median(absRes2_fac, 1);
pos2_dep = quantile(absRes2_dep, 0.75, 1) - median(absRes2_dep, 1);
pos2_null = quantile(absRes2_null, 0.75, 1) - median(absRes2_null, 1);

hold on
errorbar(log10Q2, median(absRes2_fac, 1), neg2_fac, pos2_fac, 'Color', org,'LineWidth',2)
errorbar(log10Q2+0.05, median(absRes2_dep, 1), neg2_dep, pos2_dep, 'Color', blu,'LineWidth',2)
errorbar(log10Q2+0.1, median(absRes2_null, 1), neg2_null, pos2_null, 'Color', grn,'LineWidth',2)
xticks(log10Q2)
xticklabels({'Inf','-7','-6','-5','-4'})
xlim([-0.5 + min(log10Q2)  0.5 + max(log10Q2)])
legend('Facilitation', 'Depression', 'No Plasticity')
hold off

