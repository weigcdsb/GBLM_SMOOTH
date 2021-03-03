addpath(genpath('C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH'));
load('crcns_ssc-3_dataset23.mat');

%% data 1
realData.dt = 1e-3;
realData.pre_spk_times = Tlist{pre(2)};
realData.post_spk_times = Tlist{post};

realData.vecN = round(max([realData.pre_spk_times; realData.post_spk_times])/realData.dt);
realData.pre_spk_vec = zeros(1,realData.vecN);
realData.pre_spk_vec(round(realData.pre_spk_times/realData.dt))=1;
realData.post_spk_vec = zeros(1,realData.vecN);
realData.post_spk_vec(round(realData.post_spk_times/realData.dt))=1;

[fit1, ~, Qopt1] = ...
    tune_smooth_gblm_1d_grad2(realData.pre_spk_vec, realData.post_spk_vec,...
    'iter',15, 'hist_tau', 0.01);


idx = 1:size(fit1.beta0);

subplot(3, 1, 1)
plot(idx, fit1.beta0,...
    idx, fit1.beta0 + squeeze(fit1.W(1, 1, :)), 'b:',...
    idx, fit1.beta0 - squeeze(fit1.W(1, 1, :)), 'b:');
title('beta_0')

subplot(3, 1, 2)
plot(idx, fit1.wt_long,...
    idx, fit1.wt_long + squeeze(fit1.W(2, 2, :)), 'b:',...
    idx, fit1.wt_long - squeeze(fit1.W(2, 2, :)), 'b:');
title('wt_{long}')

se_mod_fn = sqrt(diag(fit1.stp_basis'*fit1.covB*fit1.stp_basis));

subplot(3, 1, 3)
hold on
plot(1 + fit1.stp_basis'*fit1.wt_short_param, 'r', 'LineWidth', 3)
plot(1 + fit1.stp_basis'*fit1.wt_short_param + se_mod_fn,'r:', 'LineWidth', 2)
plot(1 + fit1.stp_basis'*fit1.wt_short_param - se_mod_fn,'r:', 'LineWidth', 2)
hold off


se_wt_short = zeros(size(fit1.stp_X, 1), 1);
for k = 1:size(fit1.stp_X, 1)
    se_wt_short(k) = sqrt(fit1.stp_X(k,:)*fit1.covB*(fit1.stp_X(k,:))');
end

plot(idx, fit1.wt_short,...
    idx, fit1.wt_short + se_wt_short, 'b:',...
    idx, fit1.wt_short - se_wt_short, 'b:');
title('wt_{short}')

plot(idx, fit1.wt_long.*fit1.wt_short)
save('tmp2_hist_beta.mat');
