
realData.dt = 1e-3;
realData.pre_spk_times = Tlist{pre(1)};
realData.post_spk_times = Tlist{post};

realData.vecN = round(max([realData.pre_spk_times; realData.post_spk_times])/realData.dt);
realData.pre_spk_vec = zeros(1,realData.vecN);
realData.pre_spk_vec(round(realData.pre_spk_times/realData.dt))=1;
realData.post_spk_vec = zeros(1,realData.vecN);
realData.post_spk_vec(round(realData.post_spk_times/realData.dt))=1;



[~, ~, Qopt1] = ...
    tune_smooth_gblm_1d_grad(realData.pre_spk_vec, realData.post_spk_vec,...
    'iter',15, 'hist_tau', 0.01, 'hist_beta', 0, 'doFit', false);

[fit1,~] = smooth_gblm(realData.pre_spk_vec, realData.post_spk_vec,...
    'iter',15, 'Q', Qopt1,...
    'hist_tau', 1e-3, 'hist_beta', -2, 'synParams', fit1.synParams);



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




d = corr_fast_v3(realData.pre_spk_times,realData.post_spk_times(:,1),-.025,.025,64);
t = linspace(-.025,.025,64);
t = t+mean(diff(t))/2;
bar(t,d,1,'EdgeColor','none')
box off; set(gca,'TickDir','out')

% [fit1.syn,fit1.deltaT,fit1.synParams,~]= ...
%     synapse_xcorr_conv({realData.pre_spk_times,realData.post_spk_times},...
%     fit1.syn_hyper_params);
fit1 = rmfield(fit1,'synParams');
fit1 = synConEst(realData,fit1);

disp(fit1.synParams.syn_params(1))
disp(fit1.synParams.syn_params(2))

fit1.synParams.syn_params(1)
fit1.synParams.syn_params(2) = 0.01;

fit1.syn = @(ts) max(0,ts-fit1.synParams.syn_params(1))/...
    fit1.synParams.syn_params(2).*exp(1-max(0,ts-fit1.synParams.syn_params(1))/...
    fit1.synParams.syn_params(2));

edges = linspace(-.02,.02,101);
n = histc(fit1.deltaT(:,1),edges);
bar(edges+mean(diff(edges))/2,n,1,'EdgeColor','none');
hold on
x0 = linspace(-.02,.02,1001);
plot(x0,exp(fit1.syn(x0))*mean(n), 'r')
hold off
box off

tau = 1e-3;
x0 = linspace(0,1,1/realData.dt);
kern_hist = [0 exp(-x0/tau)];
plot(kern_hist)
hist = filter(kern_hist,1,realData.post_spk_vec');
plot(hist(1:10000));
hold on
plot(fit1.hist(1:10000))
hold off

fit1.hist_tau
