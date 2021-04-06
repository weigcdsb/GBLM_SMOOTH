addpath(genpath('C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\helper'));
addpath(genpath('C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\data'));
addpath(genpath('C:\Users\gaw19004\Documents\GitHub\SMOOTH_GBLM_revision'));

load('crcns_ssc-3_dataset23.mat');

%% data
rng(1)
realData.dt = 1e-3;

% facilitation
realData.pre_spk_times = Tlist{22};
realData.post_spk_times = Tlist{281};

% depression
% realData.pre_spk_times = Tlist{136};
% realData.post_spk_times = Tlist{75};


% discard first few minitues?
lb = 0*60;
realData.pre_spk_times = realData.pre_spk_times(realData.pre_spk_times >= lb) - lb;
realData.post_spk_times = realData.post_spk_times(realData.post_spk_times >= lb) - lb;
realData.T = max([realData.pre_spk_times; realData.post_spk_times]);

realData.vecN = round(realData.T/realData.dt);
realData.pre_spk_vec = zeros(1,realData.vecN);
realData.pre_spk_vec(round(realData.pre_spk_times/realData.dt))=1;
realData.post_spk_vec = zeros(1,realData.vecN);
realData.post_spk_vec(round(realData.post_spk_times/realData.dt))=1;

%% check efficacy
syn_hyper_params.coupling_timescale = 0.05;
syn_hyper_params.bin_width = (0.05)*realData.dt;
syn_hyper_params.baseline_nsplines = 4;
[syn,deltaT,synParams,~]= ...
    synapse_xcorr_conv({realData.pre_spk_times,realData.post_spk_times},...
    syn_hyper_params);

Tpre = realData.pre_spk_times;
Tpost = realData.post_spk_times;

isi = [Inf; diff(Tpre)];
quantiles = prctile(isi,linspace(0,100,25));
qx = quantiles(1:end-1)+diff(quantiles)/2;
efficacy = zeros(1, length(quantiles)-1);
x0 = linspace(-.02,.02,1001);
effWindL_pre = range(x0(syn(x0) > 0.01)); % window length to look at efficacy

for q=1:length(quantiles)-1
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    
    d1base = corr_fast_v3(qspk,Tpost,-effWindL_pre-0.005, -0.005,20);
    d1 = corr_fast_v3(qspk,Tpost,min(x0(syn(x0) > 0.01)),...
        max(x0(syn(x0) > 0.01)),20);
    
    n1 = zeros(1,round(realData.T/realData.dt));
    n1(ceil(qspk/realData.dt))=1;
    
    efficacy(q) = (sum(d1) - sum(d1base))/length(qspk);
end

% facilitation
stp_Nq = 6; stp_Nm = 400; stp_Ns = 10;

% depression
% stp_Nq = 6; stp_Nm = 150; stp_Ns = 10;

stp_basis = getBasis('rcos', stp_Nq, stp_Nm, stp_Ns,0);
hold on
plot(qx*1000, efficacy, 'o', 'Color', [0, 0.4470, 0.7410],...
    'LineWidth', 2, 'markerfacecolor', [0, 0.4470, 0.7410])
plot(stp_basis'*max(abs(efficacy)))
hold off

%% fit data

% select Q: use first 5min
[~, ~, Q] = tune_smooth_gblm_2d_grad2(realData.pre_spk_vec(1:(5*60/realData.dt)),...
    realData.post_spk_vec(1:(5*60/realData.dt)),...
    'hist_tau', 0.01, 'Q0', [1e-5 1e-5],...
    'doFit', false,...
    'QUB', [3*1e-3 1e-3],...
    'stp_Nq', stp_Nq, 'stp_Nm', stp_Nm, 'stp_Ns', stp_Ns,...
    'synParams', synParams,'stp_tau',0.1);

% full model
[fit, ~] = smooth_gblm2(realData.pre_spk_vec, realData.post_spk_vec,...
    'iter',30, 'hist_tau', 0.01, 'Q', diag(Q),...
    'stp_Nq', stp_Nq, 'stp_Nm', stp_Nm, 'stp_Ns',...
    stp_Ns, 'synParams', synParams,'stp_tau',0.1);

lam_full = exp(fit.beta0 +...
    fit.wt_long.*fit.wt_short.*fit.Xc +...
    fit.hist*fit.hist_beta)*realData.dt;

% null model
lam_null = mean(realData.post_spk_vec);


% sub-model 1: wt_long = constant, wt_short = 1
fit1 = noSTP(realData, fit.hist_tau, fit.hist_beta,...
    fit.stp_Nq, fit.stp_Nm, fit.stp_Ns,...
    diag([fit.Q(1, 1) 0]), fit.synParams);
lam1 = exp(fit1.beta0 +...
    fit1.wt_long.*fit.Xc +...
    fit.hist*fit.hist_beta)*realData.dt;


% sub-model 2: wt_long = time-varying, wt_short = 1
fit2 = noSTP(realData, fit.hist_tau, fit.hist_beta,...
    fit.stp_Nq, fit.stp_Nm, fit.stp_Ns,...
    fit.Q, fit.synParams);
lam2 = exp(fit2.beta0 +...
    fit2.wt_long.*fit.Xc +...
    fit.hist*fit.hist_beta)*realData.dt;

% sub-model 3: wt_long = constant, wt_short = on
fit3 = smooth_gblm(realData.pre_spk_vec, realData.post_spk_vec,...
    'iter',30, 'hist_tau', 0.01, 'hist_beta', fit.hist_beta,...
    'Q', diag([fit.Q(1, 1) 0]),...
    'stp_Nq', stp_Nq, 'stp_Nm', stp_Nm, 'stp_Ns', stp_Ns, 'synParams', fit.synParams,...
    'stp_tau',0.1);
lam3 = exp(fit3.beta0 +...
    fit3.wt_long.*fit3.wt_short.*fit.Xc +...
    fit.hist*fit.hist_beta)*realData.dt;

save('facilitation.mat');
% save('depression.mat');
%% load data

% load('facilitation.mat')
% load('depression.mat');

%% likelihood

% llhd & log2-likelihood: drop the constant term
llhd = @(lam) sum(-lam + log((lam+(lam==0))).*(realData.post_spk_vec'));

[llhd(lam_full) llhd(lam_null)]/realData.T
[llhd(lam_full) llhd(lam_null)...
    llhd(lam1) llhd(lam2) llhd(lam3)]'/realData.T

% transmission
postIdx = @(delta) min(find(realData.pre_spk_vec == 1) + delta, realData.vecN);
l2hd_trans = @(lam, post_spk) sum(log2(exp(-lam)) + log2((lam+(lam==0))).*(post_spk'));
idx = [];
for j = 1:10; idx = [idx postIdx(j)]; end

[l2hd_trans(lam_full(idx), realData.post_spk_vec(idx)),...
    l2hd_trans(lam_null, realData.post_spk_vec(idx))]/sum(realData.pre_spk_vec == 1)

[l2hd_trans(lam_full(idx), realData.post_spk_vec(idx))-l2hd_trans(lam_null, realData.post_spk_vec(idx))]/sum(realData.pre_spk_vec == 1)


[l2hd_trans(lam_full(idx), realData.post_spk_vec(idx)),...
    l2hd_trans(lam_null, realData.post_spk_vec(idx)),...
    l2hd_trans(lam1(idx), realData.post_spk_vec(idx)),...
    l2hd_trans(lam2(idx), realData.post_spk_vec(idx)),...
    l2hd_trans(lam3(idx), realData.post_spk_vec(idx))]'/sum(realData.pre_spk_vec == 1)



%% parameter trace
% plotFolder = 'C:\Users\gaw19004\Documents\GitHub\SMOOTH_GBLM_revision\plots\depression';
plotFolder = 'C:\Users\gaw19004\Documents\GitHub\SMOOTH_GBLM_revision\plots\facilitation';

cd(plotFolder)

lam = exp(fit.beta0 + fit.wt_long.*fit.wt_short.*fit.Xc + fit.hist*fit.hist_beta)*realData.dt;
% mse = mean((realData.post_spk_vec' - lam).^2);
% mse
idx = 1:size(fit.beta0);

% baseline: useless
baseLine = figure;
hold on
plot(idx, fit.beta0, 'r', 'LineWidth', 3)
plot(idx, fit.beta0 + sqrt(squeeze(fit.W(1, 1, :))), 'r:', 'LineWidth', 2)
plot(idx, fit.beta0 - sqrt(squeeze(fit.W(1, 1, :))), 'r:', 'LineWidth', 2)
ylim([min(fit.beta0)-1 max(fit.beta0)+1])
xlim([0 (60*60)/realData.dt])
xticks((0:5:60)*60*1e3)
xticklabels(arrayfun(@num2str, 0:5:60, 'UniformOutput', 0))
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(baseLine,'PaperUnits','inches','PaperPosition',[0 0 10 3])
saveas(baseLine, '0_baseLine.svg')
saveas(baseLine, '0_baseLine.png')

% ltp
ltp = figure;
idx = 1:size(fit.beta0);
hold on
plot(idx, fit.wt_long, 'r', 'LineWidth', 3)
plot(idx, fit.wt_long + sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
plot(idx, fit.wt_long - sqrt(squeeze(fit.W(2, 2, :))), 'r:', 'LineWidth', 2)
ylim([0 4])
xlim([0 (60*60)/realData.dt])
xticks((0:5:60)*60*1e3)
xticklabels(arrayfun(@num2str, 0:5:60, 'UniformOutput', 0))
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(ltp,'PaperUnits','inches','PaperPosition',[0 0 10 3])
saveas(ltp, '1_ltp.svg')
saveas(ltp, '1_ltp.png')

% stp
se_wt_short = zeros(size(fit.stp_X, 1), 1);
for k = 1:size(fit.stp_X, 1)
    se_wt_short(k) = sqrt(fit.stp_X(k,:)*fit.covB*(fit.stp_X(k,:))');
end

stp = figure;
hold on
plot(idx, fit.wt_short, 'r', 'LineWidth', 2)
plot(idx, fit.wt_short + se_wt_short, 'r:', 'LineWidth', 2)
plot(idx, fit.wt_short - se_wt_short, 'r:', 'LineWidth', 2)
hold off
ylim([-1, 2])
xlim([0 (60*60)/realData.dt])
xticks((0:5:60)*60*1e3)
xticklabels(arrayfun(@num2str, 0:5:60, 'UniformOutput', 0))
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(stp,'PaperUnits','inches','PaperPosition',[0 0 10 3])
saveas(stp, '2_stp.svg')
saveas(stp, '2_stp.png')



% modification function
se_mod_fn = sqrt(diag(fit.stp_basis'*fit.covB*fit.stp_basis));
modFun = figure;
semilogx(1 + [fit.stp_basis'*fit.wt_short_param;zeros(1000, 1)], 'r', 'LineWidth', 3)
hold on
plot(1 + fit.stp_basis'*fit.wt_short_param + se_mod_fn,'r:', 'LineWidth', 2)
plot(1 + fit.stp_basis'*fit.wt_short_param - se_mod_fn,'r:', 'LineWidth', 2)
hold off
yline(1, 'k--', 'LineWidth', 2);
% line(xlim(),[1 1])
% ep = 300; % depression
ep = 800; % facilitation
xlim([2 ep])
% tickMark = [2 5 10 25 50 100 200]; % depression
tickMark = [2 5 10 25 50 100 200 400 800]; % facilitation
xticks(tickMark)
xticklabels(arrayfun(@num2str, tickMark, 'UniformOutput', 0))
ylim([0.7, 1.4])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(modFun,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(modFun, '3_modFun.svg')
saveas(modFun, '3_modFun.png')

% histogram for ISI
isiHis = figure;
histogram(log10(diff(realData.pre_spk_times/realData.dt)),'EdgeColor','none')
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
xlim([log10(2) log10(ep)])
xticks(log10(tickMark))
xticklabels(arrayfun(@num2str, tickMark, 'UniformOutput', 0))
ylim([0 2000])
box off

set(isiHis,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(isiHis, '4_isiHis.svg')
saveas(isiHis, '4_isiHis.png')

%% overall cross-correlogram
[d,~] = corr_fast_v3(Tpre, Tpost,-.02,.02,102);
[d_fit, lag_fit] = xcorr(realData.pre_spk_vec, lam, 20);

tvec = linspace(-0.02,0.02,102);
tvec = tvec+mean(diff(tvec))/2;

corrOverall = figure;
hold on
bar(tvec(1:end-1)*1e3,d(1:end-1),1,'k','EdgeColor','none');
plot(-lag_fit, d_fit*mean(diff(tvec))/realData.dt, 'r', 'LineWidth',3)
ylim([0 700])
xlim([-.01 .02]*1e3);
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(corrOverall,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(corrOverall, '5_corrOverall.svg')
saveas(corrOverall, '5_corrOverall.png')

%% "contribution", analogue of efficacy
% (sum(xcorr_syn)-sum(baseline))/n_post
effWindL = range(x0(fit.syn(x0) > 0.01)); % window length to look at efficacy
% dbase = corr_fast_v3(Tpre,Tpost,-effWindL-0.005, -0.005,20);
dbase = corr_fast_v3(Tpre,Tpost,-max(x0(syn(x0) > 0.01)),-min(x0(syn(x0) > 0.01)),20);
d = corr_fast_v3(Tpre,Tpost,min(x0(syn(x0) > 0.01)),...
    max(x0(syn(x0) > 0.01)),20);

(sum(d) - sum(dbase))/length(Tpre)
(sum(d) - sum(dbase))/length(Tpost)
%% cross-correlogram: LTP
cd(strcat(plotFolder, '\splitT'))
maxd = -Inf;
QT = 4;

for q=1:QT
    qspk = Tpre(Tpre >= (q-1)*(1/QT)*realData.T & Tpre < q*(1/QT)*realData.T);
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    if maxd<max(d),maxd=max(d); end
end

for q=1:QT
    qspk = Tpre(Tpre >= (q-1)*(1/QT)*realData.T & Tpre < q*(1/QT)*realData.T);
    qvec = zeros(1,realData.vecN);
    qvec(round(qspk/realData.dt))=1;
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    tvec = linspace(-0.025,0.025,64);
    tvec = tvec+mean(diff(tvec))/2;
    [dfit,lag_fit] = xcorr(qvec, lam, 25);
    
    corrT = figure;
    bar(tvec(1:end-1),d(1:end-1),1,'k','EdgeColor','none');
    hold on
    plot(-lag_fit*realData.dt, dfit*mean(diff(tvec))/realData.dt, 'r', 'LineWidth',3)
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    hold off
%     ylim([0 ceil(maxd/50)*50])
    ylim([0 400])
    xlim([-0.01 0.02])
    
    set(corrT,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    saveas(corrT, strcat('corrT_Q', string(q), '.svg'))
    saveas(corrT, strcat('corrT_Q', string(q), '.png'))
end


%% cross-correlogram: STP
cd(strcat(plotFolder, '\splitISI'))
maxd = -Inf;
isi = [Inf; diff(Tpre)];
quantiles = prctile(isi,linspace(0,100,5));

for q=1:length(quantiles)-1
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    d = corr_fast_v3(qspk,Tpost,-.025,.025,50);
    if maxd<max(d),maxd=max(d); end
end

for q=1:length(quantiles)-1
    
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    qvec = zeros(1,realData.vecN);
    qvec(round(qspk/realData.dt))=1;
    d = corr_fast_v3(qspk,Tpost,-.025,.025,50);
    tvec = linspace(-0.025,0.025,50);
    tvec = tvec+mean(diff(tvec))/2;
    [dfit,lag_fit] = xcorr(qvec, lam, 25);
    
    corrISI = figure;
    bar(tvec(1:end-1),d(1:end-1),1,'k','EdgeColor','none');
    hold on
    plot(-lag_fit*realData.dt, dfit*mean(diff(tvec))/realData.dt, 'r', 'LineWidth',3)
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    hold off
%     ylim([0 ceil(maxd/50)*50])
    ylim([0 600])
    xlim([-0.01 0.02])
    
    set(corrISI,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    saveas(corrISI, strcat('corrISI_Q', string(q), '.svg'))
    saveas(corrISI, strcat('corrISI_Q', string(q), '.png'))
end



%% efficacy trace 1: T
cd(plotFolder)
% bootstrap to spikes may not be appropriate for efficacy over T
% but should be OK over ISI.

windowL = 5*60;
stepL = 1*60;
bootBlock = 10; % number of blocks
innerWindowL = windowL/5; % fix the window length (sampling)
nSub = 10;
nBoot = 1000;

starts = 0:stepL:(realData.T - windowL);
efficacyT = zeros(1, length(starts));
effTInnerSD = zeros(1, length(starts));
effTInnerMean = zeros(1, length(starts));
effBootMean = zeros(1, length(starts));
effTBootSD = zeros(1, length(starts));
model_efficacyT = zeros(1, length(starts));

wt = (1 + fit.stp_X*fit.wt_short_param).*fit.wt_long;
for s = 1:length(starts)
    Tpre_trunc = Tpre(Tpre >= starts(s) & Tpre <= starts(s) + windowL);
    qvec = zeros(1,realData.vecN);
    qvec(round(Tpre_trunc/realData.dt))=1;
    
    [dbase,dbaseM] = corr_fast_v3(Tpre_trunc,Tpost,-effWindL-0.005, -0.005,20);
    [d,dM] = corr_fast_v3(Tpre_trunc,Tpost,min(x0(fit.syn(x0) > 0.01)),...
        max(x0(fit.syn(x0) > 0.01)),20);
    efficacyT(s) = (sum(d) - sum(dbase))/length(Tpre_trunc);
    
    % iid bootstrap: spiking bootstap
    base_set = histc(dbaseM(:,2),1:length(Tpre_trunc));
    post_set = histc(dM(:,2),1:length(Tpre_trunc));
    b = bootstrp(1000,'sum',[post_set-base_set]);
    effTInnerMean(s) = mean(b/length(Tpre_trunc));
    effTInnerSD(s) = std(b/length(Tpre_trunc));
    
    % block bootstrap: 10-block
    effBoot = zeros(1, nBoot);
    idx = (starts(s)/realData.dt+1):((starts(s) + windowL)/realData.dt);
    for b = 1:nBoot
        blockSample = datasample(1:bootBlock, bootBlock);
        Tpre_trunc_spk2 = [];
        Tpost_trunc_spk2 = [];
        for bl = 1:bootBlock
            idxBlock = idx(((blockSample(bl)-1)*round(length(idx)/bootBlock)+ 1):...
                (blockSample(bl)*round(length(idx)/bootBlock)));
            Tpre_trunc_spk2 = [Tpre_trunc_spk2, realData.pre_spk_vec(idxBlock)];
            Tpost_trunc_spk2 = [Tpost_trunc_spk2, realData.post_spk_vec(idxBlock)];
        end
        Tpre_trunc2 = find(Tpre_trunc_spk2 > 0)*realData.dt;
        Tpost_trunc2 = find(Tpost_trunc_spk2 > 0)*realData.dt;
        dbaseBoot = corr_fast_v3(Tpre_trunc2,Tpost_trunc2,-effWindL-0.005, -0.005,20);
        dBoot = corr_fast_v3(Tpre_trunc2,Tpost_trunc2,min(x0(fit.syn(x0) > 0.01)),...
            max(x0(fit.syn(x0) > 0.01)),20);
        effBoot(b) = (sum(dBoot) - sum(dbaseBoot))/length(Tpre_trunc2);
    end
    effTBootSD(s) = std(effBoot);
    effBootMean(s) = mean(effBoot);
    
    % model efficacy: use correlation
    [dfit,lag_fit] = xcorr(qvec, lam, 20);
    model_efficacyT(s) = ...
        sum( dfit(-lag_fit*realData.dt >= min(x0(fit.syn(x0) > 0.01)) &...
        -lag_fit*realData.dt <= max(x0(fit.syn(x0) > 0.01)))) -...
        sum(dfit(-lag_fit*realData.dt >= -effWindL-0.005 &...
        -lag_fit*realData.dt <= -0.005));
    model_efficacyT(s) = model_efficacyT(s)/length(Tpre_trunc);
    
    
end

model_efficacyT = (model_efficacyT)/nanmean(model_efficacyT)*nanmean(efficacyT);

effPlot_bootSpk = figure;
hold on
plot(starts + windowL/2, efficacyT,...
    'b', 'LineWidth', 2);
plot(starts + windowL/2, efficacyT + effTInnerSD,...
    'b:', 'LineWidth', 2)
plot(starts + windowL/2, efficacyT - effTInnerSD,...
    'b:', 'LineWidth', 2)
plot(starts + windowL/2, model_efficacyT,...
    'Color', 'r', 'LineWidth', 3);
ylim([-0.05 0.2])
xlim([0 60*60])
xticks((0:5:60)*60)
xticklabels(arrayfun(@num2str, 0:5:60, 'UniformOutput', 0))
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(effPlot_bootSpk,'PaperUnits','inches','PaperPosition',[0 0 10 3])
saveas(effPlot_bootSpk, '6_effPlot_bootSpk.svg')
saveas(effPlot_bootSpk, '6_effPlot_bootSpk.png')

effPlot_bootT = figure;
hold on
plot(starts + windowL/2, efficacyT,...
    'b', 'LineWidth', 2);
plot(starts + windowL/2, efficacyT + effTBootSD,...
    'b:', 'LineWidth', 2)
plot(starts + windowL/2, efficacyT - effTBootSD,...
    'b:', 'LineWidth', 2)
plot(starts + windowL/2, model_efficacyT,...
    'Color', 'r', 'LineWidth', 3);
ylim([-0.05 0.2])
xlim([0 60*60])
xticks((0:5:60)*60)
xticklabels(arrayfun(@num2str, 0:5:60, 'UniformOutput', 0))
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(effPlot_bootT,'PaperUnits','inches','PaperPosition',[0 0 10 3])
saveas(effPlot_bootT, '7_effPlot_bootT.svg')
saveas(effPlot_bootT, '7_effPlot_bootT.png')

%% efficacy trace 2: ISI
isi = [Inf; diff(Tpre)];
quantiles = prctile(isi,linspace(0,100,25));
qx = quantiles(1:end-1)+diff(quantiles)/2;
efficacy = zeros(1, length(quantiles)-1);
model_efficacy_isi = zeros(1, length(quantiles)-1);
model_efficacy2 = zeros(1, length(quantiles)-1);

for q=1:length(quantiles)-1
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    qvec = zeros(1,realData.vecN);
    qvec(round(qspk/realData.dt))=1;
    
    [d1base,dbaseM] = corr_fast_v3(qspk,Tpost,-effWindL-0.005, -0.005,20);
    [d1,dM] = corr_fast_v3(qspk,Tpost,min(x0(fit.syn(x0) > 0.01)),...
        max(x0(fit.syn(x0) > 0.01)),20);
    efficacy(q) = (sum(d1) - sum(d1base))/length(qspk);
    
    
    base_set = histc(dbaseM(:,2),1:length(qspk));
    post_set = histc(dM(:,2),1:length(qspk));
    b = bootstrp(1000,'sum',[post_set-base_set]);
    bootMean(q) = mean(b/length(qspk));
    bootSD(q) = std(b/length(qspk));
    
    % model efficacy 1: use correlation
    [dfit,lag_fit] = xcorr(qvec, lam, 20);
    model_efficacy_isi(q) = ...
        sum( dfit(-lag_fit*realData.dt >= min(x0(fit.syn(x0) > 0.01)) &...
        -lag_fit*realData.dt <= max(x0(fit.syn(x0) > 0.01)))) -...
        sum(dfit(-lag_fit*realData.dt >= -effWindL-0.005 &...
        -lag_fit*realData.dt <= -0.005));
    model_efficacy_isi(q) = model_efficacy_isi(q)/length(qspk);
    
    
end

effPlot_isi = figure;
semilogx(qx*1000, bootMean, 'o', 'Color', 'b',...
    'LineWidth', 3)
hold on
errorbar(qx*1000, bootMean,bootSD,'o','CapSize',0, 'LineWidth', 1.5)
plot(qx*1000, (model_efficacy_isi)/nanmean(model_efficacy_isi)*nanmean(efficacy),...
    'Color', 'r', 'LineWidth', 3)
xlim([2 ep])
xticks(tickMark)
xticklabels(arrayfun(@num2str, tickMark, 'UniformOutput', 0))
ylim([-0.05 0.15])
hold off
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(effPlot_isi,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(effPlot_isi, '8_effPlot_isi.svg')
saveas(effPlot_isi, '8_effPlot_isi.png')

