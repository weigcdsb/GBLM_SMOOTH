
%% pre-fitting plots
plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\stdp\informal';
cd(plotFolder)

Tpre = data.pre_spk_times;
Tpost = data.post_spk_times(:, 1);
data.vecN = length(data.pre_spk_vec);

% g
g = figure;
cla
plot(linspace(0,length(sim.g)*sim.dt/60,length(sim.g)),sim.g,'k')
box off
saveas(g, '1_g.png')

% overall crosscorrelogram
[d,~] = corr_fast_v3(Tpre, Tpost,-.02,.02,102);
tvec = linspace(-0.02,0.02,102);
tvec = tvec+mean(diff(tvec))/2;

crosscorr = figure;
hold on
bar(tvec(1:end-1)*1e3,d(1:end-1),1,'k','EdgeColor','none');
xlim([-.01 .02]*1e3);
box off
hold off
saveas(crosscorr, '2_crosscorr.png')

% efficacy
isi = [Inf; diff(Tpre)];
quantiles = prctile(isi,linspace(0,100,25));
qx = quantiles(1:end-1)+diff(quantiles)/2;
efficacy = zeros(1, length(quantiles)-1);

for q=1:length(quantiles)-1
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    
    dbase = corr_fast_v3(qspk,Tpost,-.008,-0.004,20);
    d = corr_fast_v3(qspk,Tpost,0.004,.008,20);
    efficacy(q) = (sum(d) - sum(dbase))/length(qspk);
end

effi = figure;
clf
hold on
plot(qx*1000, efficacy, 'o', 'Color', [0, 0.4470, 0.7410],...
    'LineWidth', 2, 'markerfacecolor', [0, 0.4470, 0.7410])
box off
hold off
saveas(effi, '3_effi.png')

% firing rate
fr = figure;
data.vecN = length(data.pre_spk_vec);
hold on
plot(linspace(0,sim.T,data.vecN), filter(ones(2000,1),1,data.pre_spk_vec)/2, 'b', 'LineWidth', 1.5);
plot(linspace(0,sim.T,data.vecN), filter(ones(2000,1),1,data.post_spk_vec)/2, 'k', 'LineWidth', 1.5);
hold off
% ylim([0 50])
xlim([0 sim.T])
box off
saveas(effi, '4_fr.png')

%% fitting results
lam = exp(fit.beta0 + fit.wt_long.*fit.wt_short.*fit.Xc +...
    fit.hist*fit.hist_beta)*data.dt;

plot(1 + fit.stp_basis'*fit.wt_short_param, 'r', 'LineWidth', 2)
plot(fit.beta0)
plot(fit.wt_long)
plot(fit.wt_short)
plot(lam)

fit.synParams.syn_params(1)
fit.synParams.syn_params(2)

%% plots
plotFolder = 'C:\Users\gaw19004\Documents\GitHub\GBLM_SMOOTH\plot\stdp\formal';
cd(plotFolder)

% firing rate
firRate = figure;
hold on
plot(linspace(0,sim.T,data.vecN), filter(ones(2000,1),1,data.pre_spk_vec)/2, 'b', 'LineWidth', 1.5);
plot(linspace(0,sim.T,data.vecN), filter(ones(2000,1),1,data.post_spk_vec)/2, 'k', 'LineWidth', 1.5);
plot(linspace(0,sim.T,data.vecN), filter(ones(2000,1),1,lam)/2, 'r', 'LineWidth', 1.5);
hold off
xlim([0 sim.T])
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off

set(firRate,'PaperUnits','inches','PaperPosition',[0 0 5 3])
saveas(firRate, '1_firRate.svg')
saveas(firRate, '1_firRate.png')

% overall crosscorrelogram
[d,~] = corr_fast_v3(Tpre, Tpost,-.02,.02,102);
[d_fit, lag_fit] = xcorr(data.pre_spk_vec, lam, 20);

tvec = linspace(-0.02,0.02,102);
tvec = tvec+mean(diff(tvec))/2;

corrOverall = figure;
hold on
bar(tvec(1:end-1)*1e3,d(1:end-1),1,'k','EdgeColor','none');
plot(-lag_fit*data.dt*1e3, d_fit*mean(diff(tvec))/data.dt, 'r', 'LineWidth',3)
xlim([-.01 .02]*1e3);
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(corrOverall,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(corrOverall, '2_corrOverall.svg')
saveas(corrOverall, '2_corrOverall.png')

% efficacy
isi = [Inf; diff(Tpre)];
quantiles = prctile(isi,linspace(0,100,25));
qx = quantiles(1:end-1)+diff(quantiles)/2;

efficacy = zeros(1, length(quantiles)-1);
model_efficacy = zeros(1, length(quantiles)-1);
syn_kern = fit.syn(linspace(0, 1, 1/data.dt));
sk = exp(syn_kern);
sk=sk(sk>1);

wt = (1+fit.stp_X*fit.wt_short_param).*fit.wt_long;
for q=1:length(quantiles)-1
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    
    dbase = corr_fast_v3(qspk,Tpost,-.008,-0.004,20);
    d = corr_fast_v3(qspk,Tpost,0.004,.008,20);
    
    n = zeros(1,data.T/data.dt);
    n(ceil(qspk/data.dt))=1;
    model_efficacy(q) = mean(exp(wt(n>0)));
    
    efficacy(q) = (sum(d) - sum(dbase))/length(qspk);
end

effPlot = figure;
clf
hold on
plot(qx*1000, efficacy, 'o', 'Color', [0, 0.4470, 0.7410],...
    'LineWidth', 2, 'markerfacecolor', [0, 0.4470, 0.7410])
plot(qx*1000, (model_efficacy)/nanmean(model_efficacy)*nanmean(efficacy),...
    'Color', [0, 0.4470, 0.7410], 'LineWidth', 3)
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(effPlot,'PaperUnits','inches','PaperPosition',[0 0 4 3])
saveas(effPlot, '3_effPlot.svg')
saveas(effPlot, '3_effPlot.png')

% split cross-correlogram (STP)
cd(strcat(plotFolder, '\corrISI'))
maxd = -Inf;
isi = [Inf; diff(Tpre)];
quantiles = prctile(isi,linspace(0,100,5));

for q=1:length(quantiles)-1
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    if maxd<max(d),maxd=max(d); end
end

for q=1:length(quantiles)-1
    
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    qvec = zeros(1,data.vecN);
    qvec(round(qspk/data.dt))=1;
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    tvec = linspace(-0.025,0.025,64);
    tvec = tvec+mean(diff(tvec))/2;
    [dfit,lag_fit] = xcorr(qvec, lam, 25);
    
    corrISI = figure;
    bar(tvec(1:end-1),d(1:end-1),1,'k','EdgeColor','none');
    hold on
    plot(-lag_fit*data.dt, dfit*mean(diff(tvec))/data.dt, 'r', 'LineWidth',3)
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    hold off
    %     ylim([0 ceil(maxd/50)*50])
    ylim([0 200])
    xlim([-0.01 0.02])
    
    set(corrISI,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    saveas(corrISI, strcat('corrISI_Q', string(q), '.svg'))
    saveas(corrISI, strcat('corrISI_Q', string(q), '.png'))
end

% split cross-correlogram (LTP)
cd(strcat(plotFolder, '\corrT'))
maxd = -Inf;

for q=1:4
    qspk = Tpre(Tpre >= (q-1)*0.25*sim.T & Tpre < q*0.25*sim.T);
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    if maxd<max(d),maxd=max(d); end
end

for q=1:4
    qspk = Tpre(Tpre >= (q-1)*0.25*sim.T & Tpre < q*0.25*sim.T);
    qvec = zeros(1,data.vecN);
    qvec(round(qspk/data.dt))=1;
    d = corr_fast_v3(qspk,Tpost,-.025,.025,64);
    tvec = linspace(-0.025,0.025,64);
    tvec = tvec+mean(diff(tvec))/2;
    [dfit,lag_fit] = xcorr(qvec, lam, 25);
    
    corrT = figure;
    bar(tvec(1:end-1),d(1:end-1),1,'k','EdgeColor','none');
    hold on
    plot(-lag_fit*data.dt, dfit*mean(diff(tvec))/data.dt, 'r', 'LineWidth',3)
    set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
    box off
    hold off
    ylim([0 ceil(maxd/50)*50])
    xlim([-0.01 0.02])
    
    set(corrT,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    saveas(corrT, strcat('corrT_Q', string(q), '.svg'))
    saveas(corrT, strcat('corrT_Q', string(q), '.png'))
end





