
% rng(0);

data.T = 300;
data.dt = 0.001;
data.vecN = round(data.T/data.dt);

sim.seed = 0;
sim.pPreSpike = 10*data.dt;
sim.alpha_dt = 0.002;
sim.alpha_tau = 0.001;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1; % timeconstant for decay
sim.stp_B = [1 1 1 1 1]'*-.05;
sim.hist_tau = .01;
sim.hist_beta = 0;

% sim.beta0 = ones(data.vecN,1)*2;
% sim.wt_long = ones(data.vecN,1)*2;
sim.beta0 = 2+detrend(cumsum(randn(data.vecN,1)*10e-4));
sim.wt_long = 4+detrend(cumsum(randn(data.vecN,1)*10e-3));

[data,sim] = sim_model(data,sim);



fit.stp_X = sim.stp_X;
fit.stp_Nq = size(fit.stp_X,2);
fit.hist = sim.hist;
fit.hist_beta = sim.hist_beta;

fit.syn_hyper_params.coupling_timescale = 0.05;
fit.syn_hyper_params.bin_width = (0.05)*data.dt;
fit.syn_hyper_params.baseline_nsplines = 4; 

[fit.syn,deltaT,fit.synParams,~]= synapse_xcorr_conv({data.pre_spk_times,data.post_spk_times},fit.syn_hyper_params);
% [fit.syn,deltaT,fit.synParams,~]= synapse_xcorr({data.pre_spk_times,data.post_spk_times},fit.syn_hyper_params);

% use oracle/fixed...
% fit.synParams.syn_params(1) = sim.alpha_dt;
% fit.synParams.syn_params(2) = sim.alpha_tau;
% fit.syn = @(ts) max(0,ts-fit.synParams.syn_params(1))/fit.synParams.syn_params(2).*exp(1-max(0,ts-fit.synParams.syn_params(1))/fit.synParams.syn_params(2));

fit.syn_kern = fit.syn(linspace(0, .1, .1/data.dt));
fit.Xc = filter(fit.syn_kern, 1, data.pre_spk_vec');


%% Verify that alpha function is well-fit and that simulation has a clear monosynaptic connection

figure(4)
subplot(1,3,1)
edges = linspace(-.02,.02,101);
n = histc(deltaT(:,1),edges);
bar(edges+mean(diff(edges))/2,n,1,'EdgeColor','none');
hold on
x0 = linspace(-.02,.02,1001);
plot(x0,exp(fit.syn(x0))*mean(n))
plot(x0,exp(sim.syn(x0))*mean(n))
hold off
box off

subplot(1,3,2)
[aagram,~] = corrFast(data.pre_spk_times,data.pre_spk_times,min(edges),max(edges),length(edges));
[~,i]=max(aagram);
aagram(i)=NaN;
bar(edges+mean(diff(edges))/2,aagram,1,'EdgeColor','none');
box off

subplot(1,3,3)
[apgram,~] = corrFast(data.post_spk_times,data.post_spk_times,min(edges),max(edges),length(edges));
[~,i]=max(apgram);
apgram(i)=NaN;
bar(edges+mean(diff(edges))/2,apgram,1,'EdgeColor','none');
box off

[[sim.alpha_dt; sim.alpha_tau] fit.synParams.syn_params(1:2)]*1000

%% Fit with fixed-Q

fit.F = eye(2);
fit.Q = eye(2)*10e-6;
fit.doFiltOnly = false;
[fit,fit_trace] = fit_model(data, fit,'iter',5);

%% Optimize Q using predictions from filtering
% assumes that p(y|Qb,Qw) ~ p(y|Qb)p(y|Qw)
% should make this section a function eventually

nq = 16;
Qvec = logspace(-9,-3,nq);
qllhd=[]; qllhd_pred=[];
clear fitQ
for q=1:length(Qvec)
    fprintf('Qbeta0 %02i/%02i...',q,length(Qvec))
    fit.Q = [Qvec(q) 0; 0 0];
    fit.doFiltOnly=true;
    [fit,fit_trace] = fit_model(data, fit,'iter',1);
    qllhd(q)=fit.llhd;
    qllhd_pred(q)=fit.llhd_pred;

    fitQ(q)=fit;
end
[~,qi]=max(qllhd_pred);

Qvec = logspace(-9,-3,nq);
qwllhd=[]; qwllhd_pred=[];
clear fitQw
for q=1:length(Qvec)
    fprintf('Qwtlong %02i/%02i...',q,length(Qvec))
    fit.Q = [Qvec(qi) 0; 0 Qvec(q)];
    fit.doFiltOnly=true;
    [fit,fit_trace] = fit_model(data, fit,'iter',1);
    qwllhd(q)=fit.llhd;
    qwllhd_pred(q)=fit.llhd_pred;

    fitQw(q)=fit;
end
[~,qj]=max(qwllhd_pred);

figure(3)
subplot(1,2,1)
semilogx(Qvec,qllhd,Qvec,qllhd_pred)
ylabel('log likelihood')
xlabel('Q')
title('beta0')
% ylim([qllhd_pred(1) max(qllhd)])
legend({'b_{k|k}','b_{k|k-1}'})
subplot(1,2,2)
semilogx(Qvec,qwllhd,Qvec,qwllhd_pred)
xlabel('Q')
title('wtlong')
ylim([2*qwllhd_pred(1) .5*qwllhd_pred(1)])

fit.Q = [Qvec(qi) 0; 0 Qvec(qj)];
fit.doFiltOnly = false;
[fit,fit_trace] = fit_model(data, fit,'iter',10);

%% Show results

plot_simfit(sim,fit)
plot_simfit(sim,fitQ(3))
plot_simfit(sim,fitQw(qj))
plot_simfit(sim,fit_trace(5))

%%

plot(sim.wt_long.*(1+sim.stp_X*sim.stp_B))
hold on
plot(fit.wt_long.*fit.wt_short)
hold off
