
rng(0);

data.T = 300;
data.dt = 0.001;
data.vecN = round(data.T/data.dt);

sim.pPreSpike = 10*data.dt;
sim.pPreSpike = ones(data.vecN,1)*10*data.dt;
sim.pPreSpike(round(data.vecN/3):round(data.vecN*2/3)) = 0;
sim.alpha_dt = 0.002;
sim.alpha_tau = 0.001;
sim.stp_Nq = 5;
sim.stp_Nm = 450;
sim.stp_Ns = 50;
sim.stp_tau= 1; % timeconstant for decay
sim.stp_B = [-1 -2 -3 0 0]'*.0;
sim.hist_tau = .01;
sim.hist_beta = -5;

sim.beta0 = ones(data.vecN,1)*2;
% sim.beta0(round(data.vecN/3):round(data.vecN*2/3)) = -1;
sim.wt_long = ones(data.vecN,1)*2;
% sim.beta0 = 2+detrend(cumsum(randn(data.vecN,1)*10e-4));
% sim.wt_long = 2+detrend(cumsum(randn(data.vecN,1)*10e-3));

[data,sim] = sim_model(data,sim);



fit.stp_X = sim.stp_X;
fit.stp_Nq = size(fit.stp_X,2);
fit.hist = sim.hist;
fit.hist_beta = sim.hist_beta;

fit.syn_hyper_params.coupling_timescale = 0.05;
fit.syn_hyper_params.bin_width = (0.05)*data.dt;
fit.syn_hyper_params.baseline_nsplines = 4; 

[fit.syn,deltaT,fit.synParams,~]= synapse_xcorr({data.pre_spk_times,data.post_spk_times},fit.syn_hyper_params);
fit.syn_kern = fit.syn(linspace(0, .1, .1/data.dt));
fit.Xc = filter(fit.syn_kern, 1, data.pre_spk_vec');

%%
figure(4)
edges = linspace(-.02,.02,101);
n = histc(deltaT(:,1),edges);
bar(edges+mean(diff(edges))/2,n,1,'EdgeColor','none');
hold on
x0 = linspace(-.02,.02,1001);
plot(x0,exp(fit.syn(x0))*mean(n))
hold off

%%

fit.F = eye(2);
fit.Q = [0.000001 0; 0 0.0001];
fit.doOracleSTP = true;
fit.oracle_stp_B = sim.stp_B;
[fit,fit_trace] = fit_model(data, fit,'iter',1);


%%

plot_simfit(sim,fit)

%%

figure(1)
subplot(4,1,1)
plot(filter(ones(2000,1),1,data.pre_spk_vec)/2)
hold on
plot(filter(ones(2000,1),1,data.post_spk_vec)/2)
hold off
ylabel('Firing Rates')
subplot(4,1,2)
plot(sim.beta0)
hold on
plot(fit.beta0)
plot(fit.beta0+sqrt(squeeze(fit.W(1,1,:))),'r:');
plot(fit.beta0-sqrt(squeeze(fit.W(1,1,:))),'r:');
hold off
ylabel('beta0')

subplot(4,1,3)
plot(sim.wt_long)
hold on
plot(fit.wt_long)
plot(fit.wt_long+sqrt(squeeze(fit.W(2,2,:))),'r:');
plot(fit.wt_long-sqrt(squeeze(fit.W(2,2,:))),'r:');
hold off
ylabel('wt-long')

subplot(4,1,4)
plot(squeeze(fit.W(2,2,:)))
ylabel('W-wt-long')