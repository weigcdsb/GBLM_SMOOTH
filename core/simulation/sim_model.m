function [data,sim] = sim_model(data,sim)

rng(sim.seed);

% Presynaptic Spike Times
if length(sim.pPreSpike)==1 % homogeneous Poisson
    mean_isi = data.dt/sim.pPreSpike;
    data.pre_spk_times = cumsum(exprnd(mean_isi,sim.T/mean_isi*5,1));
    data.pre_spk_times = data.pre_spk_times(data.pre_spk_times<sim.T);
else % inhomogeneous Poisson
    data.pre_spk_times = find(poissrnd(sim.pPreSpike)>0);
    data.pre_spk_times = data.pre_spk_times+rand(size(data.pre_spk_times))-.5;
    data.pre_spk_times = sort(data.pre_spk_times*data.dt);
end
data.pre_spk_vec = zeros(1,sim.vecN);
data.pre_spk_vec(round(data.pre_spk_times/data.dt))=1;

% Synaptic Connection
sim.syn = @(ts) max(0,ts-sim.alpha_dt)/sim.alpha_tau.*exp(1-max(0,ts-sim.alpha_dt)/sim.alpha_tau);
syn_kern = sim.syn(linspace(0, 1, 1/data.dt));
Xc = filter(syn_kern, 1, data.pre_spk_vec');

% Short-Term Plasticity
isi = diff(find(data.pre_spk_vec>0)*data.dt);
sim.stp_basis = getBasis('rcos', sim.stp_Nq, sim.stp_Nm, sim.stp_Ns,0);
Bm = padarray(sim.stp_basis,[0 max(round(isi/data.dt))],'post');
Bm_dt = Bm(:,round(isi/data.dt));
s=zeros(sim.stp_Nq,sim.vecN);
for m =1:sim.stp_Nq
    s(m,data.pre_spk_vec>0)=[0 Bm_dt(m,:)];
end
x0 = linspace(0,1,1/data.dt);
kern_stp = exp(-x0/sim.stp_tau);
sim.stp_X = filter(kern_stp,1,s');
wt_short = 1 + sim.stp_X*sim.stp_B;

% Generate Post-Synaptic Spikes
kern_hist = exp(-x0/sim.hist_tau);
hist_filt = [mean(sim.beta0) kern_hist*sim.hist_beta];
data.post_spk_times = sim_glm(hist_filt,data.dt,sim.T,1,sim.beta0-mean(sim.beta0)+sim.wt_long.*wt_short.*Xc);
data.post_spk_vec = zeros(1,sim.vecN);
data.post_spk_vec(round(data.post_spk_times/data.dt))=1;

kern_hist = [0 kern_hist];
sim.hist = filter(kern_hist,1,data.post_spk_vec');
sim.lam = exp(sim.beta0 + sim.wt_long.*wt_short.*Xc + sim.hist*sim.hist_beta);

end
