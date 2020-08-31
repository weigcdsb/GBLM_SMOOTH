function [data,sim] = sim_model_stdp(sim)

rng(sim.seed)

% Presynaptic Spike Times
if length(sim.pPreSpike)==1 % homogeneous Poisson
    mean_isi = sim.dt/sim.pPreSpike;
    data.pre_spk_times = cumsum(exprnd(mean_isi,sim.T/mean_isi*5,1));
    data.pre_spk_times = data.pre_spk_times(data.pre_spk_times<sim.T);
else % inhomogeneous Poisson
    data.pre_spk_times = find(poissrnd(sim.pPreSpike)>0);
    data.pre_spk_times = data.pre_spk_times+rand(size(data.pre_spk_times))-.5;
    data.pre_spk_times = sort(data.pre_spk_times*sim.dt);
end

data.pre_spk_vec = getSpkMat(data.pre_spk_times,sim.dt,sim.T,1);

% coupling term
% sim.mprops.basis{2} = getBasis('rcos',sim.mprops.nfilt,sim.mprops.delay,50,0);
% sim.b=sim.b+randn(size(sim.b))/norm(sim.b)/10;

% Synaptic Connection
sim.syn = @(ts) max(0,ts-sim.alpha_dt)/sim.alpha_tau.*exp(1-max(0,ts-sim.alpha_dt)/sim.alpha_tau);
syn_kern = sim.syn(linspace(0, 1, 1/sim.dt));
Xc = filter(syn_kern, 1, data.pre_spk_vec');
wt = sim.wt_long(1);


% history filter
x0 = linspace(0,1,1/sim.dt);
kern_hist = exp(-x0/sim.hist_tau);
sim.post_alph = [log(sim.postBaseRate) kern_hist*sim.hist_beta];

% Generate postsynaptic spikes
[firingsPost,sim.g] = simLNP_stdp_v2(sim.post_alph,sim.dt,sim.T,wt,Xc,...
    data.pre_spk_vec,sim.stdp_params);

data.post_spk_times = firingsPost(:,1);
data.post_spk_vec = getSpkMat(data.post_spk_times,sim.dt,sim.T,1);
data.dt = sim.dt;
data.T = sim.T;

end
