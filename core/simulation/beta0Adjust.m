function beta0 = beta0Adjust(lamObj, sim)
% approximately

rng(sim.seed);

if length(sim.pPreSpike)==1 % homogeneous Poisson
    mean_isi = sim.dt/sim.pPreSpike;
    pre_spk_times = cumsum(exprnd(mean_isi,sim.T/mean_isi*5,1));
    pre_spk_times = pre_spk_times(pre_spk_times<sim.T);
else % inhomogeneous Poisson
    pre_spk_times = find(poissrnd(sim.pPreSpike)>0);
    pre_spk_times = pre_spk_times+rand(size(pre_spk_times))-.5;
    pre_spk_times = sort(pre_spk_times*sim.dt);
end
pre_spk_vec = zeros(1,sim.vecN);
pre_spk_vec(round(pre_spk_times/sim.dt))=1;


% Synaptic Connection
sim.syn = @(ts) max(0,ts-sim.alpha_dt)/sim.alpha_tau.*exp(1-max(0,ts-sim.alpha_dt)/sim.alpha_tau);
syn_kern = sim.syn(linspace(0, 1, 1/sim.dt));
Xc = filter(syn_kern, 1, pre_spk_vec');

% Short-Term Plasticity
isi = diff(find(pre_spk_vec>0)*sim.dt);
sim.stp_basis = getBasis('rcos', sim.stp_Nq, sim.stp_Nm, sim.stp_Ns,0);
Bm = padarray(sim.stp_basis,[0 max(round(isi/sim.dt))],'post');
Bm_dt = Bm(:,round(isi/sim.dt));
s=zeros(sim.stp_Nq,sim.vecN);
for m =1:sim.stp_Nq
    s(m,pre_spk_vec>0)=[0 Bm_dt(m,:)];
end
x0 = linspace(0,1,1/sim.dt);
kern_stp = exp(-x0/sim.stp_tau);
sim.stp_X = filter(kern_stp,1,s');
wt_short = 1 + sim.stp_X*sim.stp_B;

beta0_amp = log(lamObj/(mean(exp(sim.wt_long.*wt_short.*Xc))));
beta0 = repmat(ceil(beta0_amp*2)/2, 1, sim.vecN)';

end