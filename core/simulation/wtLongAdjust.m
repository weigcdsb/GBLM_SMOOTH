function wt_long = wtLongAdjust(meanwt, sim)

rng(sim.seed);


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
data.pre_spk_vec = zeros(1,sim.vecN);
data.pre_spk_vec(round(data.pre_spk_times/sim.dt))=1;


% Short-Term Plasticity
isi = diff(find(data.pre_spk_vec>0)*sim.dt);
sim.stp_basis = getBasis('rcos', sim.stp_Nq, sim.stp_Nm, sim.stp_Ns,0);
Bm = padarray(sim.stp_basis,[0 max(round(isi/sim.dt))],'post');
Bm_dt = Bm(:,round(isi/sim.dt));
s=zeros(sim.stp_Nq,sim.vecN);
for m =1:sim.stp_Nq
    s(m,data.pre_spk_vec>0)=[0 Bm_dt(m,:)];
end
x0 = linspace(0,1,1/sim.dt);
kern_stp = exp(-x0/sim.stp_tau);
sim.stp_X = filter(kern_stp,1,s');
wt_short = 1 + sim.stp_X*sim.stp_B;

longAmp_dep = meanwt/mean(wt_short);
wt_long = repmat(longAmp_dep, 1, sim.T/sim.dt)';

end