function e = getE(seed, pPreSpike, T, dt, Nq_true, Nm, Nst)

rng(seed);
preSpike = binornd(1, pPreSpike, 1, round(T/dt));
isi = diff(find(preSpike>0)*dt);
N = length(preSpike);

% wt_short
Bm0 = getBasis('rcos', Nq_true, Nm, Nst,0);
Bm = padarray(Bm0,[0 5000]);
Bm=Bm(:,5001:end);
Bm_dt = Bm(:,round(isi*1000));
s=zeros(Nq_true,N);
for m =1: Nq_true
    s(m,preSpike>0)=[0 Bm_dt(m,:)];
end

tau= .2; % timeconstant for decay
x0 = linspace(0,1,1000);
kern_stp = exp(-x0/tau);
e = filter(kern_stp,1,s');

end