function lamCehck(T, dt, pPreSpike, beta0, wt_long,...
t_alpha, tau_alpha, trueParam, Nq_true, Nm, Nst, seed)

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
wt_short = 1 + e*trueParam;

% alpha function
syn = @(ts) max(0,ts-t_alpha)/tau_alpha.*exp(1-max(0,ts-t_alpha)/tau_alpha);
x0 = linspace(0, 1, 1000);
kern_c = syn(x0);
Xc = filter(kern_c, 1, preSpike');

% generate post spike
lam = exp(beta0 + wt_long.*wt_short.*Xc);

plot(lam*dt);

end