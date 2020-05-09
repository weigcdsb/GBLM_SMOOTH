function wt_long_dep = wtLongAdjust(T, dt, pPreSpike, wt_long_fac,...
trueParam_seq, Nq_true, Nm, Nst, seed)

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

wt_short_fac = 1 + e*trueParam_seq(:,1); % facilitation
wt_short_dep = 1 + e*trueParam_seq(:,2); % depression

longAmp_dep = mean(wt_long_fac.*wt_short_fac)/mean(wt_short_dep);
wt_long_dep = repmat(longAmp_dep, 1, T/dt)';

end