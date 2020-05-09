function llhdCheck(T, dt, pPreSpike, beta0, wt_long,...
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
postSpike = spike_poiss2(T, dt, lam);



lam_oracle = exp(beta0 + wt_long.*wt_short.*Xc);
llhd_oracle = sum(-lam_oracle*dt + log((lam_oracle+(lam_oracle==0))*dt).*(postSpike'));

lam_null = mean(postSpike)/dt;
llhd_null   = sum(-lam_null*dt + log((lam_null+(lam_null==0))*dt).*(postSpike'));

lam_sat = (postSpike)'/dt;
llhd_sat   = sum(-lam_sat*dt + log((lam_sat+(lam_sat==0))*dt).*(postSpike'));

line(xlim(),[1 1]*llhd_oracle,'Color','r');
hold on;
line(xlim(),[1 1]*llhd_null,'Color','k');
line(xlim(),[1 1]*llhd_sat,'Color','b');
hold off;


end