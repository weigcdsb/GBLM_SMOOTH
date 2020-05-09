function [preSpike, postSpike, isi,...
beta0, wt_short, wt_long,...    
W_beta0_final, W_wt_long_final, beta0_final, wt_long_final,...
wt_short_final, wt_short_param_final, covB_final...
, static, dev, llhd, k] = simulation2(T, dt, pPreSpike, beta0, wt_long,...
t_alpha, tau_alpha, trueParam, Nq_true, Nq_fit, Nm, Nst, seed, iter, Q)

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

population{1} = find(preSpike ~= 0 )*dt;
population{2} = find(postSpike ~= 0 )*dt;

[dynamic, static, varSmooth, llhd, k, dev, covB_final] = ...
    smooth_gblm(population, preSpike, postSpike, Nq_fit, Nm, Nst, Q,...
    'iter', iter);

[W_beta0_final, W_wt_long_final, beta0_final, wt_long_final,...
    wt_short_final, wt_short_param_final] = findk(k,...
    varSmooth, dynamic);


end
