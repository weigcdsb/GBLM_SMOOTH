function fit = synConEst(data,fit)

isi = diff(find(data.pre_spk_vec>0)*data.dt);
fit.stp_basis = getBasis('rcos', fit.stp_Nq, fit.stp_Nm, fit.stp_Ns,0);
Bm = padarray(fit.stp_basis,[0 max(round(isi/data.dt))],'post');
Bm_dt = Bm(:,round(isi/data.dt));
s=zeros(fit.stp_Nq, data.vecN);
for m =1:fit.stp_Nq
    s(m,data.pre_spk_vec>0)=[0 Bm_dt(m,:)];
end
x0 = linspace(0,1,1/data.dt);
kern_stp = exp(-x0/fit.stp_tau);
fit.stp_X = filter(kern_stp,1,s');

% alpha function
if (~isfield(fit, 'synParams'))
    [fit.syn,fit.deltaT,fit.synParams,~]= synapse_xcorr_conv({data.pre_spk_times,data.post_spk_times},fit.syn_hyper_params);
else
    fit.syn = @(ts) max(0,ts-fit.synParams.syn_params(1))/fit.synParams.syn_params(2).*exp(1-max(0,ts-fit.synParams.syn_params(1))/fit.synParams.syn_params(2));
end

fit.syn_kern = fit.syn(linspace(0, .1, .1/data.dt));
fit.Xc = filter(fit.syn_kern, 1, data.pre_spk_vec');

% history filter for refractory effect
kern_hist = [0 exp(-x0/fit.hist_tau)];
fit.hist = filter(kern_hist,1,data.post_spk_vec');


end