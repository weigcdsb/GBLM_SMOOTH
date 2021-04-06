function fit = noSTP(data, hist_tau, hist_beta,...
    stp_Nq, stp_Nm, stp_Ns, Q, synParams)


fit.hist_tau = hist_tau;
fit.hist_beta = hist_beta;

% parameters for GBLM for STP
fit.stp_Nq = stp_Nq;
fit.stp_Nm = stp_Nm;
fit.stp_Ns = stp_Ns;
fit.stp_tau= 1;

% parameters for synaptic connection
fit.syn_hyper_params.coupling_timescale = 0.05;
fit.syn_hyper_params.bin_width = (0.05)*data.dt;
fit.syn_hyper_params.baseline_nsplines = 4;

% parameters for smoothing
fit.F = eye(2);
fit.Q = Q;

% Synaptic Connection
fit.synParams = synParams;
fit = synConEst(data,fit);

fit.W = zeros(2, 2, data.vecN);
offset = log(data.dt) + fit.hist*fit.hist_beta;
alph = glmfit([fit.Xc],data.post_spk_vec,'poisson','Offset',offset);
fit.wt_long = ones(data.vecN, 1)*alph(2);
fit.beta0 = ones(data.vecN,1)*alph(1);

W0 = eye(2)*0.1;
b0 = [fit.beta0(1, :) fit.wt_long(1, :)]';
X = [ones(data.vecN,1) fit.Xc];
offset = fit.hist*fit.hist_beta;
[b, W, ~, ~] = ppasmoo_poissexp(data.post_spk_vec, X, b0, W0, fit.F, fit.Q, offset);

fit.beta0 = b(1,:)';
fit.wt_long = b(2,:)';
fit.W = W;

end