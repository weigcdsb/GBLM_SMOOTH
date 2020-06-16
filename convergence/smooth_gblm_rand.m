function [FIT,FIT_trace] = smooth_gblm_rand(preSpk, postSpk, varargin)

% data set up
data.dt = 0.001;
data.pre_spk_vec = preSpk;
data.post_spk_vec = postSpk;
data.vecN = length(data.pre_spk_vec);
data.pre_spk_times = find(preSpk ~= 0 )*data.dt;
data.post_spk_times = find(postSpk ~= 0 )*data.dt;

% parameters for history filter
fit.hist_tau = .01;
fit.hist_beta = -2;

% parameters for GBLM for STP
fit.stp_Nq = 5;
fit.stp_Nm = 450;
fit.stp_Ns = 50;
fit.stp_tau= 1;

% parameters for synaptic connection
fit.syn_hyper_params.coupling_timescale = 0.05;
fit.syn_hyper_params.bin_width = (0.05)*data.dt;
fit.syn_hyper_params.baseline_nsplines = 4;

% parameters for smoothing
fit.iter = 30;
fit.toleranceValue= 1e-6;
fit.F = eye(2);
fit.Q = eye(2)*1e-5;
fit.doFiltOnly = false;

rstart = 6;
startSeed = 1;

if (~isempty(varargin))
    c = 1 ;
    while c <= length(varargin)
        switch varargin{c}
            case {'hist_tau'}
                fit.hist_tau = varargin{c+1};
            case {'hist_beta'}
                fit.hist_beta = varargin{c+1};
            case {'stp_Nq'}
                fit.stp_Nq  = varargin{c+1};
            case {'stp_Nm'}
                fit.stp_Nm = varargin{c+1};
            case {'stp_Ns'}
                fit.stp_Ns = varargin{c+1};
            case {'stp_tau'}
                fit.stp_tau = varargin{c+1};
            case {'syn_hyper_params.coupling_timescale'}
                fit.syn_hyper_params.coupling_timescale = varargin{c+1};
            case {'syn_hyper_params.bin_width'}
                fit.syn_hyper_params.bin_width = varargin{c+1};
            case {'syn_hyper_params.baseline_nsplines'}
                fit.syn_hyper_params.baseline_nsplines = varargin{c+1};
            case {'iter'}
                fit.iter = varargin{c+1};
            case {'toleranceValue'}
                fit.toleranceValue = varargin{c+1};
            case {'F'}
                fit.F = varargin{c+1};
            case {'Q'}
                fit.Q = varargin{c+1};
            case {'doFiltOnly'}
                fit.doFiltOnly = varargin{c+1};
            case {'doOracleSTP'}
                fit.doOracleSTP = varargin{c+1};
            case {'oracle_stp_B'}
                fit.oracle_stp_B = varargin{c+1};
            case {'longTermFirst'}
                Lfirst = varargin{c+1};
        end % switch
        c = c + 2;
    end % for
end % if

%% Synaptic Connection
fit = synConEst(data,fit);

%% Smoothing_GBLM
% initialization
rng(startSeed)

FIT.wt_short_param = unifrnd(-0.5, 0.5, fit.stp_Nq, rstart);
FIT.wt_short = 1 + fit.stp_X*FIT.wt_short_param;
FIT.wt_long = ones(data.vecN, 1)*unifrnd(-3, 3, 1, rstart);
FIT.beta0 = ones(data.vecN, 1)*unifrnd(-4, 4, 1, rstart);

FIT.W = zeros(2, 2, data.vecN, rstart);
FIT.dev= ones(1, rstart)*NaN;
FIT.llhd= ones(1, rstart)*NaN;
FIT.maxIter = ones(1, rstart)*NaN;
FIT.covB = zeros(fit.stp_Nq, fit.stp_Nq, rstart)*NaN;

FIT_trace(1).wt_short_param = FIT.wt_short_param;
FIT_trace(1).wt_long = FIT.wt_long;
FIT_trace(1).beta0 = FIT.beta0;
FIT_trace(1).dev = FIT.dev;

for r = 1:rstart
    fprintf('randomStart %02i...',r)
    fitr = fit;
    fitr.wt_short_param =  FIT_trace(1).wt_short_param(:, r);
    fitr.wt_short = 1 + fitr.stp_X*fitr.wt_short_param;
    fitr.wt_long = FIT_trace(1).wt_long(:, r);
    fitr.beta0 = FIT_trace(1).beta0(:, r);
    fitr.W = zeros(2, 2, data.vecN);
    fitr.dev= NaN;
    fitr.llhd= NaN;
    fitr.maxIter = NaN;
    fitr.covB = zeros(fitr.stp_Nq)*NaN;
    dev_p = fitr.dev;
    
    for k = 2:fitr.iter
        fprintf('Iter %02i...',k)
        
        if ~Lfirst
            fitr = shortEffect(fitr, data);
            [fitr, warnFlag] = longEffect(fitr, data);
            if warnFlag; break; end
            
        else
            [fitr, warnFlag] = longEffect(fitr, data);
            if warnFlag; break; end
            fitr = shortEffect(fitr, data);
            
        end
        
        [fitr, FIT_trace] = llhdDev(fitr, FIT_trace, data, k, r);
        dev_c = fitr.dev;
        if k>5 && abs((dev_c-dev_p)/dev_p)<fitr.toleranceValue;break;end
        dev_p = fitr.dev;
    end
    fprintf('\n')
    
    fitr.maxIter = k;
    offset = fitr.beta0 + log(data.dt) + fitr.hist*fitr.hist_beta;
    W_short = repmat(fitr.wt_long.*fitr.Xc, 1, fitr.stp_Nq).*fitr.stp_X;
    
    [~,~,stat] = glmfit(W_short,data.post_spk_vec,'poisson','Offset',offset,'constant','off');
    fitr.covB = stat.covb;
    
    lam = exp(fitr.beta0 + fitr.wt_long.*...
        fitr.wt_short.*fitr.Xc + fitr.hist*fitr.hist_beta)*data.dt;
    fitr.llhd = sum(-lam + log((lam+(lam==0))).*(data.post_spk_vec'));
    
    FIT_trace(k).wt_short_param(:, r) = fitr.wt_short_param;
    FIT_trace(k).wt_long(:, r) = fitr.wt_long;
    FIT_trace(k).beta0(:, r) = fitr.beta0;
    FIT_trace(k).dev(r) = fitr.dev;
    
    FIT.wt_short_param(:, r) = fitr.wt_short_param;
    FIT.wt_short(:, r) = fitr.wt_short;
    FIT.wt_long(:, r) = fitr.wt_long;
    FIT.beta0(:, r) = fitr.beta0;
    
    FIT.W(:, :, :, r) = fitr.W;
    FIT.dev(r) = fitr.dev;
    FIT.llhd(r)= fitr.llhd;
    FIT.maxIter(r) = fitr.maxIter;
    FIT.covB(:, :, r) = fitr.covB;
    
    clear fitr;
    
    
end


end


function [fit, warnFlag] = longEffect(fit, data)

warnFlag = false;
W0 = eye(2)*0.1;
b0 = [fit.beta0(1) fit.wt_long(1)]';
X = [ones(data.vecN,1) (fit.wt_short.*fit.Xc)];
offset = fit.hist*fit.hist_beta;

[b, W, ~, ~] = ppasmoo_poissexp(data.post_spk_vec, X, b0, W0, fit.F, fit.Q, offset);

[warnmsg, msgid] = lastwarn;
if strcmp(msgid,'MATLAB:illConditionedMatrix')
    disp(warnmsg);
    warnFlag = true;
end

fit.beta0 = b(1,:)';
fit.wt_long = b(2,:)';
fit.W = W;

end


function fit = shortEffect(fit, data)
offset = fit.beta0 + log(data.dt) + fit.hist*fit.hist_beta;
W_short = repmat(fit.wt_long.*fit.Xc, 1, fit.stp_Nq).*fit.stp_X;

[alph,~,~] = glmfit([fit.Xc.*fit.wt_long W_short],...
    data.post_spk_vec,'poisson','Offset',offset);

fit.wt_short_param = alph(3:end)/alph(2);
fit.wt_short = 1 + fit.stp_X*fit.wt_short_param;
fit.beta0 = fit.beta0 + alph(1);
fit.wt_long = fit.wt_long*alph(2);

end

function [fit, FIT_trace] = llhdDev(fit, FIT_trace, data, k, r)

lam = exp(fit.beta0 + fit.wt_long.*...
    fit.wt_short.*fit.Xc + fit.hist*fit.hist_beta)*data.dt;
fit.llhd = sum(-lam + log((lam+(lam==0))).*(data.post_spk_vec'));
fit.dev = 2*sum(data.post_spk_vec'.*(log((data.post_spk_vec'+(data.post_spk_vec'==0)))...
    - log((lam+(lam==0))))- (data.post_spk_vec' - lam));

FIT_trace(k).wt_short_param(:, r) = fit.wt_short_param;
FIT_trace(k).wt_long(:, r) = fit.wt_long;
FIT_trace(k).beta0(:, r) = fit.beta0;
FIT_trace(k).dev(r) = fit.dev;

end
