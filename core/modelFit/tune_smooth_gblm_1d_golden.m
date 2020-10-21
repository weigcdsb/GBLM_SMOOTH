function [fit, fit_trace, Qopt] =...
    tune_smooth_gblm_1d_golden(preSpk, postSpk, varargin)

doFit = true;
QLB = 1e-9;
QUB = 1e-4;
QIter = 25;

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
iter = 30;
fit.toleranceValue= 1e-6;
fit.F = eye(2);

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
                iter = varargin{c+1};
            case {'toleranceValue'}
                fit.toleranceValue = varargin{c+1};
            case {'F'}
                fit.F = varargin{c+1};
            case {'QLB'}
                QLB = varargin{c+1};
            case {'QUB'}
                QUB = varargin{c+1};
            case {'tolX'}
                tolX = varargin{c+1};
            case {'QIter'}
                QIter = varargin{c+1};
            case {'doFit'}
                doFit = varargin{c+1};
            case {'synParams'}
                fit.synParams = varargin{c+1};
        end % switch
        c = c + 2;
    end % for
end % if

if exist('tolX','var') == 0; tolX = QLB*1e-1;end

%% Synaptic Connection
fit = synConEst(data,fit);

options = optimset('MaxIter', QIter, 'TolX', tolX);
f = @(Q) helper_1d(Q, fit, data, 0, 1);
QbEst = fminbnd(f, QLB, QUB, options);

f = @(Q) helper_1d(Q, fit, data, QbEst, 2);
QwEst = fminbnd(f, QLB, QUB, options);

Qopt = [QbEst 0; 0 QwEst];

fit.Q = Qopt;
fit_trace = fit;
if doFit
    fit.Q = Qopt;
    fit.doFiltOnly = false;
    fit.noSTP = false;
    fit.iter = iter;
    [fit,fit_trace] = loopCore(data, fit);
end


end


function neg_llhd_pred = helper_1d(Q, fit, data, Qother, idx)

if idx == 1
    fit.Q = diag([Q, Qother]);
elseif idx == 2
    fit.Q = diag([Qother, Q]);
end

fit.doFiltOnly=true;
fit.noSTP = true;
fit.iter = 1;
[fit, ~] = loopCore(data, fit);
neg_llhd_pred = -fit.llhd_pred;

end








