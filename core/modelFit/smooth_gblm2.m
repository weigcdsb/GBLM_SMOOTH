function [fit,fit_trace] = smooth_gblm2(preSpk, postSpk, varargin)
% SMOOTH_GBLM estimates LTP, STP and baseline effects simultaneously. The
% baseline and LTP are estimated by adaptive smoothing, while the STP is
% estimated by GBLM. The data and fitting parameters need to be passed into
% the function. The structure array 'data' should contain pre- and post-
% synaptic spike trains and spike time. The delta t is fixed as 1e-3.All
% parameters for fitting are contained in the structure array called 'fit':
%
%   'hist_tau':     time constant for decay in refractory history filter.
%                   The default is 0.01.
%   'hist_beta':    amplitude of refractory history filter. The default is
%                   -2
%
%   'stp_Nq':       number of basis vectors for GBLM. The default is 5.
%   'stp_Nm':       two times stp_Nm is the peak (i.e. center) of the last
%                   raised cosine basis vectors. The default is 450.
%   'stp_Ns':       offset for nonlinear stretching of x axis. The defaut
%                   is 50.
%   'stp_tau':      time constant for decay in STP. The default is 1.
%
%   'syn_hyper_params.coupling_timescale':
%                   time scale of coupling when estimating the alpha
%                   function (for synpatic connection). The defalt is 0.05.
%   'syn_hyper_params.bin_width':
%                   bin width for estimating the alpha function. The
%                   default is 5*1e-5.
%   'syn_hyper_params.baseline_nsplines':
%                   number of basis for estimation alaph function. The
%                   default is 4.
%
%   'iter':         number of iteration for smoothing-gblm. The default is
%                   30.
%   'toleranceValue':
%                   criteria for convergence. Use dev_c and dev_p to denote
%                   deviance for current and previous iterations. When
%                   abs((dev_c-dev_p)/dev_p) < toleranceValue, stop the
%                   iteration. The default is 1e-6.
%   'F':            transition matrix for adaptive filtering. The default
%                   is identity matrix to avoid bias.
%   'Q':            covariance matrix for Gaussian white noise. The default
%                   is eye(2)*1e-5.
%   'doFilterOnly': if true, do filtering without smoothing. The default is
%                   false.
%   'doOracle':     use true parameters. The default is false. If true,
%                   must provide at least one of 'oracle_beta0',
%                   'oracle_wt_long' or 'oracle_stp_B'.
%   'oracle_beta0': true baseline
%   'oracle_wt_long':
%                   true LTP
%   'oracle_stp_B': true parameters for STP.

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
fit.doOracle = false;

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
            case {'doOracle'}
                fit.doOracle = varargin{c+1};
            case {'oracle_beta0'}
                fit.oracle_beta0 = varargin{c+1};
            case {'oracle_wt_long'}
                fit.oracle_wt_long = varargin{c+1};
            case {'oracle_stp_B'}
                fit.oracle_stp_B = varargin{c+1};
            case {'synParams'}
                fit.synParams = varargin{c+1};
        end % switch
        c = c + 2;
    end % for
end % if


%% Synaptic Connection
fit = synConEst(data,fit);

disp(fit.synParams.syn_params(1))
disp(fit.synParams.syn_params(2))

%% Smoothing_GBLM
if fit.doOracle
    if sum(isfield(fit,{'oracle_beta0','oracle_wt_long','oracle_stp_B'})) == 0
        fit_trace = NaN;
        fprintf('Must provide at least one oracle value');
    else
        [fit, fit_trace] = loopCore_oracle(data, fit);
    end
    
else
    [fit, fit_trace] = loopCore2(data, fit);
end

end