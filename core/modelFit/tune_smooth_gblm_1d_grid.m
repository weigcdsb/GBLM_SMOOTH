function [fit, fit_trace, Qvec,...
    qbllhd_pred, qwllhd_pred, Qopt] =...
    tune_smooth_gblm_1d_grid(preSpk, postSpk, varargin)
% TUNE_SMOOTH_GBLM is a Q tuned version for SMOOTH_GBLM. This function
% optimizes Q by seraching over a grid of Q. the default lower bound of Q
% is 1e-9 and the default upper bound of Q is 1e-3. The default number of
% searched Q points is 10. To turn on the plot, set llhdPlot = true.

doFit = true;
nq = 10;
QLB = 1e-9;
QUB = 1e-3;
llhdPlot = false;

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
            case {'nq'}
                nq = varargin{c+1};
            case {'QLB'}
                QLB = varargin{c+1};
            case {'QUB'}
                QUB = varargin{c+1};
            case {'llhdPlot'}
                llhdPlot = varargin{c+1};
            case {'doFit'}
                doFit = varargin{c+1};
        end % switch
        c = c + 2;
    end % for
end % if

%% Synaptic Connection
fit = synConEst(data,fit);

%% Q tune
Qvec = logspace(log10(QLB), log10(QUB), nq);

qbllhd_pred = ones(1, nq)*NaN;
for q=1:nq
    fprintf('Qbeta0 %02i/%02i...', q, nq)
    fit.Q = [Qvec(q) 0; 0 0];
    fit.doFiltOnly=true;
    fit.noSTP = true;
    fit.iter = 1;
    [fit, ~] = loopCore(data, fit);
    qbllhd_pred(q)=fit.llhd_pred;
end
[~, qb_indx] = max(qbllhd_pred);

fprintf('\n')
qwllhd_pred = ones(1, nq)*NaN;
for q=1:nq
    fprintf('Qwtlong %02i/%02i...', q, nq)
    fit.Q = [Qvec(qb_indx) 0; 0 Qvec(q)];
    fit.doFiltOnly=true;
    fit.noSTP = true;
    fit.iter = 1;
    [fit,~] = loopCore(data, fit);
    qwllhd_pred(q)=fit.llhd_pred;
end
[~, qw_indx] = max(qwllhd_pred);

fprintf('\n')
Qopt = [Qvec(qb_indx) 0; 0 Qvec(qw_indx)];
fit_trace = fit;
if doFit
    fit.Q = Qopt;
    fit.doFiltOnly = false;
    fit.noSTP = false;
    fit.iter = iter;
    [fit,fit_trace] = loopCore(data, fit);
end


if llhdPlot
    subplot(1,2,1)
    semilogx(Qvec,qbllhd_pred, 'b')
    ylabel('log likelihood')
    xlabel('Q')
    title('beta0')
    subplot(1,2,2)
    semilogx(Qvec,qwllhd_pred, 'b')
    xlabel('Q')
    title('wtlong')
end



end
