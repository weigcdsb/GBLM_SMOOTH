function [wt_short_param_track_fac, wt_short_param_track_dep, wt_short_param_track_null,...
    tAlpha_fac, tauAlpha_fac,...
    tAlpha_dep, tauAlpha_dep,...
    tAlpha_null, tauAlpha_null, e] =...
    mseAnalysisInner(T, pPreSpike, Q_wt_long_seq, seed)

dt = 0.001;
t_alpha = 0.005; tau_alpha = 0.005;

Nq_true = 5; Nq_fit = 5;
Nm = 450; Nst = 50;
iter = 30;

%% attain e
e = getE(seed, pPreSpike, T, dt, Nq_true, Nm, Nst);

%% analysis for parameters control

trueParam_seq = ...
    [[0 1 2 3 4]'*(0.06)...
    [0 -1 -2 -3 -4]'*(0.06)...
    [0 0 0 0 0]'];

% 1. modify wt_long to make mean(wt_long.*wt_short) the same
wt_long_fac = ones(1, T/dt)';
wt_long_dep = wtLongAdjust(T, dt, pPreSpike, wt_long_fac, trueParam_seq, Nq_true, Nm, Nst, seed);
wt_long_null = ones(1, T/dt)';

% 2. set mean(lam*dt) = 10Hz for post, to get beta0
lamObj = 10;

beta0_fac = beta0Adjust(T, dt, lamObj, pPreSpike, wt_long_fac,...
t_alpha, tau_alpha, trueParam_seq(:, 1), Nq_true, Nm, Nst, seed);
beta0_dep = beta0Adjust(T, dt, lamObj, pPreSpike, wt_long_dep,...
t_alpha, tau_alpha, trueParam_seq(:, 2), Nq_true, Nm, Nst, seed);
beta0_null = beta0Adjust(T, dt, lamObj, pPreSpike, wt_long_null,...
t_alpha, tau_alpha, trueParam_seq(:, 3), Nq_true, Nm, Nst, seed);

%%

nQ = length(Q_wt_long_seq);
Q = zeros(2, 2, nQ);

wt_short_param_track_fac = zeros(Nq_fit, nQ);
wt_short_param_track_dep = zeros(Nq_fit, nQ);
wt_short_param_track_null = zeros(Nq_fit, nQ);

for d = 1:nQ
    
    disp(d);
    Q(:, :, d) = diag([0 Q_wt_long_seq(d)]);
    
    
    %1. facilitation
    [~, ~, ~,...
        ~, ~, ~,...
        ~, ~, ~, ~,...
        ~, wt_short_fac_param_final, ~...
        , static, ~, ~, ~] = simulation2(T, dt, pPreSpike, beta0_fac, wt_long_fac,...
        t_alpha, tau_alpha, trueParam_seq(:, 1), Nq_true, Nq_fit, Nm, Nst, seed, iter, Q(:, :, d));
    
    
    wt_short_param_track_fac(:, d) = wt_short_fac_param_final;
    
    tAlphFit = static.tAlpha;
    tauAlphFit = static.tauAlpha;
    
    if d == 1
        tAlpha_fac = tAlphFit;
        tauAlpha_fac = tauAlphFit;
    end
    
    
    %2. depression
    [~, ~, ~,...
        ~, ~, ~,...
        ~, ~, ~, ~,...
        ~, wt_short_dep_param_final, ~...
        , static, ~, ~, ~] = simulation2(T, dt, pPreSpike, beta0_dep, wt_long_dep,...
        t_alpha, tau_alpha, trueParam_seq(:, 2), Nq_true, Nq_fit, Nm, Nst, seed, iter, Q(:, :, d));
    
    wt_short_param_track_dep(:, d) = wt_short_dep_param_final;
    
    tAlphFit = static.tAlpha;
    tauAlphFit = static.tauAlpha;
    
    if d == 1
        tAlpha_dep = tAlphFit;
        tauAlpha_dep = tauAlphFit;
    end
    
    %3. no plasticity
    [~, ~, ~,...
        ~, ~, ~,...
        ~, ~, ~, ~,...
        ~, wt_short_null_param_final, ~...
        , static, ~, ~, ~] = simulation2(T, dt, pPreSpike, beta0_null, wt_long_null,...
        t_alpha, tau_alpha, trueParam_seq(:, 3), Nq_true, Nq_fit, Nm, Nst, seed, iter, Q(:, :, d));
    
    wt_short_param_track_null(:, d) = wt_short_null_param_final;
    
    tAlphFit = static.tAlpha;
    tauAlphFit = static.tauAlpha;
    
    if d == 1
        tAlpha_null = tAlphFit;
        tauAlpha_null = tauAlphFit;
    end
    
end

end





