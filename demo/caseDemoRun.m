function caseDemoRun(T, beta0, wt_long, trueParam, seed, Q,...
    saveDir, fileName, varargin)

dt = 0.001; pPreSpike = 6*dt;
t_alpha = 0.005; tau_alpha = 0.005;

Nq_true = 5; Nq_fit = 5;
Nm = 450; Nst = 50;
iter = 30;

if (~isempty(varargin))
    c = 1 ;
    while c <= length(varargin)
        switch varargin{c}
            case {'dt'}
                dt = varargin{c+1};
                
            case {'pPreSpike'}
                pPreSpike = varargin{c+1};
                
            case {'t_alpha'}
                t_alpha = varargin{c+1};
                
            case {'tau_alpha'}
                tau_alpha = varargin{c+1};
                
            case {'Nq_true'}
                Nq_true = varargin{c+1};
                
            case {'Nq_fit'}
                Nq_fit = varargin{c+1};
                
            case {'Nm'}
                Nm = varargin{c+1};
                
            case {'Nst'}
                Nst = varargin{c+1};
                
            case {'iter'}
                iter = varargin{c+1};
        end % switch
        c = c + 2;
    end % for
end % if



tic;
[preSpike, postSpike, isi,...
beta0, wt_short, wt_long,...    
W_beta0_final, W_wt_long_final, beta0_final, wt_long_final,...
wt_short_final, wt_short_param_final, covB_final...
, static, dev, llhd, k] = simulation2(T, dt, pPreSpike, beta0, wt_long,...
t_alpha, tau_alpha, trueParam, Nq_true, Nq_fit, Nm, Nst, seed, iter, Q);
elpsT = toc;

save(fullfile(saveDir, fileName));


end