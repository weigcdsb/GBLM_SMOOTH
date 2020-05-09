function [dynamic, static, diagnose] =...
    smooth_gblm_conv(population, s_pre, s_post, Nq, Nm, Nst, Q, varargin)

dt = 1e-3;
N = length(s_post);
isi = diff(find(s_pre>0)*dt);

iter = 30;
toleranceValue=.02;

longFirst = true;
rstart = 5; % number of random randomstart
grandStartSeed = 123;

hyper_params.coupling_timescale = 0.05;
hyper_params.bin_width = (0.05)*dt;
hyper_params.baseline_nsplines = 4; 
hyper_params.thr =.05; % useless

if (~isempty(varargin))
    c = 1 ;
    while c <= length(varargin)
        switch varargin{c}
            case {'iter'}
                iter = varargin{c+1};
            case {'toleranceValue'}
                toleranceValue = varargin{c+1};
            case {'hyper_params'}
                hyper_params = varargin{c+1};
            case {'longFirst'}
                longFirst = varargin{c+1};
            case {'rstart'}
                rstart = varargin{c+1};
            case {'grandStartSeed'}
                grandStartSeed = varargin{c+1};
        end % switch
        c = c + 2;
    end % for
end % if

%% 1. get the shape of alpha function (cross-correlation)

% population{2} = population{2}+(rand(size(population{2}))- 1)*dt/2;
[syn,~,synParams,~]= synapse_xcorr(population,hyper_params);

%% 2. convolution to get Xc
x0 = linspace(0, 1, 1000);
kern_c = syn(x0);
Xc = filter(kern_c, 1, s_pre');

%% 3. get basis for wt_short (LM_regression)
Bm0 = getBasis('rcos',Nq, Nm, Nst,0);
Bm = padarray(Bm0,[0 5000]);
Bm=Bm(:,5001:end);
Bm_dt = Bm(:,round(isi*1000));
s=zeros(Nq,N);
for m =1: Nq
    s(m,s_pre>0)=[0 Bm_dt(m,:)];
end

%% 4. transient effect of stp carried by exp function
tau= .2; % timeconstant for decay
x0 = linspace(0,1,1000);
kern_stp = exp(-x0/tau);
e = filter(kern_stp,1,s');

%% 5. history for postSpike
postFilt = [0 ones(1, 5)];
postHist = filter(postFilt, 1, s_post');

postHist2 = abs(postHist - 1);
postHist2(postHist2 == 0) = exp(-16);


%% 6. model fitting
rng(grandStartSeed);
rStartSeed_seq = randperm(10*rstart, rstart);

BETA0 = zeros(N, 5, rstart);
WT_SHORT = zeros(N, 5, rstart);
WT_SHORT_PARAM = zeros(Nq, 5, rstart);
WT_LONG = zeros(N, 5, rstart);

W_BETA0 = zeros(N, rstart);
W_WT_LONG = zeros(N, rstart);

K = zeros(1, rstart);
LLHD = zeros(iter, rstart);
DEV = zeros(iter, rstart);
COVB = zeros(Nq + 1, Nq + 1, rstart);

% for hpc
% poolobj = gcp('nocreate');
% delete(poolobj);
% num_proc=str2num(getenv('SLURM_NTASKS'));
% parpool('local',num_proc);

% for this laptop
poolobj = gcp('nocreate');
delete(poolobj);
parpool('local', 2);

parfor m = 1:rstart
    
    [dyn, varS, LLHD(:, m), K(m), DEV(:, m), COVB(:, :, m)] = ...
        GBLMloop(rStartSeed_seq(m), longFirst, N, dt, Nq, iter, e, Q, s_post, Xc, postHist2, toleranceValue);
    
    BETA0(:, :, m) = dyn.beta0;
    WT_SHORT(:, :, m) = dyn.wt_short;
    WT_SHORT_PARAM(:, :, m) = dyn.wt_short_param;
    WT_LONG(:, :, m) = dyn.wt_long;
    
    W_BETA0(:, m) = varS.W_beta0;
    W_WT_LONG(:, m) = varS.W_wt_long; 

end

dynamic.seed_seq = rStartSeed_seq;
dynamic.beta0 = BETA0;
dynamic.wt_short = WT_SHORT;
dynamic.wt_short_param = WT_SHORT_PARAM;
dynamic.wt_long = WT_LONG;
dynamic.W_beta0 = W_BETA0;
dynamic.W_wt_long = W_WT_LONG;
dynamic.covB = COVB;

diagnose.k = K;
diagnose.llhd = LLHD;
diagnose.dev = DEV;

static.tAlpha = synParams.syn_params(1);
static.tauAlpha = synParams.syn_params(2);

end