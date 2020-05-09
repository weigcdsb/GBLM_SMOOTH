function [dynamic, static, varSmooth, llhd, k, dev, covB] =...
    smooth_gblm(population, s_pre, s_post, Nq, Nm, Nst, Q, varargin)

dt = 1e-3;
N = length(s_post);
isi = diff(find(s_pre>0)*dt);

iter = 30;
toleranceValue=.02;

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

% initialization
% save the track for debugging, may not include for final version
wt_short = ones(N, iter);
wt_short_param = zeros(Nq, iter);
wt_long = ones(N, iter);
beta0 = repmat(log(mean(s_post)/dt), N, iter);
dev = zeros(1, iter);
W_track = zeros(2, 2, N, iter);
llhd = zeros(1, iter);

for k = 2:iter
    
    % 1. update beta0, wt_long_amp & wt_long_loc: adaptive smoothing
    W0 = eye(2)*0.1;
    F = eye(2);
    b0 = [beta0(1, k-1) wt_long(1, k-1)]';
    
    y = s_post';
    X = [ones(N,1) (wt_short(:, k-1).* Xc)];
    [b, W, ~] = ppasmoo2_poissexp(y, X, b0, W0, F, Q);
    bw = b';
    
    % if warning, break;
    [warnmsg, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:illConditionedMatrix')
        disp(warnmsg);
        beta0(:, k) = beta0(:, k-1);
        wt_long(:, k) = wt_long(:, k-1);
        W_track(:, :, :, k) = W_track(:, :, :, k-1);
        wt_short_param(:, k) = wt_short_param(:, k-1);
        wt_short(:, k) = wt_short(:, k-1);
        beta0(:,k) = beta0(:,k-1);
        wt_long(:,k) = wt_long(:,k-1);
        llhd(:, k) = llhd(:, k-1);
        dev(k) = dev(k-1);
        break;
    end
    
    % force wt_long_amp to be non-negative
    beta0(:, k) = bw(:, 1);
    temp = bw(:, 2);
    temp(temp <= 0) = 0;
    wt_long(:, k) = temp;
    W_track(:, :, :, k) = W;
    
    % 2. update wt_short: glm
    offset = beta0(:, k) + log(dt) + log(postHist2);
    W_short = repmat(wt_long(:, k) .* Xc, 1, Nq) .*e;
    
    [alph,dev(k),~] = glmfit([Xc.*wt_long(:,k) W_short],s_post,'poisson','Offset',offset);
    
    wt_short_param(:, k) = alph(3:end);
    wt_short(:, k) = 1 + e*wt_short_param(:, k)/alph(2);
    
    beta0(:,k) = beta0(:, k) + alph(1);
    wt_long(:,k) = wt_long(:,k)*alph(2);
    
    lam = exp(beta0(:, k) + wt_long(:, k).*wt_short(:, k).*Xc + log(postHist2))*dt;
    llhd(:, k) = sum(-lam + log((lam+(lam==0))).*(s_post'));
    
    if abs(dev(k)-dev(k-1))<toleranceValue && k>5;break;end
    
end

% re-fit the model to get var-cov matrix
offset_f = beta0(:, k) + wt_long(:, k).*Xc + log(dt) + log(postHist2);
W_short_f = repmat(wt_long(:, k) .* Xc, 1, Nq) .*e;

[alph_f,dev(k),stats_f] = glmfit(W_short_f,s_post,'poisson','Offset',offset_f);

wt_short_param(:, k) = alph_f(2:end);
wt_short(:, k) = 1 + e*wt_short_param(:, k);
beta0(:,k) = beta0(:, k) + alph_f(1);

covB = stats_f.covb;

lam = exp(beta0(:, k) + wt_long(:, k).*wt_short(:, k).*Xc + log(postHist2))*dt;
llhd(:, k) = sum(-lam + log((lam+(lam==0))).*(s_post'));


dynamic.beta0 = beta0;
dynamic.wt_short = wt_short;
dynamic.wt_short_param = wt_short_param;
dynamic.wt_long = wt_long;

static.tAlpha = synParams.syn_params(1);
static.tauAlpha = synParams.syn_params(2);

varSmooth = W_track;

end