function [dynamic, varSmooth, llhd, k, dev, covB] = ...
    GBLMloop(seed, longFirst, N, dt, Nq, iter, e, Q, s_post, Xc, postHist2, toleranceValue)

rng(seed);

% initialization
wt_short_param = zeros(Nq, iter);
wt_short_param(:, 1) = unifrnd(-0.5, 0.5, 5, 1);
wt_short = 1 + e*wt_short_param;

wt_long = ones(N, iter);
wt_long(:, 1) = ones(N, 1)*unifrnd(-2, 2);

beta0 = zeros(N, iter);
beta0(:, 1) = ones(N, 1)*unifrnd(-5, 5);

W_track = zeros(2, 2, N, iter);
dev = zeros(1, iter);
llhd = zeros(1, iter);

if longFirst == true
    for k = 2:iter
        % 1. long-term effect
        [beta0, wt_long, W_track] =...
            longEffect(k, N, Q, beta0, wt_long, W_track, wt_short, s_post, Xc);
        [warnmsg, msgid] = lastwarn;
        if strcmp(msgid,'MATLAB:illConditionedMatrix')
            disp(warnmsg);
            [beta0, wt_long, W_track, wt_short_param, wt_short, llhd, dev] =...
                breakAssign(k, beta0, wt_long, W_track, wt_short_param,...
                wt_short, llhd, dev);
            break;
        end
        % 2. short-term effect
        [beta0, wt_short_param, wt_short, wt_long] =...
            shortEffect(k, dt, Nq, s_post, beta0, wt_long, wt_short_param,...
            wt_short, postHist2, Xc, e);
        
        lam = exp(beta0(:, k) + wt_long(:, k).*wt_short(:, k).*Xc + log(postHist2))*dt;
        dev(k) = 2*sum(s_post'.*(log((s_post'+(s_post'==0))) - log((lam+(lam==0))))...
            - (s_post' - lam));
        llhd(k) = sum(-lam + log((lam+(lam==0))).*(s_post') - 1);
        
        if abs(dev(k)-dev(k-1))<toleranceValue && k>5;break;end
    end
else
    for k = 2:iter
        % 1. short-term effect
        [beta0, wt_short_param, wt_short, wt_long] =...
            shortEffect(k, dt, Nq, s_post, beta0, wt_long, wt_short_param,...
            wt_short, postHist2, Xc, e);
        % 2. long-term effect
        [beta0, wt_long, W_track] =...
            longEffect(k, N, Q, beta0, wt_long, W_track, wt_short, s_post, Xc);
        [warnmsg, msgid] = lastwarn;
        if strcmp(msgid,'MATLAB:illConditionedMatrix')
            disp(warnmsg);
            [beta0, wt_long, W_track, wt_short_param, wt_short, llhd, dev] =...
                breakAssign(k, beta0, wt_long, W_track, wt_short_param,...
                wt_short, llhd, dev);
            break;
        end
        
        lam = exp(beta0(:, k) + wt_long(:, k).*wt_short(:, k).*Xc + log(postHist2))*dt;
        dev(k) = 2*sum(s_post'.*(log((s_post'+(s_post'==0))) - log((lam+(lam==0))))...
            - (s_post' - lam));
        llhd(k) = sum(-lam + log((lam+(lam==0))).*(s_post') - 1);
        
        if abs(dev(k)-dev(k-1))<toleranceValue && k>5;break;end
    end
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
llhd(:, k) = sum(-lam + log((lam+(lam==0))).*(s_post') - 1);

% to show independence of random start
lam = exp(beta0(:, 1) + wt_long(:, 1).*wt_short(:, 1).*Xc + log(postHist2))*dt;
dev(1) = 2*sum(s_post'.*(log((s_post'+(s_post'==0))) - log((lam+(lam==0))))...
    - (s_post' - lam));
llhd(1) = sum(-lam + log((lam+(lam==0))).*(s_post') - 1);

varSmooth.W_beta0 = squeeze( W_track(1, 1, :, k));
varSmooth.W_wt_long = squeeze( W_track(2, 2, :, k));

dynamic.beta0 = beta0(:, [1 2 3 4 k]);
dynamic.wt_short = wt_short(:, [1 2 3 4 k]);
dynamic.wt_short_param = wt_short_param(:, [1 2 3 4 k]);
dynamic.wt_long = wt_long(:, [1 2 3 4 k]);

end