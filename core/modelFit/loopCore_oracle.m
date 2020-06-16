function [fit, fit_trace] = loopCore_oracle(data, fit)

% initialization
fit.wt_short_param = ones(fit.stp_Nq,1)*NaN;
fit.W = zeros(2, 2, data.vecN);
fit.llhd_pred=NaN;
fit.dev=NaN;
fit.covB = zeros(fit.stp_Nq)*NaN;

offset = log(data.dt) + fit.hist*fit.hist_beta;
alph = glmfit([fit.Xc],data.post_spk_vec,'poisson','Offset',offset);
fit.wt_long = ones(data.vecN, 1)*alph(2);
fit.beta0 = ones(data.vecN,1)*alph(1);
fit.wt_short = ones(data.vecN, 1);

lam = exp(fit.beta0 + fit.wt_long.*fit.wt_short.*fit.Xc + fit.hist*fit.hist_beta)*data.dt;
fit.llhd = sum(-lam + log((lam+(lam==0))).*(data.post_spk_vec'));

dev_prev=fit.dev;
fit_trace(1)=fit;
c=2;

for k = 1:fit.iter
    fprintf('Iter %02i...',k)
    
    if sum(isfield(fit, {'oracle_beta0', 'oracle_wt_long'})) == 2
        fit.beta0 = fit.oracle_beta0;
        fit.wt_long = fit.oracle_wt_long;
        fit.W = zeros(2, 2, data.vecN);
    elseif isfield(fit, 'oracle_beta0')
        fit.beta0 = fit.oracle_beta0;
        fit.F = 1;
        if k==1; fit.Q = fit.Q(2, 2); end
        W0 = 0.1;
        b0 = fit.wt_long(1, :);
        X = fit.wt_short.*fit.Xc;
        offset = fit.beta0 + fit.hist*fit.hist_beta;
        [b, W, fit, breakFlag] = filtOrSmooth(data, fit, W0, b0, X, offset);
        if breakFlag; break; end
        fit.wt_long = b';
        fit.W = zeros(2, 2, data.vecN);
        fit.W(2, 2, :) = W;
        
    elseif isfield(fit, 'oracle_wt_long')
        fit.wt_long = fit.oracle_wt_long;
        fit.F = 1;
        if k==1; fit.Q = fit.Q(1, 1); end
        W0 = 0.1;
        b0 = fit.beta0(1, :);
        X = ones(data.vecN,1);
        offset = fit.wt_long.*fit.wt_short.*fit.Xc + fit.hist*fit.hist_beta;
        [b, W, fit, breakFlag] = filtOrSmooth(data, fit, W0, b0, X, offset);
        if breakFlag; break; end
        fit.beta0 = b';
        fit.W = zeros(2, 2, data.vecN);
        fit.W(1, 1, :) = W;
        
    else
        W0 = eye(2)*0.1;
        b0 = [fit.beta0(1, :) fit.wt_long(1, :)]';
        X = [ones(data.vecN,1) (fit.wt_short.*fit.Xc)];
        offset = fit.hist*fit.hist_beta;
        [b, W, fit, breakFlag] = filtOrSmooth(data, fit, W0, b0, X, offset);
        if breakFlag; break; end
        fit.beta0 = b(1,:)';
        fit.wt_long = b(2,:)';
        fit.W = W;
    end
    
    lam = exp(fit.beta0 + fit.wt_long.*fit.wt_short.*fit.Xc + fit.hist*fit.hist_beta)*data.dt;
    fit.llhd = sum(-lam + log((lam+(lam==0))).*(data.post_spk_vec'));
    
    fprintf('%.02f...',-2*fit.llhd);
    fit_trace(c)=fit;
    c=c+1;
    
    
    % Update STP Parameters
    if isfield(fit,'oracle_stp_B')
        fit.wt_short = 1+fit.stp_X*fit.oracle_stp_B;
    else
        offset = fit.beta0 + log(data.dt) + fit.hist*fit.hist_beta;
        W_short = repmat(fit.wt_long.*fit.Xc, 1, fit.stp_Nq).*fit.stp_X;
        
        if isfield(fit, 'oracle_wt_long')
            [alph,fit.dev,~] = glmfit(W_short,data.post_spk_vec,'poisson','Offset',offset);
            fit.wt_short_param = alph(2:end);
            fit.wt_short = 1 + fit.stp_X*fit.wt_short_param;
            if ~isfield(fit, 'oracle_beta0'); fit.beta0 = fit.beta0 + alph(1); end
            
        else
            [alph,fit.dev,~] = glmfit([fit.Xc.*fit.wt_long W_short],data.post_spk_vec,'poisson','Offset',offset);
            fit.wt_short_param = alph(3:end)/alph(2);
            fit.wt_short = 1 + fit.stp_X*fit.wt_short_param;
            if ~isfield(fit, 'oracle_beta0'); fit.beta0 = fit.beta0 + alph(1); end
            fit.wt_long = fit.wt_long*alph(2);
        end
        
    end
    
    lam = exp(fit.beta0 + fit.wt_long.*fit.wt_short.*fit.Xc + fit.hist*fit.hist_beta)*data.dt;
    fit.llhd = sum(-lam + log((lam+(lam==0))).*(data.post_spk_vec'));
    
    if isfield(fit,'oracle_stp_B')
        fit.dev = 2*sum(data.post_spk_vec'.*(log((data.post_spk_vec'+(data.post_spk_vec'==0)))...
            - log((lam+(lam==0))))- (data.post_spk_vec' - lam));
    end
    
    fprintf('%.02f...',-2*fit.llhd);
    
    fit_trace(c)=fit;
    c=c+1;
    fprintf('\n');
    
    if k>5 && abs((fit.dev-dev_prev)/dev_prev)<fit.toleranceValue;break;end
    dev_prev=fit.dev;
    
end

% re-fit the GBLM to get var-cov matrix
if isfield(fit,'oracle_stp_B')
    fit.covB = zeros(fit.stp_Nq);
else
    offset = fit.beta0 + log(data.dt) + fit.hist*fit.hist_beta;
    W_short = repmat(fit.wt_long.*fit.Xc, 1, fit.stp_Nq).*fit.stp_X;
    
    [~,~,stat] = glmfit(W_short,data.post_spk_vec,'poisson','Offset',offset,'constant','off');
    fit.covB=stat.covb;
    
    lam = exp(fit.beta0 + fit.wt_long.*fit.wt_short.*fit.Xc + fit.hist*fit.hist_beta)*data.dt;
    fit.llhd = sum(-lam + log((lam+(lam==0))).*(data.post_spk_vec'));
    fit_trace(c)=fit;
    
end

end


function [b, W, fit, breakFlag] = filtOrSmooth(data, fit, W0, b0, X, offset)

breakFlag = false;

if isfield(fit,'doFiltOnly') && fit.doFiltOnly
    [b, W, lam_pred] = ppafilt_poissexp(data.post_spk_vec, X, b0, W0, fit.F, fit.Q, offset);
    lam_pred = lam_pred'*data.dt;
    fit.llhd_pred = sum(-lam_pred + log((lam_pred+(lam_pred==0))).*(data.post_spk_vec'));
else
    [b, W, lam_pred, ~] = ppasmoo_poissexp(data.post_spk_vec, X, b0, W0, fit.F, fit.Q, offset);
    lam_pred = lam_pred'*data.dt;
    fit.llhd_pred = sum(-lam_pred + log((lam_pred+(lam_pred==0))).*(data.post_spk_vec'));
end

[warnmsg, msgid] = lastwarn;
if strcmp(msgid,'MATLAB:illConditionedMatrix')
    disp(warnmsg);
    breakFlag = true;
end

end
