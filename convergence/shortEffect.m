function [beta0, wt_short_param, wt_short, wt_long] =...
    shortEffect(k, dt, Nq, s_post, beta0, wt_long, wt_short_param,...
            wt_short, postHist2, Xc, e)

offset = beta0(:, k) + log(dt) + log(postHist2);
W_short = repmat(wt_long(:, k) .* Xc, 1, Nq) .*e;

[alph,~,~] = glmfit([Xc.*wt_long(:,k) W_short],s_post,'poisson','Offset',offset);

wt_short_param(:, k) = alph(3:end);
wt_short(:, k) = 1 + e*wt_short_param(:, k)/alph(2);

beta0(:,k) = beta0(:, k) + alph(1);
wt_long(:,k) = wt_long(:,k)*alph(2);

end