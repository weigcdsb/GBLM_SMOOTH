function [W_beta0_final, W_wt_long_final, beta0_final, wt_long_final,...
    wt_short_final, wt_short_param_final] = findk(k, varSmooth, dynamic)

W_final = varSmooth(:, :, :, k);
W_beta0_final = squeeze(W_final(1, 1, :));
W_wt_long_final = squeeze(W_final(2, 2, :));


beta0_final = dynamic.beta0(:, k);
wt_long_final = dynamic.wt_long(:, k);


wt_short_final = dynamic.wt_short(:, k);
wt_short_param_final = dynamic.wt_short_param(:, k);


end