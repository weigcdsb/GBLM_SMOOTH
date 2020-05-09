function simPlot(beta0_final, W_beta0, b0,...
    wt_long_final, W_wt_long,...
    wt_long, wt_short_final, wt_short, covB_final, e)


idx = 1:size(beta0_final);

subplot(3, 1, 1)
plot(idx, beta0_final,...
    idx, beta0_final + sqrt(W_beta0), 'b:',...
    idx, beta0_final - sqrt(W_beta0), 'b:',...
    idx, b0, 'r');
title('beta_0')

subplot(3, 1, 2)
plot(idx, wt_long_final,...
    idx, wt_long_final + sqrt(W_wt_long), 'b:',...
    idx, wt_long_final - sqrt(W_wt_long), 'b:',...
    idx, wt_long, 'r');
title('wt_{long}')

vb = covB_final(2:end, 2:end);
se_wt_short = zeros(size(e, 1), 1);
for k = 1:size(e, 1)
    se_wt_short(k) = sqrt(e(k,:)*vb*(e(k,:))');
end



subplot(3, 1, 3)
plot(idx, wt_short_final,...
    idx, wt_short_final + se_wt_short, 'b:',...
    idx, wt_short_final - se_wt_short, 'b:',...
    idx, wt_short, 'r');
title('wt_{short}')

end