function [beta0, wt_long, W_track, wt_short_param, wt_short, llhd, dev] =...
    breakAssign(k, beta0, wt_long, W_track, wt_short_param,...
    wt_short, llhd, dev)

beta0(:, k) = beta0(:, k-1);
wt_long(:, k) = wt_long(:, k-1);
W_track(:, :, :, k) = W_track(:, :, :, k-1);
wt_short_param(:, k) = wt_short_param(:, k-1);
wt_short(:, k) = wt_short(:, k-1);
llhd(:, k) = llhd(:, k-1);
dev(k) = dev(k-1);

end