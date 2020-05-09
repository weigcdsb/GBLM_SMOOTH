function [beta0, wt_long, W_track] =...
    longEffect(k, N, Q, beta0, wt_long, W_track, wt_short, s_post, Xc)

% W0 = W_track(:, :, 1, k-1);
W0 = eye(2, 2)*0.1;
F = eye(2);
b0 = [beta0(1, k-1) wt_long(1, k-1)]';

y = s_post';
X = [ones(N,1) (wt_short(:, k-1).* Xc)];
[b, W, ~] = ppasmoo2_poissexp(y, X, b0, W0, F, Q);
bw = b';

% if warning, break;
[~, msgid] = lastwarn;
if strcmp(msgid,'MATLAB:illConditionedMatrix')
    return
end


beta0(:, k) = bw(:, 1);
% temp = bw(:, 2);
% if mean(temp) >= 0
%     LTPSign = 1;
% else
%     LTPSign = -1;
% end
% temp(temp <= 0) = 0;
% wt_long(:, k) = LTPSign*temp;

wt_long(:, k) = bw(:, 2);
W_track(:, :, :, k) = W;

end













