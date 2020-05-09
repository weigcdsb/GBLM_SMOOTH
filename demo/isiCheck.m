function isiCheck(trueParam, Nq_true, Nm, Nst)

% wt_short
Bm0 = getBasis('rcos', Nq_true, Nm, Nst,0);
plot(1 + Bm0'*trueParam, 'r');

end