function isiCheck2(Nq_true, Nq_fit, Nm, Nst, trueParam, wt_short_param_final,covB_final)

Bm0_true = getBasis('rcos', Nq_true, Nm,Nst,0);
Bm0_fit = getBasis('rcos', Nq_fit, Nm,Nst,0);

seM = sqrt(diag(Bm0_fit' * covB_final(2:end, 2:end) * Bm0_fit));

plot(1 +Bm0_true'*trueParam, 'r');
hold on;
plot(1 + Bm0_fit'*wt_short_param_final, 'b');
plot(1 + Bm0_fit'*wt_short_param_final - seM, 'b:');
plot(1 + Bm0_fit'*wt_short_param_final + seM, 'b:');

hold off;

end