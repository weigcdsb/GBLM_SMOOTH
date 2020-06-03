function plot_simfit(sim,fit)


figure(1)
subplot(3,1,1)
plot(sim.beta0)
hold on
plot(fit.beta0)
plot(fit.beta0+sqrt(squeeze(fit.W(1,1,:))),'r:');
plot(fit.beta0-sqrt(squeeze(fit.W(1,1,:))),'r:');
hold off
ylabel('beta0')

subplot(3,1,2)
plot(sim.wt_long)
hold on
plot(fit.wt_long)
plot(fit.wt_long+sqrt(squeeze(fit.W(2,2,:))),'r:');
plot(fit.wt_long-sqrt(squeeze(fit.W(2,2,:))),'r:');
hold off
ylabel('wt-long')

subplot(3,1,3)
plot(1+sim.stp_X*sim.stp_B)
hold on
plot(fit.wt_short)
hold off
ylabel('wt-short')


se_mod_fn = sqrt(diag(sim.stp_basis'*fit.covB*sim.stp_basis));

figure(2)
plot(sim.stp_basis'*sim.stp_B)
hold on
plot(sim.stp_basis'*fit.wt_short_param)
plot(sim.stp_basis'*fit.wt_short_param + se_mod_fn,'r:')
plot(sim.stp_basis'*fit.wt_short_param - se_mod_fn,'r:')
hold off