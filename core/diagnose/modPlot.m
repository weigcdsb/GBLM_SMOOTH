function modPlot(sim, fit)

se_mod_fn = sqrt(diag(sim.stp_basis'*fit.covB*sim.stp_basis));

plot(sim.stp_basis'*sim.stp_B, 'k')
hold on
plot(sim.stp_basis'*fit.wt_short_param, 'r')
plot(sim.stp_basis'*fit.wt_short_param + se_mod_fn,'r:')
plot(sim.stp_basis'*fit.wt_short_param - se_mod_fn,'r:')
hold off

end