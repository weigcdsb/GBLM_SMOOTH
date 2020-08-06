function modPlot(sim, fit, width, width_dash)

se_mod_fn = sqrt(diag(sim.stp_basis'*fit.covB*sim.stp_basis));

plot(1 + sim.stp_basis'*sim.stp_B, 'k', 'LineWidth', width)
hold on
plot(1 + sim.stp_basis'*fit.wt_short_param, 'r', 'LineWidth', width)
plot(1 + sim.stp_basis'*fit.wt_short_param + se_mod_fn,'r:', 'LineWidth', width_dash)
plot(1 + sim.stp_basis'*fit.wt_short_param - se_mod_fn,'r:', 'LineWidth', width_dash)
hold off

end