
function dw = getModFn(stdp_params,prepost)

if ~isfield(stdp_params,'type'),
    stdp_params.type = 'dexp';
end

if strcmp(stdp_params.type,'ahebb')
    % Hebbian with negative component
    dw = stdp_params.A_plus*normpdf(prepost,0,stdp_params.tau_plus)/normpdf(0,0,stdp_params.tau_plus) ...
        -stdp_params.A_minus*normpdf(prepost,0,stdp_params.tau_minus)/normpdf(0,0,stdp_params.tau_minus);
else
    % Basic STDP
    dw = stdp_params.A_plus*exp(prepost/stdp_params.tau_plus);
    dw(prepost>0) = -stdp_params.A_minus*exp(-prepost(prepost>0)/stdp_params.tau_minus);    
end