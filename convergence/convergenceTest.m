% addpath(genpath('/home/gaw19004/showConvergence/helper'));

addpath(genpath('D:/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('D:/GitHub/GBLM_SMOOTH/core'));

%%
clc;clear all;close all;
T = 500;
dt = 0.001; pPreSpike = 5*dt;
t_alpha = 0.005; tau_alpha = 0.005;
Q = 1e-7;
seed = 123;

Nq_true = 5; Nq_fit = 5;
Nm = 450; Nst = 50;
iter = 30;

trueParam = [0 1 2 3 4]'*(0.06);
wt_long = ones(1, T/dt)';
lamObj = 10;

beta0 = beta0Adjust(T, dt, lamObj, pPreSpike, wt_long,...
t_alpha, tau_alpha, trueParam, Nq_true, Nm, Nst, seed);

rng(seed);
preSpike = binornd(1, pPreSpike, 1, round(T/dt));
isi = diff(find(preSpike>0)*dt);
N = length(preSpike);

% wt_short
Bm0 = getBasis('rcos', Nq_true, Nm, Nst,0);
Bm = padarray(Bm0,[0 5000]);
Bm=Bm(:,5001:end);
Bm_dt = Bm(:,round(isi*1000));
s=zeros(Nq_true,N);
for m =1: Nq_true
    s(m,preSpike>0)=[0 Bm_dt(m,:)];
end

tau= .2; % timeconstant for decay
x0 = linspace(0,1,1000);
kern_stp = exp(-x0/tau);
e = filter(kern_stp,1,s');
wt_short = 1 + e*trueParam;

% alpha function
syn = @(ts) max(0,ts-t_alpha)/tau_alpha.*exp(1-max(0,ts-t_alpha)/tau_alpha);
x0 = linspace(0, 1, 1000);
kern_c = syn(x0);
Xc = filter(kern_c, 1, preSpike');

% generate post spike
lam = exp(beta0 + wt_long.*wt_short.*Xc);
postSpike = spike_poiss2(T, dt, lam);

population{1} = find(preSpike ~= 0 )*dt;
population{2} = find(postSpike ~= 0 )*dt;

%%
rstart = 10;

[dynamic_L1, static_L1, diagnose_L1] =...
    smooth_gblm_conv(population, preSpike, postSpike, Nq_fit, Nm, Nst, Q, ...
    'longFirst', true, 'rstart', rstart);

[dynamic_S1, static_S1, diagnose_S1] =...
    smooth_gblm_conv(population, preSpike, postSpike, Nq_fit, Nm, Nst, Q, ...
    'longFirst', false, 'rstart', rstart);

save('convergence.mat');
%%
clc;clear all;close all;
load('convergence.mat')

hold on
plot(diagnose_L1.dev, 'r')
plot(diagnose_S1.dev, 'b:')
hold off

% well, hard to see anything
% plot the interested estimations iteration by iteration

invesVar = dynamic_S1.beta0;   % interested estimation. e.g. S1 = short first, beta0 
invesIter = 1;  % interested iteration

hold on
for k = 1:rstart
    plot(invesVar(:, invesIter, k));
end
hold off

% modification function
hold on
for k = 1:rstart
    plot(1 + Bm0'*dynamic_S1.wt_short_param(:, invesIter, k))
end
hold off




