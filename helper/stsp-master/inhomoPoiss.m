% random rate generator

% gets :
% T : time
% nknots : number of knots
% fr : average firing rate

function [tspike,lambda] = inhomoPoiss(T,dt,nknots,fr)
% nknots = 50;
% T=50;
% fr=10;
% dt=.0001;
y=randn(1,nknots);
xi=linspace(0,T,T/dt);
x=linspace(0,T,nknots);
yi=spline(x,y,xi); % cubic spline interpolation, query points = xi
lambda = exp(yi)/mean(exp(yi))*fr; % make log(lambda) normally distributed,
% use mean(exp(yi)) for normalization 
tspike=spike_poiss2(T,dt,lambda);
% sum(tspike)

%% plot
% stem(tspike,'.')
% hold on
% plot(lambda/max(lambda),'r','LineWidth',2)
