
%% Generate presynaptic spikes...
gmod.Npre = 1;
gmod.dt = 0.004;
gmod.isFitG = 0;
smod.dt = gmod.dt;
gmod.pre_alph = log(5);
gmod.maxT = 2*60*60;   % in s

% firings = simLNP_trsX(gmod.pre_alph*ones(gmod.Npre,1),gmod.dt,gmod.maxT);
% smod.S = getSpkMat({firings(:,1)},smod.dt,gmod.maxT,1);

% Fast simulation using exprnd... (this works since it's homogeneous)
firings = exprnd(1/exp(gmod.pre_alph),floor(exp(gmod.pre_alph)*gmod.maxT*2),gmod.Npre);
firings = cumsum(firings);
Tlist=[];
for i=1:gmod.Npre
    tmp = firings(:,i);
    Tlist{i} = tmp(tmp<gmod.maxT);
end
smod.S = getSpkMat(Tlist,smod.dt,gmod.maxT,1);

%% COUPLING PARAMETERS
gmod.mprops.nfilt = 5;
gmod.mprops.delay = 100/(smod.dt*1000);
gmod.mprops.basis{2} = getBasis('rcos',gmod.mprops.nfilt,gmod.mprops.delay,20,0);
% gmod.mprops.basis{2} = gmod.mprops.basis{2}(1,:);
gmod.b=[0 0.1 0.9 0.4 0.1]'/1.75;
gmod.b=gmod.b+randn(size(gmod.b))/norm(gmod.b)/10;

figure(1)
subplot(2,3,1)
t0 = (1:size(gmod.mprops.basis{2},2))*smod.dt*1000;
plot(t0,gmod.mprops.basis{2}')
hold on
plot(t0,exp(gmod.mprops.basis{2}'*gmod.b),'k')
hold off
drawnow
box off

%% SPIKE-HISTORY PARAMETERS
gmod.a = [-3 -0.01 -0.01 0.04 0.015]';
gmod.a=gmod.a+randn(size(gmod.a))/norm(gmod.a)/10;

figure(1)
subplot(2,3,2)
plot(t0,gmod.mprops.basis{2}')
hold on
plot(t0,exp(gmod.mprops.basis{2}'*gmod.a),'k')
hold off
drawnow
box off

gmod.post_alph = [log(5); gmod.mprops.basis{2}'*gmod.a]';
gmod.post_alph(2)=-100;
% gmod.post_alph(2:end) = gmod.post_alph(2:end)-mean(gmod.post_alph(2:end));
% gmod.post_alph = log(8);

%% STDP PARAMETERS

stdp_params.noise = 0; % in s
stdp_params.tau_forgetting = 60; % in s

stdp_params.type = 'dexp';
stdp_params.tau_plus = 20/1000; % in s
stdp_params.tau_minus = 20/1000; % in s
% stdp_params.A_plus = 0.001;
stdp_params.A_plus = 0.005;
% % stdp_params.A_minus = 1*stdp_params.A_plus;
% % stdp_params.A_minus = 1.39*stdp_params.A_plus;
stdp_params.A_minus = 1.38*stdp_params.A_plus;
% 
% stdp_params.type = 'ahebb';
% stdp_params.tau_plus = 20/1000; % in s
% stdp_params.tau_minus = 50/1000; % in s
% stdp_params.A_plus = 0.009;
% stdp_params.A_minus = 0.475*stdp_params.A_plus;

stdp_params.g_max = 50;
stdp_params.g_init = 1;

subplot(2,3,3)
plot(linspace(-0.1,0.1,100),getModFn(stdp_params,linspace(-0.1,0.1,100)));
box off

%% Check for E(dw)=0
% t0 = (1:size(gmod.mprops.basis{2},2))/1000;
% fw = stdp_params.A_plus*exp(-t0/stdp_params.tau_plus);
% psp = poisspdf(1,exp(gmod.post_alph(1) + gmod.mprops.basis{2}'*gmod.b)*smod.dt);
% ewp = fw*psp;
% 
% subplot(2,3,3)
% plot(t0,psp,t0,fw)
% hold on
% 
% t0 = -(1:size(gmod.mprops.basis{2},2))/1000;
% fw = -stdp_params.A_minus*exp(t0/stdp_params.tau_minus);
% psm = poisspdf(1,exp(repmat(gmod.post_alph(1),length(t0),1))*smod.dt);
% ewm = fw*ps;
% 
% plot(t0,psm,t0,fw)
% hold off
% 
% % stdp_params.A_minus = -ewp/ewm*stdp_params.A_minus;

%% Generate postsynaptic spikes...

X = getX(smod.S,gmod.mprops.basis{2},0,1,0)';
%%
[firingsPost,g] = simLNP_stdp_v2(gmod.post_alph,gmod.dt,gmod.maxT,gmod.b,X,smod.S,stdp_params);

smod.Tlist{1} = firingsPost(:,1);
for i=1:gmod.Npre
    smod.Tlist{i+1} = Tlist{i};
end
smod.S = getSpkMat(smod.Tlist,smod.dt,gmod.maxT,1);

%%
figure(1)
subplot(2,3,4:5)
cla
plot(linspace(0,length(g)*smod.dt/60,length(g)),g,'k')
box off
% %%
% % double check mod-fn
% subplot(2,3,6)
% Xstdp = getRelSpkTiming_full_indie(smod.Tlist{stdpPair(1)},smod.Tlist{stdpPair(2)},stdp_kern_vec,smod.dt,max(trRange));
% 
% %% clf
% dg = circshift(diff(g)',0);
% for i=1:size(Xstdp,2)
%     idx = (Xstdp(1:100000,i)>0);
%     plot(i,dg(idx),'ko')
%     hold on
% end
% hold off

%%
% gy=g(1:size(Xstdp,1));
% XstdpEF = expfilt(Xstdp,1/60/1000);
% bbb = glmfit(XstdpEF,gy');
% subplot(2,3,4:5)
% hold on
% plot(linspace(0,length(gy)/1000,length(gy)),[Xstdp(:,1)*0+1 XstdpEF]*bbb,'r')
% % plot(linspace(0,length(gy)/1000,length(gy)),gy,'r')
% hold off
% % %% Plot some nice things...
% 
% trange=[154 156];
% trange=[445 447];
% clf
% % subplot(2,1,1)
% tmp = [[firings firings*0-1]; [firingsPost]];
% drawRaster(tmp(find(tmp(:,1)>trange(1) & tmp(:,1)<trange(2)),:))
% xlim([150 155])
% box off
% set(gca,'XTick',[])
% set(gca,'YTick',[])
% % subplot(2,1,2)
% hold on
% plot(linspace(trange(1),trange(2),length(trange(1)/smod.dt:trange(2)/smod.dt)),X(trange(1)/smod.dt:trange(2)/smod.dt,:)*gmod.b.*g(trange(1)/smod.dt:trange(2)/smod.dt)','r')
% plot(linspace(trange(1),trange(2),length(trange(1)/smod.dt:trange(2)/smod.dt)),X(trange(1)/smod.dt:trange(2)/smod.dt,:)*gmod.b)
% axis tight
% hold off


%% A little post analysis
% 
% % If the simulation hit a boundary truncate
% ddown=min(find(g==stdp_params.g_max | g==0));
% if isempty(ddown),ddown=length(g); end
% 
% figure(1); clf
% if gmod.isFitG
%     % Make sure that the simulated synaptic strength matches spike pairs
%     stdpPair=[1 2];
%     trRange=[1 gmod.maxT/2/smod.dt];
%     tsRange=[gmod.maxT/2/smod.dt+1 gmod.maxT/smod.dt];
%     stdp_kern_vec = linspace(-0.1,0.1,15);
%     XstdpKernFull = getRelSpkTiming_full_indie(smod.Tlist{stdpPair(1)},smod.Tlist{stdpPair(2)},stdp_kern_vec,smod.dt,max(size(smod.S,2),max(tsRange)));
% 
%     d=([ones(size(XstdpKernFull(1:ddown,:),1),1) cumsum(XstdpKernFull(1:ddown,:))])\g(1:ddown)';
%     ghat = ([ones(size(XstdpKernFull(1:ddown,:),1),1) cumsum(XstdpKernFull(1:ddown,:))])*d;
%     subplot(1,4,1:3)
%     plot(ghat(1:ddown))
%     hold on
%     plot(g(1:ddown),'r')
%     hold off
%     corrcoef(ghat(1:ddown),g(1:ddown))
%     subplot(1,4,4)
%     plot(stdp_kern_vec(1:end-1)+diff(stdp_kern_vec)/2,d(2:end))
% else
%     ds = floor(linspace(1,length(g),10000));
%     plot(g(ds))
% end
% drawnow
% 
% ddown/length(g)
% trRange=[1 floor(ddown/2)];
% tsRange=[floor(ddown/2)+1 ddown];