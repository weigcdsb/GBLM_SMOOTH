function [firings,g] = simLNP_stdp(alph,delta,lim,beta,X,S,stdp_params)

C = size(alph,1);
Ms = (size(alph,2)-1)/C;
if ~isfinite(Ms), Ms=0; end

if nargin<4, beta=[]; end

maxrows = ceil(lim/delta);

% Preallocate for speed
firings = zeros(maxrows,2);
fcount = zeros(C,1);
I = [1 zeros(1,C*Ms)]';
t = 1;

u = zeros(C,1); % clocks
tau = exprnd(1,C,1); % next spike time

% STDP parameters
g = zeros(size(S,1),size(X,1));
g(:,1) = stdp_params.g_init;

% store old spikes
buff_size = 10;
prespk = ones(buff_size,size(S,1))*-Inf;
posspk = ones(buff_size,C)*-Inf;

ind = (1:C)';
ttp = ones(C,1);
smax =1;
cont = true;

tic
while cont
    cont = (t*delta)<lim;

    % Update conditional intensities
    if isempty(beta)
        u = u + delta*exp(alph*I);
    else
        u = u + delta*exp(alph*I + sum(X(t,:).*beta'.*g(:,t)'));
    end
    fired = u>tau;
    if sum(fired)
        if (smax+sum(fired) > size(firings,1))
            firings = [firings; zeros(maxrows,2)]; % allocate some more space
        end

        firings(smax:(smax+sum(fired)-1),:) = [t*ttp(1:sum(fired)) ind(fired)];
        fcount = fcount+fired;
        smax = sum(fcount);

        tau(fired) = exprnd(1,sum(fired),1);
        u(fired) = 0;
    end

    % Update synaptic strengths
    prespk = prespk+delta;
    posspk = posspk+delta;
    if ~isempty(beta)
        g(:,t+1) = g(:,t)-delta*(g(:,t)-1)/stdp_params.tau_forgetting;
        % Update for pre-syn spike
        if any(S(:,t))
            prespk(:,S(:,t)>0) = circshift(prespk(:,S(:,t)>0),1);
            prespk(:,S(:,t)>0) = 0;
            g(S(:,t)>0,t+1) = g(S(:,t)>0,t)+sum(getModFn(stdp_params,posspk+randn(size(posspk))*stdp_params.noise))';
        end
        % Update for post-syn spike
        if any(fired)
            posspk(:,fired) = circshift(posspk(:,fired),1);
            posspk(:,fired) = 0;
            g(:,t+1) = g(:,t)+sum(getModFn(stdp_params,-prespk+randn(size(prespk))*stdp_params.noise))';
        end
        % Cap out of bounds weights
        if any(g(:,t+1)>stdp_params.g_max | g(:,t+1)<0)
            g(g(:,t+1)>stdp_params.g_max,t+1) = stdp_params.g_max;
            g(g(:,t+1)<0,t+1) = 0;
            cont=false;
        end
    end

    %     I = circshift(I',1)';
    I(2:end) = I(1:end-1);
    I(1) = 1;
    if Ms>0,
        I(2:Ms:C*Ms+1) = fired;
    end

    if (mod(t,100000) == 0)
        fprintf('t:%04i >> %04i spikes (%04i min)\n',t*delta,round(mean(fcount)),min(fcount));
        % Sometings wrong...
        if min(fcount) == 0
            cont = false;
        end
    end

    t = t+1;
end

g = g(:,1:end-1);
firings = firings(1:sum(fcount)-1,:);
firings(:,1) = firings(:,1)*delta;
fprintf('t:%04i >> %04i spikes\n',t*delta,min(fcount));
toc