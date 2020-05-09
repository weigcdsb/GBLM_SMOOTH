function spike=spike_poiss2(T,dt,rr)
% just for simulation of a spike train
N = round(T/dt);
spike=zeros(1,N);
t=linspace(0,T,N);

for i=6:N
    if rand<1-exp(-rr(i)*dt) && sum(spike(i-5:i-1))==0;
        % P(X = 0) = exp(-lambda(i)*dt) = exp(-rr(i)*dt);
        % 1- exp(-rr(i)*dt) = P(X ~= 0)
        % sum(spike(i-5:i-1))==0, take refractory of neurons into
        % consideration;
       spike(i)=1;
    end
end