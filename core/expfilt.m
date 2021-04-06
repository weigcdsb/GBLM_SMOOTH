function y = expfilt(x,tau,overwrite)

y=x*0;
y(1,:) = x(1,:);
if nargin<3 || overwrite==false
    for t=2:size(x,1)
        y(t,:) = (1-tau)*y(t-1,:)+x(t,:);
    end
%      y=y*tau;
else
    for t=2:size(x,1)
        if x(t,:)>0
            y(t,:) = x(t,:);
        else
            y(t,:) = (1-tau)*y(t-1,:);
        end
    end
end