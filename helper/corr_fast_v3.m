% Fast cross/auto-correlation for point processes...
%
%  st1  Sorted spike times for neuron 1
%  st2  Sorted spike times for neuron 2
% 
%  Ta   Lower bound for differences
%  Tb   Upper bound for differences
%  bin  histogram construction
%
%  d        histogram
%  deltaT   matrix of differences [diff n1 n2]
%
%
% Based on ccc.m from Memming Park
% CC-BY Ian Stevenson

function [d,deltaT] = corr_fast_v3(st1,st2,Ta,Tb,bin)

n1 = length(st1);
n2 = length(st2);

if n1 == 0 || n2 == 0
    d = [];
    deltaT = [];
    return;
elseif ~issorted(st1) || ~issorted(st2)
    warn('Input spike times should be sorted')
end

% rough estimate of # of time difference required (assuming independence)
maxTTT = (Tb-Ta);
eN = ceil(prod([n1 n2]) * maxTTT * 2 / min(st1(end), st2(end)));
deltaT = zeros(10 * eN, 3);

% Compute all the time differences within the range [Ta,Tb]
lastStartIdx = 1;
k = 1;
for n = 1:n1
    incIdx = 0;
    for m = lastStartIdx:n2
        timeDiff = st2(m) - st1(n);
        if timeDiff >= Ta
            if incIdx==0
                incIdx = m;
            end
            if timeDiff <= Tb
                deltaT(k,:) = [timeDiff n m];
                k = k + 1;
            else % this is the ending point
                break;
            end
        end
    end
    if incIdx>0
        lastStartIdx = incIdx;
    end
end
deltaT = deltaT(1:(k-1),:);
edges = linspace(Ta,Tb,bin);
d = histc(deltaT(:,1),edges);