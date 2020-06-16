function crosscorrDiag(fit, data, varargin)

sim = struct;
if (~isempty(varargin))
    c = 1 ;
    while c <= length(varargin)
        switch varargin{c}
            case {'sim'}
                sim = varargin{c+1};
        end % switch
        c = c + 2;
    end % for
end % if
strIndx = isempty(fieldnames(sim));

subplot(1,3,1)
edges = linspace(-.02,.02,101);
n = histc(fit.deltaT(:,1),edges);
bar(edges+mean(diff(edges))/2,n,1,'EdgeColor','none');
hold on
x0 = linspace(-.02,.02,1001);
plot(x0,exp(fit.syn(x0))*mean(n), 'r')
if strIndx == 0;plot(x0,exp(sim.syn(x0))*mean(n), 'b');end;
hold off
box off

subplot(1,3,2)
[aagram,~] = corrFast(data.pre_spk_times,data.pre_spk_times,min(edges),max(edges),length(edges));
[~,i]=max(aagram);
aagram(i)=NaN;
bar(edges+mean(diff(edges))/2,aagram,1,'EdgeColor','none');
box off

subplot(1,3,3)
[apgram,~] = corrFast(data.post_spk_times,data.post_spk_times,min(edges),max(edges),length(edges));
[~,i]=max(apgram);
apgram(i)=NaN;
bar(edges+mean(diff(edges))/2,apgram,1,'EdgeColor','none');
box off

if strIndx == 0
    A = [[sim.alpha_dt; sim.alpha_tau] fit.synParams.syn_params(1:2)]*1000;
    array2table(A, 'VariableNames', {'True', 'fit'})
end


end