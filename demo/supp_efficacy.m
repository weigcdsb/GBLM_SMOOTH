
isi = [Inf; diff(Tpre)];
quantiles = prctile(isi,linspace(0,100,25));
qx = quantiles(1:end-1)+diff(quantiles)/2;


efficacy = zeros(1, length(quantiles)-1);
model_efficacy = zeros(1, length(quantiles)-1);
syn_kern = sim.syn(linspace(0, 1, 1/data.dt));
sk = exp(syn_kern);
sk=sk(sk>1);

wt = (1+sim.stp_X*sim.stp_B).*sim.wt_long;
for q=1:length(quantiles)-1
    qspk = Tpre(isi>=quantiles(q) & isi<quantiles(q+1));
    
    dbase = corr_fast_v3(qspk,Tpost,-.008,-0.004,20);
    d = corr_fast_v3(qspk,Tpost,0.004,.008,20);
    
    n = zeros(1,data.T/data.dt);
    n(ceil(qspk/data.dt))=1;
    model_efficacy(q) = mean(exp(wt(n>0)));
    
    
    efficacy(q) = (sum(d) - sum(dbase))/length(qspk);
end

effPlot = figure;
clf
hold on
plot(qx*1000, efficacy, 'o', 'Color', [0, 0.4470, 0.7410],...
    'LineWidth', 2, 'markerfacecolor', [0, 0.4470, 0.7410])
plot(qx*1000, (model_efficacy)/nanmean(model_efficacy)*nanmean(efficacy),...
    'Color', [0, 0.4470, 0.7410], 'LineWidth', 3)
set(gca,'FontSize',15, 'LineWidth', 1.5,'TickDir','out')
box off
hold off

set(effPlot,'PaperUnits','inches','PaperPosition',[0 0 4 3])