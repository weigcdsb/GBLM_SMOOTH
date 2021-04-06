addpath(genpath('D:\GitHub\GBLM_SMOOTH'));

%%
Tpre = data.pre_spk_times;
Tpost = data.post_spk_times;
dbase = corr_fast_v3(Tpre,Tpost,-.008,-0.004,20);
d = corr_fast_v3(Tpre,Tpost,0.004,.008,20);

efficacy = (sum(d) - sum(dbase))/length(Tpre);
efficacy


