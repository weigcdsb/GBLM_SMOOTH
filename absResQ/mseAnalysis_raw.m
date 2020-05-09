% addpath(genpath('/home/gaw19004/18/helper'));

addpath(genpath('D:/GitHub/GBLM_SMOOTH/helper'));
addpath(genpath('D:/GitHub/GBLM_SMOOTH/core'));

%%
clc;clear all;close all;
saveName = 'mseResults_raw.mat';

T = 600;
dt = 0.001;
pPreSpike = 5*dt;
Nq_fit = 5;
grandSeed = 816;
rng(grandSeed);
seed_seq = randperm(5000, 1000);

Q_wt_long_seq = [0 logspace(-7, -3, 5)];

tAlpha_track_fac = zeros(length(seed_seq), 1);
tAlpha_track_dep = zeros(length(seed_seq), 1);
tAlpha_track_null = zeros(length(seed_seq), 1);

tauAlpha_track_fac = zeros(length(seed_seq), 1);
tauAlpha_track_dep = zeros(length(seed_seq), 1);
tauAlpha_track_null = zeros(length(seed_seq), 1);

short_param_track_fac = zeros(Nq_fit, length(Q_wt_long_seq), length(seed_seq));
short_param_track_dep = zeros(Nq_fit, length(Q_wt_long_seq), length(seed_seq));
short_param_track_null = zeros(Nq_fit, length(Q_wt_long_seq), length(seed_seq));

% for hpc
% poolobj = gcp('nocreate');
% delete(poolobj);
% num_proc=str2num(getenv('SLURM_NTASKS'));
% parpool('local',num_proc);

% for this laptop
poolobj = gcp('nocreate');
delete(poolobj);
parpool('local',2);

parfor s = 1:length(seed_seq)
    
    [short_param_track_fac(:, :, s), short_param_track_dep(:, :, s), short_param_track_null(:, :, s),...
    tAlpha_track_fac(s), tauAlpha_track_fac(s),...
    tAlpha_track_dep(s), tauAlpha_track_dep(s),...
    tAlpha_track_null(s), tauAlpha_track_null(s), ~] =...
    mseAnalysisInner(T, pPreSpike, Q_wt_long_seq, seed_seq(s));

end

save(saveName);