function postData = deleteSing(preData)

facIndx = find(squeeze(sum(sum(preData.short_param_track_fac(:, :, :)==0, 1), 2)) == 0);
depIndx = find(squeeze(sum(sum(preData.short_param_track_dep(:, :, :)==0, 1), 2)) == 0);
nullIndx = find(squeeze(sum(sum(preData.short_param_track_null(:, :, :)==0, 1), 2)) == 0);

indx = intersect(intersect(facIndx, depIndx), nullIndx);

postData = preData;
% postData.mse_dep = preData.mse_dep(indx,:);
% postData.mse_fac = preData.mse_fac(indx,:);
% postData.mse_null = preData.mse_null(indx,:);
postData.seed_seq = preData.seed_seq(indx);

postData.short_param_track_dep = preData.short_param_track_dep(:, :, indx);
postData.short_param_track_fac = preData.short_param_track_fac(:, :, indx);
postData.short_param_track_null = preData.short_param_track_null(:, :, indx);

postData.tAlpha_track_dep = preData.tAlpha_track_dep(indx);
postData.tAlpha_track_fac = preData.tAlpha_track_fac(indx);
postData.tAlpha_track_null = preData.tAlpha_track_null(indx);

postData.tauAlpha_track_dep = preData.tauAlpha_track_dep(indx, 1);
postData.tauAlpha_track_fac = preData.tauAlpha_track_fac(indx, 1);
postData.tauAlpha_track_null = preData.tauAlpha_track_null(indx, 1);


end