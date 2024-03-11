function [nanvec, x, y, m, mediation_covariates, additionalM] = nanremove(x, y, m, mediation_covariates, additionalM)
% Dummy function for mediation.m. We don't need to check for missing data
% as vbm.m already does this...and this way we don't need to download all
% of the canlab helper functions.

nanvec = [];
x = x;
y=y;
m=m;
mediation_covariates = mediation_covariates;