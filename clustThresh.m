function [oSv,oSm,ncSurv,pcSurv] = clustThresh(s,cl,cclMxP,cThr,mid)
% Thresholds a voxelwise map based on cluster size of permuted maps. A future
% update will permit permutation analysis at the voxel level without doing
% this kind of cluster correction (i.e., here we apply a statistical
% threshold and look at max cluster size to correct for cluster size).
% Inputs ----------------------------------------------------------------
% -s is n x 1 where n is number of samples and col is a statistical map
% (thresholded or unthresholded)
% -cThr is the threshold to apply for voxel correction (e.g., 0.05 means
%   cluster correction at p < 0.05)
% -cl is a cell array of clusters in s. These should come from a
%   *thresholded* version of s (i.e., voxelwise threshold not cluster
%   forming threshold)
% -cclMxP is a cell array of maximum cluster sizes in permuted data. These
%   sizes should be computed on permuted data that has the same threshold
%   applied to it as was applied to s in order to get cl
% -mid: which voxels were kept in the mask to vectorize s and p
%
% Outputs ----------------------------------------------------------------
% -oSv is the statistical map corrected for cluster size based on cluster
%   sizes in permuted data. Data is vectorized in line with input s so rows
%   match with s.
% -oSm is the same map as oSv but in 3D format consistent with the shape of
%   the original images before thresholding
% -ncSurv is the number of clusters above the cluster-forming
%   threshold specified
% -pCSurv is the proportion of clusters above the cluster-forming
%   threshold specified
%
% [oSv,oSm] = clustThresh(s,cl,cclMxP,oSize,mid)

% extract clusters for true data
ccs = cl{1}.PixelIdxList;
ccsSz = cellfun(@length, ccs);

% get threshold based on p-value (convert to percentile)
cThr = prctile(cclMxP,100-(cThr*100));

% get surviving clusters
surv = find(ccsSz >= cThr);
ncSurv = length(surv);
pcSurv = ncSurv./length(ccsSz);

% setup outputs
oSm = zeros(cl{1}.ImageSize);
if ~isempty(surv)
    [c,~,~] = intersect(ccs{surv},mid);
    oSm(c) = s(c);
else
    warning('No clusters survived correction!')
end
oSv = oSm(mid);

end



