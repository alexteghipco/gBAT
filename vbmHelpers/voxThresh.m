function [oS,oP,ccl,cclMx,ccKp] = voxThresh(s,p,pType,thr,oSize,mid,ccKp,tail)
% Threshold voxelwise maps and extract their clusters as well as cluster
% sizes
% Inputs ----------------------------------------------------------------
% -s is n x m where n is number of samples and m are statistical maps to
%   treat independently (e.g., statistical values at each voxel)
% -p is n x m where n is number of samples and m are statistical maps to
%   treat independently (p-values associated with statistical values in s)
% -pType : can be 'fdr' to fdr correct input p-values, anything else does
%   nothing
% -thr: threshold for p-values (fdr or normal depending on pType)
% -oSize:  x 1 vector corresponding to size of *raw* maps from which s
%   comes from (i.e., 3D maps, not vectorized as expected by s and p). If
%   empty, we assume the inputs are ROIs rather than voxels so you can't
%   really do cluster correction (corresponding outputs will be empty as
%   well)
% -mid: which voxels were kept in the mask to vectorize s and p 
% -ccKp: do we keep the cluster identities for each map in s? This will
%   take up a lot of space and they are not essential to keep for
%   permutation testing but maybe you have some reason to...set to true or
%   false
%
% Outputs ----------------------------------------------------------------
% -oS is the thresholded image(s). If ROI, will be ROI vector, otherwise
%   image 
% -oP is the p-values. These are maybe fdr-corrected depending on
%   your setting for pType. If ROI will be ROI vector, otherwise image
% -ccl is a cell array of cluster voxel identities for each image. The
%   number of clusters in each image can vary
% -cclMx is the maximum cluster size in each image
%
% [oS,oP,ccl,cclMx] = voxThresh(s,p,pType,thr,oSize,mid)
%
% Alex Teghipco // alex.teghipco@sc.edu
ccl = []; cclMx = []; cclKp = []; % if ROIs we won't do cluster identification....

if strcmpi(tail,'pos')
    id = find(s > 0);
    p(id) = p(id)./2;
elseif strcmpi(tail,'neg')
    id = find(s < 0);
    p(id) = p(id)./2;
end

switch pType
    case 'fdr'
        for i = 1:size(p,2)
            p(:,i) = mafdr(p(:,i),'BHFDR',true); % see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3699340/
        end
end

if ~isempty(oSize)
    oS = zeros([oSize size(p,2)]);
    oP = zeros([oSize size(p,2)]);
    for i = 1:size(p,2)
        tmp = zeros(oSize);
        tmp(mid) = s(:,i);
        oS(:,:,:,i) = tmp;
        tmp = zeros(oSize);
        tmp(mid) = p(:,i);
        oP(:,:,:,i) = tmp;
    end

    bmt = oP > thr;
    oS(bmt) = 0; % thresholded stats map
    id = oS ~= 0;
    tmp = zeros(size(oS));
    tmp(id) = 1; % binarized version

    for i = 1:size(p,2)
        cc = bwconncomp(squeeze(tmp(:,:,:,i)));
        if ccKp
            ccl{i} = cc;
        else
            ccl = [];
        end
        if cc.NumObjects > 0
            cclMx(i,1) = max(cellfun(@length, cc.PixelIdxList));
        else
            cclMx(i,1) = 0;
            warning('Voxelwise threshold resulted in an empty map')
        end
    end
else
    oS = s;
    oP = p;

    id = find(oP >= thr);
    oS(id) = 0;

    ccl = [];
    cclMx = [];
end