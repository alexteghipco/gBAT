function [oSv,oSm,oR,oP,oPth,oPthStd,ccIds,pcSurvPos,pcSurvNeg,cclMxPPos,cclMxPNeg] = vbm(x,y,m,z,aType,nPerms,pType,vThr,cThr,oSize,mid,sv,nrm,ext,par,tail)
% Perform VBM analysis on brain data with permutation testing or FDR correction
%
% Inputs ----------------------------------------------------------------
% -x is a predictor variable (n x 1 where n is number of subjects)
%
% -y is the response variable (n x 1 where n is number of
%   subjects) you want to relate to x if performing mediation analysis. If
%   performing any correlation-based analysis this should be empty
%
% -m is the brain data (n x p where n is the number of subjects and p is
%   the number of voxels). If analysis is correlation-based, we treat m as
%   the y variable that we correlate with x
%
% -z is a set of variables to control for n x c where c is a nuissance
%   variable
%
% -aType: what type of analysis are we performing?
%   * 'mediation' -- mediation analysis
%   * 'pearson' -- pearson correlation coefficient (z should be empty)
%   * 'partial' -- partial correlation coefficient
%   * 'distance' -- nonlinear distance correlation (z should be empty)
%   * 'partialDistance -- nonlinear partial distance correlation
%
% -nPerms: number of permutation (0 for none)
%
% -pType: if set to 'fdr' will do FDR correction to p-values instead of
%   cluster correction
%
% -vThr: a p-value based threshold for the map
%
% -cThr: a p-value based threshold for cluster correction
%
% -oSize: 3 x 1 vector corresponding to size of *raw* maps from which m
%   comes from (i.e., 3D maps, not vectorized as expected by m). If
%   empty, we assume the inputs are ROIs rather than voxels so you can't
%   really do cluster correction (corresponding outputs will be empty as
%   well)
%
% -mid: which voxels were kept in the mask to vectorize m
%
% -sv: save out permutations incrementally...10 times over course of
%   permutation analysis
%
% -nrm: should we normalize? We zscore if this is set to true
%
% -ext: do you have any additional arguments to pass into the mediation or
%   correlation/partial correlation? Consider using 'robust' and/or 'boot'
%   for mediation. Consider using 'rows','pairwise' for correlation/partial
%   correlation. For example: {'rows'','pairwise'}. Leave empty if no
%   argument is to be passed in (i.e., [])
%
% -par: if true we use parfor
% 
% -tail: can be 'neg' to only test for negative effects, 'pos' for positive
%   effects or 'both' for two-tail. A one-tail test halves the p-value for
%   the corresponding effects while keeping the other p-values the same.
%   This is appropriate for eg FDR correction to handle the increased type
%   1 errors for one-tail tests. 
%
% Note, for mediation analysis we save out the whole model but only
% permutation test based on the indirect effect. If you want to test based
% on some other part of the model, change variable stt.
%
% Outputs ----------------------------------------------------------------
% Real data:
% -oSv : vectorized statistical map after voxelwise or cluster/FDR correction
% -oSm : 3D statistical map after voxelwise or cluster/FDR correction (for
%   mediation it is ab path t-statistic)
% -oR : if mediation t-stat for a b c' c ab, otherwise correlation
%   coefficient or distance correlation
% -oP : if mediation, p-vlaue for a b c' c ab, otherwise correlation
%   p-value
% -oPth : path for a b c' c ab
% -oPthStd : standardized path (effect size) for a b c' c ab
% -ccIds : voxel indices of all surviving clusters in oSv/oSm
% -pcSurvPos : proportion of clusters that were significant after cluster
%   correction (for positive effects)
% -pcSurvNeg : proportion of clusters that were significant after cluster
%   correction (for negative effects)
% -cclMxPPos: null distribution of cluster sizes (for positive effects)
% -cclMxPNeg: null distribution of cluster sizes (for negative effects)
%
% [oSv,oSm,oR,oP,oPth,oPthStd,ccIds,pcSurvPos,pcSurvNeg,cclMxPPos,cclMxPNeg] = medVoxFast(x,y,m,z,aType,nPerms,pType,vThr,cThr,oSize,mid,sv,nrm,ext,par)
%
% % Alex Teghipco // alex.teghipco@sc.edu
% % To fix: 
% % 1) matrix multiplication instead of for loop for correlation and
% partial correlation
% % 2) treating m as y in some cases is convenient but strange
% % 3) varargin

% setup some outputs that aren't always defined
oPth = [];  oPthStd = []; ccIds = []; pcSurvPos = []; pcSurvNeg=[]; cclMxPPos=[]; cclMxPNeg=[]; oSv = []; oSm = [];

% check that mediation.m exists...
tmp = which('mediation.m');
if isempty(tmp) && strcmpi(aType,'mediation')
    error('You have not downloaded mediation.m! Please see readme...')
elseif ~isempty(tmp) && strcmpi(aType,'mediation')
    warning('If you encounter errors, please manually ensure you have the Canlab-core repository downloaded and added to path...')
end

% set up rng
rng('default')
rng('shuffle');

% if mediation, which part of the model do we test?
cK = 5; % which column to keep from mediation results (t-stat). Let's say we only keep 5th column (indirect effect). This can be an array of indices reflecting which paths to keep
stt = 1; % Which of the kept columns from above should we keep during permutation testing? In other words, of the retained columns from the model, which column do we want to cluster correct for? This cannot be an array (overhead is too expensive for parfor)
lsp = 10; % how many times to save out permutations over course of analysis?

% If partial distance correlation, parellalization will be better if we
% implement it in another way so set to false (i.e., default
% parallelization is across voxels)
if strcmpi(aType,'partialDistance') && par % parfor will be faster at the level of computing correlation for partial distance so fix...
    parpd = true;
    par = false;
else
    parpd = false;
end

% save out progress as we go if necessary
if sv
    svN = linspace(1,nPerms,lsp);
    disp(['Temporary file with ongoing permutations will be saved in: ' pwd ' and is called medVoxTmpSv.mat'])
else
    svN = [];
end

%
if nPerms == 0 && ~isempty(cThr)
    cThr = [];
end

% number of voxels we loop over
tot = size(m,2);

% first remove nans...
id1 = find(isnan(x));
id2 = find(isnan(y));
if ~isempty(z)
    id3 = find(isnan(z));
else
    id3 = [];
end
id4 = find(any(isnan(m), 2));
id = unique([id1 id2 id3 id4]);
x(id) = [];
y(id) = [];
if ~isempty(z)
   z(id,:) = []; 
end
m(id,:) = [];

% normalize input data
if nrm
    x = zscore(x);
    y = zscore(y);
    z = zscore(z,0,2);
    m = zscore(m,0,2);
end

% concatenate all analysis arguments with z for mediation (nuissance covariates)
if isempty(ext)
    ext = {};
end
if ~isempty(z) && strcmpi(aType,'mediation')
    ext = [ext,{'covs',z}];
end

% start working on true data...waitbar setup first
if par 
    dq = parallel.pool.DataQueue;
    afterEach(dq, @updateWaitbar);
end
h = waitbar(0, 'Working on true data...'); 

if par
    counter = 0; % for parfor broadcast
    parfor i = 1:tot
        %tic
        try
            switch aType
                case 'mediation'
                    [~, stats] = mediation(x, y, m(:,i),ext{:});
                    if isfield(stats,'t')
                        oR(i,:) = stats.t(cK);
                    else
                        oR(i,:) = stats.z(cK);
                    end
                    oPth(i,:) = stats.paths(cK);
                    oPthStd(i,:) = stats.stdpaths(cK);
                    oP(i,:) = stats.p(cK);
                case 'partial'
                    [oR(i,:),oP(i,:)] = partialcorr(x,m(:,i),z,ext{:});
                case 'pearson'
                    [oR(i,:),oP(i,:)] = corr(x,m(:,i),ext{:});
                case 'distance'
                    nanRows = any(isnan(x), 2) | any(isnan(m(:,i)), 2); % need to ensure no NaN for distance correlation otherwise failure
                    xtmp = x;
                    ytmp = m(:,i);
                    xtmp(nanRows) = [];
                    ytmp(nanRows) = [];
                    [oR(i,:),oP(i,:), ~, ~] = bcdistcorr(xtmp, ytmp);
            end
        catch % if there is an error in statistical test (e.g., no variance will do this in some cases), we write in no effect. It may be better to use NaN here in the future but that may break bwconncomp so keeping as is for now
            if strcmpi(aType,'mediation')
                oR(i,:) = zeros(1,length(cK)); %[0 0 0 0 0];
                oPth(i,:) = zeros(1,length(cK)); %[0 0 0 0 0];
                oPthStd(i,:) = zeros(1,length(cK)); %[0 0 0 0 0];
                oP(i,:) = ones(1,length(cK)); %[1 1 1 1 1];
            else
                oR(i,:) = 0;
                oP(i,:) = 1;
            end
        end
        send(dq,i) % update waitbar for parfor (see bottom func)
        %toc
    end
else
    for i = 1:tot
        disp([num2str(i) ' of ' num2str(tot)])
        %tic 
        try
            switch aType
                case 'mediation'
                    [~, stats] = mediation(x, y, m(:,i),ext{:});
                    if isfield(stats,'t')
                        oR(i,:) = stats.t(cK);
                    else
                        oR(i,:) = stats.z(cK);
                    end
                    oPth(i,:) = stats.paths(cK);
                    oPthStd(i,:) = stats.stdpaths(cK);
                    oP(i,:) = stats.p(cK);
                case 'partial'
                    [oR(i,1),oP(i,1)] = partialcorr(x,m(:,i),z,ext{:});
                case 'pearson'
                    [oR(i,1),oP(i,1)] = corr(x,m(:,i),ext{:});
                case 'distance'
                    nanRows = any(isnan(x), 2) | any(isnan(m(:,i)), 2);
                    xtmp = x;
                    ytmp = m(i,:);
                    xtmp(nanRows) = [];
                    ytmp(nanRows) = [];
                    [oR(i,1),oP(i,1), ~, ~] = bcdistcorr(xtmp, ytmp);
                case 'partialDistance'
                    [oR(i,1),oP(i,1), ~, ~, ~] = pdcPerm(x,m(:,i),z,0,parpd);
            end
        catch
            if strcmpi(aType,'mediation')
                oR(i,:) = zeros(1,length(cK)); %[0 0 0 0 0];
                oPth(i,:) = zeros(1,length(cK)); %[0 0 0 0 0];
                oPthStd(i,:) = zeros(1,length(cK)); %[0 0 0 0 0];
                oP(i,:) = ones(1,length(cK)); %[1 1 1 1 1];
            else
                oR(i,:) = 0;
                oP(i,:) = 1;
            end
        end
        waitbar(i/tot, h);
        %toc
    end
end

% we save out whole model (all paths) but now which part do we do permutation testing
% for? We focus on indirect effect (see stt variable above) but you should consider
% looking at the other estimates to get more insight (e.g., is it a complete
% mediator? etc...)
if strcmpi(aType,'mediation')
    oRr = oR(:,stt);
    oPr = oP(:,stt);
    % oPthStdr = oPthStd(:,stt);
    % oPthr = oPth(:,stt);
else
    oRr = oR;
    oPr = oP;
end

% The next thing we have to do is separate out positive and negative
% effects. Since we are cluster correcting, it's possible a large cluster
% displays adjacent positive and negative effects that get grouped
% together, which we would treat as distinct when we interpret the map...
[oRNeg,oRPos,oPNeg,oPPos] = sepPosNeg(oRr,oPr); % oSmPos

% We can now apply our voxelwise threshold (p-value based) and identify all
% clusters as well as their sizes for the two maps
[oRPos,oPPos,cclPos,~] = voxThresh(oRPos,oPPos,lower(pType),vThr,oSize,mid,true,tail);
[oRNeg,oPNeg,cclNeg,~] = voxThresh(oRNeg,oPNeg,lower(pType),vThr,oSize,mid,true,tail);

% If we FDR correct the p-values we won't cluster correct as well...
if ~strcmpi(pType,'fdr')
    % If cluster correction, start permutation analysis
    for n = 1:nPerms
        disp(['Working on permutation: ' num2str(n)])
        waitbar(n/nPerms, h, ['Running permutation analysis: ' num2str(n)]); % update same waitbar as before
        xp = x(randperm(length(x))); % permute our predictor
        if par % just as above but run for permuted data...should make this a function in future...
            parfor i = 1:tot
                try
                    switch aType
                        case 'mediation'
                            [~, stats] = mediation(xp, y, m(:,i),ext{:});
                            if isfield(stats,'t')
                                oRP(i,1) = stats.t(stt);
                            else
                                oRP(i,1) = stats.z(stt);
                            end
                            oPP(i,1) = stats.p(stt);
                        case 'partial'
                            [oRP(i,1),oPP(i,1)] = partialcorr(xp,m(:,i),z,ext{:});
                        case 'pearson'
                            [oRP(i,1),oPP(i,1)] = corr(xp,m(:,i),ext{:});
                        case 'distance'
                            nanRows = any(isnan(xp), 2) | any(isnan(m(:,i)), 2);
                            xtmp = xp;
                            ytmp = m(:,i);
                            xtmp(nanRows) = [];
                            ytmp(nanRows) = [];
                            [oRP(i,1),oPP(i,1), ~, ~] = bcdistcorr(xtmp, ytmp);
                    end
                catch
                    oRP(i,1) = 0;
                    oPP(i,1) = 1;
                end
            end
        else
            for i = 1:tot
                try
                    switch aType
                        case 'mediation'
                            [~, stats] = mediation(xp, y, m(:,i),ext{:});
                            if isfield(stats,'t')
                                oRP(i,1) = stats.t(stt);
                            else
                                oRP(i,1) = stats.z(stt);
                            end
                            oPP(i,1) = stats.p(stt);
                        case 'partial'
                            [oRP(i,1),oPP(i,1)] = partialcorr(xp,m(:,i),z,ext{:});
                        case 'pearson'
                            [oRP(i,1),oPP(i,1)] = corr(xp,m(:,i),ext{:});
                        case 'distance'
                            nanRows = any(isnan(xp), 2) | any(isnan(m(:,i)), 2);
                            xtmp = xp;
                            ytmp = m(:,i);
                            xtmp(nanRows) = [];
                            ytmp(nanRows) = [];
                            [oRP(i,1),oPP(i,1), ~, ~] = bcdistcorr(xtmp, ytmp);
                        case 'partialDistance'
                            [oRP(i,1),oPP(i,1), ~, ~, ~] = pdcPerm(xp,m(:,i),z,0,parpd);
                    end
                catch
                    oRP(i,1) = 0;
                    oPP(i,1) = 1;
                end
            end
        end

        % Just as above for the true data we separate the maps and record
        % maximum cluster sizes
        [oRPNeg,oRPPos,oPPNeg,oPPPos] = sepPosNeg(oRP,oPP);
        [~,~,~,cclMxPPos(n,1)] = voxThresh(oRPPos,oPPPos,lower(pType),vThr,oSize,mid,false,tail); 
        [~,~,~,cclMxPNeg(n,1)] = voxThresh(oRPNeg,oPPNeg,lower(pType),vThr,oSize,mid,false,tail); 

        % if permutation index overlaps with our set of 10 linearly spaced
        % values across whole permutation range, save perms temporarily
        if ~isempty(intersect(n,svN))
            try
                save('medVoxTmpSv.mat')
            catch
                save('medVoxTmpSv.mat','-v7.3')
            end
        end
    end
    delete(h) % don't need waitbar any more

    % Final analysis step--we need to apply the cluster-level threshold based on the
    % permutation analysis...the appropriate percentile for permuted map
    % cluster sizes is used to identify significant clusters in the true map
    if ~isempty(cThr)
        %[oRNegv,oRPosv,oPNegv,oPPosv] = sepPosNeg(oRr,oPr); % oSmPos
        [oSvPos,oSmPos,~,pcSurvPos] = clustThresh(oRPos(mid),cclPos,cclMxPPos,cThr,mid);
        [oSvNeg,oSmNeg,~,pcSurvNeg] = clustThresh(oRNeg(mid),cclNeg,cclMxPNeg,cThr,mid);
    else
        oSmPos = oRPos;
        oSmNeg = oRNeg;
        if ~isempty(mid)
            oSvPos = oRPos(mid);
            oSvNeg = oRNeg(mid);
        else
            oSvPos = oRPos;
            oSvNeg = oRNeg;
        end
    end

    % let's also combine our maps...essentially sum them but check that there
    % is no overlap between the maps (sanity check, there should not be...)
    oSv = combPosNeg(oSvNeg,oSvPos);
    oSm = combPosNeg(oSmNeg,oSmPos);
    if ~isempty(oSize)
        cc = bwconncomp(oSm);
        ccIds = cc.PixelIdxList;
    end
else
    if ~isempty(oSize)
        warning('Your selected options are not debugged yet...')
        oSv = combPosNeg(oRNeg,oRPos);
        oPv = combPosNeg(oPNeg,oPPos);
    end
end

% function for updating waitbar with counter...
    function updateWaitbar(~)
        counter = counter + 1;
        waitbar(counter / tot, h);
    end
end
