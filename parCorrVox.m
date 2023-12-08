function [oR,oP,oRP,oPP] = parCorrVox(x,y,z,nPerms,par,corrType)
% Perform correlation, partial correlation, distance correlation or partial
% distance correlation on brain data with permutation testing for mult.
% comp. corr. (correcting for cluster size).
%
% NOTE: we remove NaNs automatically for partial distance correlation but
% you will get NaN outputs if you have no variance in any of the cleaned
% variables.
%
% Inputs ----------------------------------------------------------------
% -x is a variable you will relate to the the brain (n x 1 where n is number
%   of subjects)
% -y is brain to loop over n x p where n is subjects and p is voxels
% -z is a set of variables to control for n x c where c is a nuissance
%   variable
% -nPerms: number of permutations
% -par: if true we use parfor
% -corrType: what kind of correlation analysis are we performing?
%   * 'pearson' -- pearson correlation coefficient (z should be empty)
%   * 'partial' -- partial correlation coefficient
%   * 'distance' -- nonlinear distance correlation (z should be empty)
%   * 'partialDistance -- nonlinear partial distance correlation
%
% Outputs ----------------------------------------------------------------
% Real data:
% -oR : Correlation coefficient between x and y controlling for z
% -oP : p-vlaue
%
% Permuted data:
% -oRP : Correlation coefficient between x and y controlling for z
% -oPP : p-vlaue
%
% [oR,oP,oRP,oPP] = parCorrVox(x,y,z,nPerms,par,corrType)
%
% Alex Teghipco // alex.teghipco@sc.edu

if strcmpi(corrType,'partialDistance') && par % parfor will be faster at the level of computing correlation for partial distance so fix...
    parpd = true;
    par = false;
else
    parpd = false;
end

% work on true data...
if par
    parfor i = 1:size(y,1)
        disp([num2str(i)  ' of ' num2str(size(y,1))])
        switch corrType
            case 'partial'
                [oR(i,1),oP(i,1)] = partialcorr(x,y(i,:)',z,'rows','pairwise');
            case 'pearson'
                [oR(i,1),oP(i,1)] = corr(x,y(i,:)','rows','pairwise');
            case 'distance'
                nanRows = any(isnan(x), 2) | any(isnan(y(i,:)'), 2);
                xtmp = x;
                ytmp = y(i,:)';
                xtmp(nanRows) = [];
                ytmp(nanRows) = [];
                [oR(i,1),oP(i,1), ~, ~] = bcdistcorr(xtmp, ytmp);
        end
    end
else
    for i = 1:size(y,1)
        disp([num2str(i)  ' of ' num2str(size(y,1))])
        switch corrType
            case 'partial'
                [oR(i,1),oP(i,1)] = partialcorr(x,y(i,:)',z,'rows','pairwise');
            case 'pearson'
                [oR(i,1),oP(i,1)] = corr(x,y(i,:)','rows','pairwise');
            case 'distance'
                nanRows = any(isnan(x), 2) | any(isnan(y(i,:)'), 2);
                xtmp = x;
                ytmp = y(i,:)';
                xtmp(nanRows) = [];
                ytmp(nanRows) = [];
                [oR(i,1),oP(i,1), ~, ~] = bcdistcorr(xtmp, ytmp);
            case 'partialDistance'
                [oR(i,1),oP(i,1), ~, ~, ~] = pdcPerm(x,y(i,:)',z,0,parpd);
        end
    end
end

% start permutations
rng('shuffle');
for n = 1:nPerms
    disp(['Permutation: ' num2str(n) ' of ' num2str(nPerms)])
    xp = x(randperm(length(x)));
    if par
        parfor i = 1:size(y,1)
            switch corrType
                case 'partial'
                    [oRP(i,n),oPP(i,n)] = partialcorr(xp,y(i,:)',z,'rows','pairwise');
                case 'pearson'
                    [oRP(i,n),oPP(i,n)] = corr(xp,y(i,:)','rows','pairwise');
                case 'distance'
                    nanRows = any(isnan(xp), 2) | any(isnan(y(i,:)'), 2);
                    xtmp = xp;
                    ytmp = y(i,:)';
                    xtmp(nanRows) = [];
                    ytmp(nanRows) = [];
                    [oRP(i,n),oPP(i,n), ~, ~] = bcdistcorr(xtmp, ytmp);
            end
        end
    else
        for i = 1:size(y,1)
            switch corrType
                case 'partial'
                    [oRP(i,n),oPP(i,n)] = partialcorr(xp,y(i,:)',z,'rows','pairwise');
                case 'pearson'
                    [oRP(i,n),oPP(i,n)] = corr(xp,y(i,:)','rows','pairwise');
                case 'distance'
                    nanRows = any(isnan(xp), 2) | any(isnan(y(i,:)'), 2);
                    xtmp = xp;
                    ytmp = y(i,:)';
                    xtmp(nanRows) = [];
                    ytmp(nanRows) = [];
                    [oRP(i,n),oPP(i,n), ~, ~] = bcdistcorr(xtmp, ytmp);
                case 'partialDistance'
                    [oRP(i,n),oPP(i,n), ~, ~, ~] = pdcPerm(xp,y(i,:)',z,0,parpd);
            end
        end
    end
end
