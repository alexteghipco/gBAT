function [oR,oP,oPth,oPthStd,oRP,oPP,oPthP,oPthStdP] = medVox(x,y,m,z,nPerms,par)
% Perform mediation analysis on brain data with permutation testing
% Inputs ----------------------------------------------------------------
% -x is a variable (n x 1 where n is number of subjects)
% -y is the variable you want to relate to x (n x 1 where n is number of subjects)
% -m is the mediator (we assume it's the brain and is p x n where n is the
%   number of subjects and p is the number of voxels)
% -z is a set of variables to control for n x c where c is a nuissance
%   variable
% -nPerms: number of permutation
% -par: if true we use parfor
%
% Outputs ----------------------------------------------------------------
% Real data:
% -oR : t-stat for a b c' c ab
% -oP : p-vlaue for a b c' c ab
% -oPth : path for a b c' c ab
% -oPthStd : standardized path (effect size) for a b c' c ab
%
% Permuted data:
% -oRP : t-stat for a b c' c ab
% -oPP : p-vlaue for a b c' c ab
% -oPthP : path for a b c' c ab
% -oPthStdP : standardized path (effect size) for a b c' c ab
%
% [oR,oP,oPth,oPthStd,oRP,oPP,oPthP,oPthStdP] = medVox(x,y,m,z,nPerms,par)
%
% % Alex Teghipco // alex.teghipco@sc.edu

% work on true data...
if par
    parfor i = 1:size(m,1)
        disp([num2str(i)  ' of ' num2str(size(m,1))])
        try
            %tmp = [x y m(i,:)' z];
            [~, stats] = mediation(x, y, m(i,:)','covs',z);
            oR(i,:) = stats.t;
            oPth(i,:) = stats.paths;
            oPthStd(i,:) = stats.stdpaths;
            oP(i,:) = stats.p;
        catch
            oR(i,:) = [0 0 0 0 0];
            oPth(i,:) = [0 0 0 0 0];
            oPthStd(i,:) = [0 0 0 0 0];
            oP(i,:) = [1 1 1 1 1];
        end
    end
else
    for i = 1:size(m,1)
        disp([num2str(i)  ' of ' num2str(size(m,1))])
        try
            [~, stats] = mediation(x, y, m(i,:)','covs',z);
            oR(i,:) = stats.t;
            oPth(i,:) = stats.paths;
            oPthStd(i,:) = stats.stdpaths;
            oP(i,:) = stats.p;
        catch
            oR(i,:) = [0 0 0 0 0];
            oPth(i,:) = [0 0 0 0 0];
            oPthStd(i,:) = [0 0 0 0 0];
            oP(i,:) = [1 1 1 1 1];
        end
    end
end

% start permutations
rng('shuffle');
for n = 1:nPerms
    disp(['Permutation: ' num2str(n) ' of ' num2str(nPerms)])
    xp = x(randperm(length(x)));
    if par
        parfor i = 1:size(m,1)
            try
                [~, stats] = mediation(xp, y, m(i,:)','covs',z);
                oRP(i,:,n) = stats.t;
                oPthP(i,:,n) = stats.paths;
                oPthStdP(i,:,n) = stats.stdpaths;
                oPP(i,:,n) = stats.p;
            catch
                oRP(i,:,n) = [0 0 0 0 0];
                oPthP(i,:,n) = [0 0 0 0 0];
                oPthStdP(i,:,n) = [0 0 0 0 0];
                oPP(i,:,n) = [1 1 1 1 1];
            end
        end
    else
        for i = 1:size(m,1)
            try
                [~, stats] = mediation(xp, y, m(i,:)','covs',z);
                oRP(i,:,n) = stats.t;
                oPthP(i,:,n) = stats.paths;
                oPthStdP(i,:,n) = stats.stdpaths;
                oPP(i,:,n) = stats.p;
            catch
                oRP(i,:,n) = [0 0 0 0 0];
                oPthP(i,:,n) = [0 0 0 0 0];
                oPthStdP(i,:,n) = [0 0 0 0 0];
                oPP(i,:,n) = [1 1 1 1 1];
            end
        end
    end
end
