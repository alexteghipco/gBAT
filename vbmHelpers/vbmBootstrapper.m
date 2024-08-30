function [oR,oP] = vbmBootstrapper(x,y,z,m,index,aType,cK,ext)
% Resample x, m, and z using the current bootstrap sample indices
xSample = x(index);
if ~isempty(y)
    ySample = y(index);
end
mSample = m(index, :);
if ~isempty(z)
    if size(z,2) > 1
        zSample = z(index, :);
    else
        zSample = z;
    end
end

% Calculate the partial correlation for the bootstrap sample
if strcmpi(aType,'partial')
    [oR,oP] = partialcorr(xSample, mSample, zSample, ext{:});
elseif strcmpi(aType,'pearson')
    [oR,oP] = corr(xSample, mSample,ext{:});
elseif strcmpi(aType,'distance')
    nanRows = any(isnan(xSample), 2) | any(isnan(mSample), 2); % need to ensure no NaN for distance correlation otherwise failure
    xtmp = xSample;
    ytmp = mSample;
    xtmp(nanRows) = [];
    ytmp(nanRows) = [];
    [oR,oP, ~, ~] = bcdistcorr(xSample, mSample);
elseif strcmpi(aType,'partialDistance')
    [oR,oP, ~, ~, ~] = pdcPerm(xSample,mSample,zSample,0,false);
elseif strcmpi(aType,'mediation')
    [~, stats] = mediation(xSample, ySample, mSample,ext{:});
    % if medCoeff
        %oR = stats.paths(cK);
    oR = stats.mean;
    % else
    %     if isfield(stats,'t')
    %         oR = stats.t(cK);
    %     else
    %         oR = stats.z(cK);
    %     end
    %end
    %oP = stats.p(cK);
    oP = stats.p;
end
