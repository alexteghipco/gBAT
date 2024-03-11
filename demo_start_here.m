%% Example 1: FDR corrected 1-tail mediation analysis
%% Extract behavioral data
beh = readtable('./abcBeh.csv');
subs = beh.Subject;
moca  = beh.TotalScore;
age = beh.Age;

%% Find matching imaging data and load...
tor = [];
for i = 1:length(subs)
    disp(num2str(i))
    try
        f = ls([pwd '\mris\smwp1' subs{i} '.nii.gz']);
        f = f(1:end);
        tmp = load_nifti([pwd '\mris\' f]);
        gm(:,:,i) = tmp.vol(:);
    catch
        gm(:,:,i) = 0;
        tor = [tor; i];
    end
end

tor2 = [];
for i = 1:length(subs)
    disp(num2str(i))
    try
        f = ls([pwd '\mris\smwp2' subs{i} '.nii.gz']);
        f = f(1:end);
        tmp = load_nifti([pwd '\mris\' f]);
        wm(:,:,i) = tmp.vol(:);
    catch
        wm(:,:,i) = 0;
        tor2 = [tor2; i];
    end
end

% sanity check for different missing brain data...
isequal(tor,tor2) % if the result of this is not 1 then you need to treat tor and tor2 separately and make two sets of behavioral variables, 1 for GM and 1 for WM...

gm = squeeze(gm)';
wm = squeeze(wm)';
gm(tor,:) = [];
wm(tor,:) = [];
subs(tor) = [];
age(tor) = [];
moca(tor) = [];
beh(tor,:) = [];

% get size of image so we can reverse image indexing
sz = size(tmp.vol);
 
% create a mask for removing some voxels (separate masks for GM and WM to
% make the analysis run faster)
s = sum(gm)./size(gm,1); % proportion of subs with value in voxel
gmid = find(s > 0.1); % you can vary this number and check your data for a well fitting mask
tmp.vol = zeros(size(tmp.vol));
tmp.vol(gmid) = 1;
save_nifti(tmp,'mask_gm.nii.gz')
s = sum(wm)./size(wm,1);
wmid = find(s > 0.1); 
tmp.vol = zeros(size(tmp.vol));
tmp.vol(wmid) = 1;
save_nifti(tmp,'mask_wm.nii.gz')


% MORE SANITY
[ctmp,iatmp,ibtmp] = intersect(tmpp.subs(tmpp.kid),app.inBeh.Subject);
corr(tmpp.age(tmpp.kid(iatmp)),app.inBeh.Age(ibtmp))
corr(tmpp.moca(tmpp.kid(iatmp)),app.inBeh.TotalScore(ibtmp))
corr(tmpp.tiv(tmpp.kid(iatmp)),app.inBeh.TIV(ibtmp)) % FIX TIV IN abcBeh.csv

for ttt = 1:length(iatmp)
    rtmp(ttt,1) = corr(tmpp.gm(tmpp.kid(iatmp(ttt)),:)',app.inBrain(ibtmp(ttt),:)'); % FIX TIV IN abcBeh.csv
end

for i = 1:length(insub)
    id = find(ismember(tivbeh.Var1,insub{i}));
    if isempty(id)
        tivfin(i,:) = NaN;
    else
        tivfin(i,:) = tivbeh.CSF_Total(id);
    end
end


% sanity check "significant" relationships by correlating variables to
% brain data...moca and age separately
[r,p] = corr(gm(:,gmid),moca);
[r2,p2] = corr(wm(:,wmid),moca);

length(find(mafdr(p,'BHFDR',true)))./length(gmid)
length(find(mafdr(p2,'BHFDR',true)))./length(wmid)

[r,p] = corr(gm(:,gmid),age);
[r2,p2] = corr(wm(:,wmid),age);
length(find(mafdr(p,'BHFDR',true)))./length(gmid)
length(find(mafdr(p2,'BHFDR',true)))./length(wmid)

% Extract total intracranial volume (nuissance variable)
tivd = readtable('TIVDataABC.xlsx');
[c,ia,ib] = intersect(tivd.Var1,subs);
subs = subs(ib);
beh = beh(ib,:);
moca = moca(ib);
age = age(ib);
wm = wm(ib,:);
gm = gm(ib,:);
tiv = tivd.TIV(ia);
wmv = tivd.WM_Total(ia);
gmv = tivd.GM_Total(ia);
csfv = tivd.CSF_Total(ia);
%cnc(1:131) = 0;
%cnc(132:length(subs)) = 1;

% save our inputs in case we need to reload
save('input_age_moca.mat','-v7.3')


%[ctmp,iatmp,ibtmp] = intersect(tmp.beh.Subject(tmp.kid),app.inBeh);

% Subjects 133 through 243 tested COVID positive and the relationship
% between age and moca is weaker in this group...we will ignore it for now.
id = 1:132; % this is ABC
%[r,p] = corr(moca(id),age(id))


%% Voxelwise mediation analysis example...
% The standard options for mediation.m will probably be most reasonable to
% use. Adding robust or bootstrapped p-value options to mediation
% dramatically increases run time (by days for this data on a top of the
% line CPU). Note also, you can try parfor but it is often not more
% efficient if we are working with many voxels...by the way, we will start
% by analyzing white matter morphometry, then move on to grey matter
[oSv_wm_med,oSm_wm_med,oR_wm_med,oP_wm_med,~,~,~,~,~,cclMxPPos_wm_med,cclMxPNeg_wm_med] = vbm(age(id),moca(id),wm(id,wmid),[tiv(id)],'mediation',0,'default',0.05,[],sz,wmid,false,false,[],false,'neg'); % [{'robust','boot'}] <-- example of args input if we wanted to do robust OLS with bootstrapped p-values for mediation

% you can ask vbm.m to output FDR corrected p-values but since we wanted
% access to both FDR corrected and uncorrected values, we do this manually
% without specifying FDR correction to vbm...

% as you can see if we run the below lines...nothing survives correction
% for multiple comparisons:(
min(mafdr(oP_wm_med,'BHFDR',true)) 
length(find(oP_wm_med < 0.05)) 

% If we run the same analysis in the grey matter, we see the same result
[oSv_gm_med,oSm_gm_med,oR_gm_med,oP_gm_med,~,~,~,~,~,cclMxPPos_gm_med,cclMxPNeg_gm_med] = vbm(age,moca,gm(:,gmid),[tiv cnc'],'mediation',0,'default',0.05,[],sz,gmid,false,false,[],false,'neg');
min(mafdr(oP_gm_med,'BHFDR',true)) % 0.17 if no 
length(find(oP_gm_med < 0.05)) % 0.7749

%% ROI-based mediation analysis example...
% First, we need to 'bin' our morphometry data by ROIs.
% load in ROI names in atlas and their corresponding values in the atlas
labs = readtable('jhu.txt');
% load in the atlas
roiin = load_nifti('rjhu.nii');

% loop over unique values in the atlas and extract voxel identities for
% each
for i = 1:height(labs)
    rois{i} = find(roiin.vol == labs.Var1(i));
end

% now get mean morphometry within each ROI. Because our atlas has GM and WM
% regions, if the ROI is in GM, we pull from the gm morphometry data and
% same for WM. We ignore CSF (e.g., ventricles). 
kid = 1:131 %1:133; % this is ABC folks (i.e., non-covid positive cohort)
for i = 1:height(labs) % loop over ROIs
    disp(num2str(i))
    if labs.Var4(i) == 1 % if GM ROI
        [c,ia,ib] = intersect(rois{i},gmid); % find grey matter voxel indices that overlap with the roi
        tmp = gm(kid,ib); % extract ABC folks' corresponding overlapping voxels
       % zc = find(all(tmp == 0, 1)); % remove any voxels that are 0 across subjects (i.e., ROIs not overlapping with the voxelwise maps for example because they intrude on another tissue type)
     %   tmp(:,zc) = []; % remove those voxels
        morph(:,i) = mean(tmp,2,'omitnan'); % get mean across voxels for the ROI-level morphometry, for each participant
    elseif labs.Var4(i) == 2 % if WM ROI
        [c,ia,ib] = intersect(rois{i},wmid); % same as above...
        tmp = wm(kid,ib);
       % zc = find(all(tmp == 0, 1));
      %  tmp(:,zc) = [];
        morph(:,i) = mean(tmp,2,'omitnan');
    end
    %rois{i} = find(roiin.vol == labs.Var1(i));
end

% try something else--no removal...
tmp = zeros(sz);
tmp(gmid) = 

% remove bad ROIs--for example, ROIs with a mean of 0 across all subjects
% (these are CSF ROIs, which our if statements ignore) or ROIs with NaNs
% (some ROIs show no voxel overlap with our map at all which our conditions
% above do not explicitly exclude)
tor = find(sum(morph) == 0);
morph(:,tor) = [];
labs(tor,:) = [];
rois(tor) = [];
tor = find(all(isnan(morph), 1));
morph(:,tor) = [];
labs(tor,:) = [];
rois(tor) = [];

% Run ROI-level analysis. NOW we can do the fancy modeling--mediation with
% robust OLS regression and bootstrapped p-values. 
 [oSv,oSm,oR,oRm,oP,oPc,oPm,oPcm,oPth,oPall,oPthStd,oPthm,oPallm,oPthStdm,ccIds,pcSurvPos,pcSurvNeg,cclMxPPos,cclMxPNeg] = vbm(age(kid),moca(kid),morph,[tiv(kid)],'mediation',0,'default',0.05,[],[],[],false,false,[{'boot','robust'}],false,'both',false,0,false)
%[oSv_wm_med,oSm_wm_med,oR_wm_med,oP_wm_med,~,~,~,~,~,cclMxPPos_wm_med,cclMxPNeg_wm_med] = vbm(age(kid),moca(kid),morph,[tiv(kid)],'mediation',0,'default',0.05,[],[],[],false,false,[{'boot','robust'}],false,'both');  % [tiv wmv gmv csfv cnc']
%[oSv_wm_med,oSm_wm_med,oR_wm_med,oP_wm_med,~,~,~,~,~,cclMxPPos_wm_med,cclMxPNeg_wm_med] = vbm(age(kid),moca(kid),morph,[tiv(kid)],'mediation',0,'default',0.05,[],[],[],false,false,[{'boot','robust'}],false,'both',true,1000);  % [tiv wmv gmv csfv cnc']

% Again, we manually FDR-correct values so we have both uncorrected and
% corrected values available. Also apply an FDR threshold 
oP_wm_med_fdr = mafdr(oP_wm_med,'BHFDR',true);
id = find(oP_wm_med_fdr > 0.001);
oR_wm_med_fdr = oR_wm_med;
oR_wm_med_fdr(id) = 0;
id = find(isnan(oR_wm_med_fdr));
oR_wm_med_fdr(id) = 0;
save('tmp.mat','-v7.3')

% Finally, we can save the ROI effects as a nifti!
template = load_nifti('/home/alex/Downloads/jhu_resampled_to_FSL_fixed.nii'); % our atlas but in standard MNI152 space instead of warped to subject template space to aid visualization
otemplate = load_nifti('/home/alex/Downloads/brainSurfer-main/brainMapsforTesting/multidim/phonological.nii.gz_specificity_z.nii.gz'); % random file to overwrite so that we know what nifti structure of output file should look like
otemplate.vol = zeros(size(otemplate.vol));
for i = 1:length(oR_wm_med_fdr)
    vid = find(template.vol == labs.Var1(i));
    otemplate.vol(vid) = oR_wm_med_fdr(i);
end
save_nifti(otemplate,'Age_Moca_Morph_ROIs.nii.gz')









