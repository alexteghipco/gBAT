%% Demo for VBM using partial correlation or mediation analysis with permutation testing for cluster correction

%% Dependencies
% 1.Please download the following packages/data and make sure they are in
%   your MATLAB path: 1. https://github.com/canlab/MediationToolbox <--- has
%   the mediation.m function we will use from Tor's lab
%
% 2. https://github.com/canlab/CanlabCore <--- silly, but they make you
%   download their whole giant lab repo for some basic helper functions for
%   mediation.m
%
% 3. Check that you have the statistics and machine learning matlab toolbox
%   installed (type ver into the command window, hit enter, and search for the
%   toolbox)
%
% 4. If you plan to use distance correlation or partial distance
%   correlation for VBM, download this package and make sure it is in your path: https://github.com/alexteghipco/partialDistanceCorrelation
%
% 5. Download toy data for this demo here and place it in the directory you
%   will be working out of (you need a USC email account):
%   https://emailsc-my.sharepoint.com/:f:/g/personal/alex_teghipco_sc_edu/EgWF8kxZPdJMn5RhKaHZmBsBg5slH7g7yaPhsNldI8kJ3w?e=PfGWDg
%
% 6. If you are uncomfortable manipulating nifti data in matlab and are
%   going to exactly follow the commands below for loading in
%   and saving nifti data you'll need load_nifti and save_nifti. You can
%   use freesurfer's version of these if you have them handy or download
%   the versions that come with brainSurfer (same code edited to work
%   better with gzipped files esp. on windows): https://github.com/alexteghipco/brainSurfer/tree/main/scripts/fs
%
%   Note, you can also use brainLoad.m to load in brain data if you have
%   e.g., have freesurfer or spm or the matlab image processing toolbox already
%   installed

%% Importing data
% load in our behavioral data
beh = readtable('behaviorConcat.xlsx');

% load in our imaging data--the attached file is cat12 gm and wm volume
% maps concatenated in one nifti (see brainLoad.m if you need help loading
% in individual niftis you have generated for each participant)
all = load_nifti('sall3mm_RH.nii.gz');
% Note, here we have a 4D nifti image where the 4th dimension is subjects
% but you can just load in your own nifti files however you want. The idea
% is that the vectorization loop below assumes your data is 4D with the 4th
% dim being subjects. Here too, however, you can vectorize the data however
% is easiest for you. Ultimately the scripts we use will need: 1) the size
% of the original images, 2) indices of the voxels that were submitted to
% the analysis (*not* subscripts), and 3) a 2D matrix of voxels by subjects
% that we perofrm analysis on. I've provided a script called brainLoad.m,
% which will try to get all of this info/variables together for you by
% having you select a folder that contains all of your subject-level
% niftis. It will look for either spm, matlab's imaging toolbox, or
% freesurfer's load_nifti to load in brain files...for a demonstration of
% how to use brainLoad.m see our live code Tutorial3.mlx in this
% repository: https://github.com/alexteghipco/StabilitySelection

% vectorize imaging data
for i = 1:size(all.vol,4)
    tmp = all.vol(:,:,1:56,i); % this indexing is peculiar to the file I loaded in, which concatenates two images along the z-dimension that we unpack here. For a normal file where 4th dim is subjects you can change this line to tmp = all.vol;
    gm(:,i) = tmp(:);
    tmp = all.vol(:,:,57:end,i); % for a normal file, you would exclude this line and the following line
    wm(:,i) = tmp(:);
end
clear all

% now mask out areas outside the brain. NOTE, for your actual analysis
% please use a real brain mask that you have generated. After smoothing,
% voxels outside the brain may have incredibly low intensity values that
% could still be included using the heuristic below...
lm = wm~=0; % find nonzero voxels across participants
nz = sum(lm,2); % count number of participants with nonzero voxels
mid = find(nz > 150); % find voxels that are represented by at least 150 participants (>20% of data in this particular case)
wm = wm(mid,:);
gm = gm(mid,:);

% define variables for analysis. We will first try to run a partial
% correlation as the VBM analysis.
x = beh.NIHSS; % stroke severity
c = [beh.Site_num_ beh.Gender_num_ beh.Age beh.IscHemorrTia beh.wint_summ beh.PriorStroke beh.lesionSize beh.BMI beh.Chol_HDL beh.Chol_LDL beh.Choles beh.BPsyst beh.BPdias beh.Triglyc beh.Smoker beh.Diabetes beh.Hypertension beh.HeartAttack beh.FibOrFluttter beh.Dyslipidemia beh.CAD beh.Antihyperintensive beh.AntiGlucose beh.AntiPlatelet beh.IVtPA beh.Dysphagia beh.Depression beh.DrugAbuse beh.WeakenssParesis beh.AlteredConciousness beh.AphasiaLanguage];
% c above is covariates we want to control for

% You can ignore this--i did not have total intracranial volume available
% so made a variable that should correlate with it...
%mwm = sum(wm)';
%gwm = sum(gm)';
%c = [c mwm+gwm];

%% Running partial correlation (or correlation or distance correlation or partial distance correlation)
% run partial correlation--here we run them for two different sets of
% images--wm volume and gm volume. Note, parCorrVox always assumes the
% variable that comes after x is the one that contains brain data
[oR_nos,oP_nos,oRP_nos,oPP_nos] = parCorrVox(x,wm,c,2000,true,'partial');
[oR2_nos,oP2_nos,oRP2_nos,oPP2_nos] = parCorrVox(x,gm,c,2000,true,'partial');
% the scripts above are memory intensive and we save permutation results as
% we go. There's reason for this as we may at a later stage of the analysis
% want to apply nonparametric voxelwise thresholds without correcting for
% cluster size (e.g., when cluster size removes all effects).

% For the reason above, I would save the data now before anything crashes.
% We should really be running 10k permutations if possible by the way so
% bump that number up if you can
save('pc.mat','-v7.3')

% Before we move on, note that you can use parCorrVox to generate Pearson
% correlation coefficients like so: 
% [oR_nos,oP_nos,oRP_nos,oPP_nos] = parCorrVox(x,wm,[],2000,true,'pearson');
%
% You can also use parCorrVox to generate nonlinear correlation
% coefficients (capturing linear and nonlinear relations) like so: 
% [oR_nos,oP_nos,oRP_nos,oPP_nos] = parCorrVox(x,wm,[],2000,true,'distance');
%
% And finally you can use parCorrVox to generate partial nonlinear correlation
% coefficients (capturing linear and nonlinear relations) like so: 
% [oR_nos,oP_nos,oRP_nos,oPP_nos] = parCorrVox(x,wm,c(:,3),2000,true,'partialDistance');
%
% NOTE: you will get NaNs for distance correlation if you include variables
% in your analysis that have no variance (after automatic NaN removal
% inside parCorrVox) so choose, for example, covariates to control for
% carefully...

% The next thing we have to do is separate out positive and negative
% effects. Since we are cluster correcting, it's possible a large cluster
% displays adjacent positive and negative effects that get grouped
% together, which we would treat as distinct when we interpret the map...
[oR_nosNeg,oR_nosPos,oP_nosNeg,oP_nosPos] = sepPosNeg(oR_nos,oP_nos);
[oRP_nosNeg,oRP_nosPos,oPP_nosNeg,oPP_nosPos] = sepPosNeg(oRP_nos,oPP_nos);

[oR2_nosNeg,oR2_nosPos,oP2_nosNeg,oP2_nosPos] = sepPosNeg(oR2_nos,oP2_nos);
[oRP2_nosNeg,oRP2_nosPos,oPP2_nosNeg,oPP2_nosPos] = sepPosNeg(oRP2_nos,oPP2_nos);

% We can now take our true map and our permuted maps and threshold them,
% then extract the resulting maps' clusters and cluster sizes. To get cluster sizes we
% have to un-vectorize our brain maps which is why voxThresh requires you
% to tell it what the size of the original maps was before you vectorized
% them ([56 68 56] for us), and why it requires you tell it which of the
% voxels in the entire image were analyzed (mid)
[oR_nosPos,oP_nosPos,cclPos,~] = voxThresh(oR_nosPos,oP_nosPos,'default',0.05,[56 68 56],mid,true); % we get these values for true data. Note the outputs for oR and oP are now thresholded based on the 0.05 p-value threshold we selected!
[oR_nosNeg,oP_nosNeg,cclNeg,~] = voxThresh(oR_nosNeg,oP_nosNeg,'default',0.05,[56 68 56],mid,true); 

[oR_nosPPos,oP_nosPPos,~,cclMxPPos] = voxThresh(oRP_nosPos,oPP_nosPos,'default',0.05,[56 68 56],mid,false); % now we get these same cluster sizes for permuted data
[oR_nosPNeg,oP_nosPNeg,~,cclMxPNeg] = voxThresh(oRP_nosNeg,oPP_nosNeg,'default',0.05,[56 68 56],mid,false); % above we saved out all of the individual cluster voxel indices. However, we don't need these for the permutations so we set the last argument to false

% Now repeat the entire process above but for gm...
[oR2_nosPos,oP2_nosPos,ccl2Pos,cclMx2Pos] = voxThresh(oR2_nosPos,oP2_nosPos,'default',0.05,[56 68 56],mid,true); 
[oR2_nosNeg,oP2_nosNeg,ccl2Neg,cclMx2Neg] = voxThresh(oR2_nosNeg,oP2_nosNeg,'default',0.05,[56 68 56],mid,true); 

[oRP2_nosPos,oPP2_nosPos,~,cclMxP2Pos] = voxThresh(oRP2_nosPos,oPP2_nosPos,'default',0.05,[56 68 56],mid,false);
[oRP2_nosNeg,oPP2_nosNeg,~,cclMxP2Neg] = voxThresh(oRP2_nosNeg,oPP2_nosNeg,'default',0.05,[56 68 56],mid,false);

% Final analysis step--we need to apply the cluster-level threshold based on the
% permutation analysis...the appropriate percentile for permuted map
% cluster sizes is used to identify significant clusters in the true map
[oSvPos,oSmPos] = clustThresh(oR_nosPos,cclPos,cclMxPPos,0.05,mid); 
[oSvNeg,oSmNeg] = clustThresh(oR_nosNeg,cclNeg,cclMxPNeg,0.05,mid);
% let's also combine our maps...essentially sum them but check that there
% is no overlap between the maps (sanity check, there should not be...)
oSv = combPosNeg(oSvNeg,oSvPos);
oSm = combPosNeg(oSmNeg,oSmPos);

% repeat for gm since our data has two sets of maps we can test...
[oSv2Pos,oSm2Pos] = clustThresh(oR2_nosPos,ccl2Pos,cclMxP2Pos,0.05,mid);
[oSv2Neg,oSm2Neg] = clustThresh(oR2_nosNeg,ccl2Neg,cclMxP2Neg,0.05,mid);
oSv2 = combPosNeg(oSv2Pos,oSv2Neg);
oSm2 = combPosNeg(oSm2Neg,oSm2Pos);

% Now that the analysis is done we can save our outputs...
all = load_nifti('sall3mm_RH.nii.gz'); % load in a template file we can paste our data into...
all.vol = oSm;
save_nifti(all,'pcorr_test_wm.nii.gz');
all.vol = oSm2;
save_nifti(all,'pcorr_test_gm.nii.gz');

% Now go visualize your maps! If you downloaded brainSurfer for load_nifti
% and save_nifti, you could use the GUI to load these in and automatically
% project them onto the pial surface

% NOTE: voxThresh and medVox/parCorrVox can both be run in cases where we
% have roi instead of voxelwise data. However, you will have to find your
% own way of writing the outputs onto a brain. For help doing this, see
% ./misc directory of brainSurfer (specifically writeToAtlas.m)

%% Running mediation analysis
% This analysis proceeds the same way in general as the above example. Here
% too we will need to separate results by positive and negative values.
% This is because we quantify effects with t-values, and these can be
% negative. Note, medVox can keep many stats but we ignore them for
% permutation testing...if you want, you can save them out to look at
% later...

% first, let's define our new variables--now the brain matrices we
% initially set up will act as the mediator variables...
clearvars -except x beh wm mid % clear for memory
y = x; % NIH severity is our response variable
x = beh.BPdias; % blood pressure will be our predictor

% let's also define some confound variables. consider that the more nan
% values you have the lower your dof so choose condounders carefully
c = [beh.lesionSize beh.BMI];

% since gm did not show strong effects for the partial correlation, let's
% run mediation analysis on wm only. Note, brain data has to be the
% mediator (third input)
[oR,oP,~,~,oRP,oPP,~,~] = medVox(x,y,wm,c,1000,true); % for speed we will run very few permutations--consider doing 10k instead...note all the tildes representing voxelwise path coefficients, std, etc that we ignore.

% Voxelwise values are stored for each path in the model...here we only
% care about the indirect effect of the mediator but you should consider
% looking at the other estimates to get more insight (e.g., is it a complete
% mediator? etc...)
oR = oR(:,end);
oP = oP(:,end);
oRP = squeeze(oRP(:,end,:));
oPP = squeeze(oPP(:,end,:));

% here is where we need to separate out our positive and negative effects
% and threshold them
[oRNeg,oRPos,oPNeg,oPPos] = sepPosNeg(oR,oP);
[oRPNeg,oRPPos,oPPNeg,oPPPos] = sepPosNeg(oRP,oPP);

% Again,voxThresh is used to threshold your maps, but it also reshapes the
% vectorized brain data to identify clusters of contiguous voxels and
% tracks their size
[oRPos,oPPos,cclPos,cclMxPos] = voxThresh(oRPos,oPPos,'default',0.05,[56 68 56],mid,true); % we get these values for true data. Note the output oR_nos is now thresholded based on the 0.05 p-value threshold we selected! 
[oRNeg,oPNeg,cclNeg,cclMxNeg] = voxThresh(oRNeg,oPNeg,'default',0.05,[56 68 56],mid,true); 

[oRPPos,oPPPos,~,cclMxPPos] = voxThresh(oRPPos,oPPPos,'default',0.05,[56 68 56],mid,false); % now we get these values for permuted data...also set saving indivudal cluster voxel identities to false for the permuted data
[oRPNeg,oPPNeg,~,cclMxPNeg] = voxThresh(oRPNeg,oPPNeg,'default',0.05,[56 68 56],mid,false);

% And finally we can use the cluster sizes of the permuted data to identify
% significant clusters in our thresholded (unpermuted) map
[oSvPos,oSmPos] = clustThresh(oRPos,cclPos,cclMxPPos,0.05,mid); 
[oSvNeg,oSmNeg] = clustThresh(oRNeg,cclNeg,cclMxPNeg,0.05,mid);
% let's also combine our maps...essentially sum them but check that there
% is no overlap between the maps (sanity check)
oSv = combPosNeg(oSvNeg,oSvPos);
oSm = combPosNeg(oSmNeg,oSmPos);

% Now that the analysis is done we can save our outputs...
all = load_nifti('sall3mm_RH.nii.gz'); % load in a template file we can paste our data into...
all.vol = oSm;
save_nifti(all,'med_test_wm.nii.gz');
all.vol = oSm2;
save_nifti(all,'med_test_gm.nii.gz');