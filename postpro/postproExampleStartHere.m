% get volume-based results on the brain
allp = 0.05;
bd = 'D:\Science\Matlab\GitHub\VBM';
fnms = {'TotalScoreFin','OrientationFin','MemoryFin','VisuoFin','ExecFin','LangFin','Attention_MoCA'};
jhuFnm = 'D:\Science\Matlab\GitHub\VBM\rjhu.nii';
xname = 'Age';
cin = load_nifti(jhuFnm);

% apply new FDR thresholds to outputs for mediation effect and generate
% brain figures
allp = [0.01 0.005 0.001 0.0001];
bd = 'D:\Science\Matlab\GitHub\VBM\';
fnms = {'TotalScoreFin','OrientationFin','MemoryFin','VisuospatgialFin','ExecFin','LangFin','Attention_MoCA'};
for i = 1:length(allp)
    for j = 1:length(fnms)
        nm{i,j} = postpro([bd filesep fnms{j} filesep 'Age_mediation_ab_coeff.nii.gz'],[bd filesep fnms{j} filesep 'Age_mediation_ab_pval.nii.gz'],allp(i),true,'neg','D:\Science\Matlab\GitHub\VBM\MemoryFin\meanSubMap_GM.nii.gz',0.1,'neg','D:\Science\Matlab\GitHub\VBM\rjhu.nii',[]); % -0.0000001 -0.015
        close all
        atlasVol2Surf(nm{i,j},'D:\Science\Matlab\GitHub\VBM\rjhu.nii','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_LH.nii.gz','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_RH.nii.gz');
        atlasVol2Surf('D:\Science\Matlab\GitHub\VBM\TotalScoreFin\Age_mediation_ab_coeff_corrected_FDR_0.005_tail_neg.nii.gz','D:\Science\Matlab\GitHub\VBM\rjhu.nii','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_LH.nii.gz','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_RH.nii.gz');
    end
end

% same with other coefficients...
allp = [0.01 0.005 0.001 0.0001];
bd = 'D:\Science\Matlab\GitHub\VBM\';
fnms = {'TotalScoreFin','OrientationFin','MemoryFin','VisuospatgialFin','ExecFin','LangFin','Attention_MoCA'};
for i = 1:length(allp)
    for j = 1:length(fnms)
        nm{i,j} = postpro([bd filesep fnms{j} filesep 'Age_mediation_a_coeff.nii.gz'],[bd filesep fnms{j} filesep 'Age_mediation_a_pval.nii.gz'],allp(i),true,'neg','D:\Science\Matlab\GitHub\VBM\MemoryFin\meanSubMap_GM.nii.gz',0.1,'neg','D:\Science\Matlab\GitHub\VBM\rjhu.nii',[-0.000001 -0.0005]); % -0.0000001 -0.015
        close all
        atlasVol2Surf(nm{i,j},'D:\Science\Matlab\GitHub\VBM\rjhu.nii','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_LH.nii.gz','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_RH.nii.gz');
    end
end
for i = 1:length(allp)
    for j = 1:length(fnms)
        nm{i,j} = postpro([bd filesep fnms{j} filesep 'Age_mediation_c''_coeff.nii.gz'],[bd filesep fnms{j} filesep 'Age_mediation_c''_pval.nii.gz'],allp(i),true,'neg','D:\Science\Matlab\GitHub\VBM\MemoryFin\meanSubMap_GM.nii.gz',0.1,'neg','D:\Science\Matlab\GitHub\VBM\rjhu.nii',[]); % -0.0000001 -0.015
        close all
        atlasVol2Surf(nm{i,j},'D:\Science\Matlab\GitHub\VBM\rjhu.nii','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_LH.nii.gz','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_RH.nii.gz');
    end
end
for i = 1:length(allp)
    for j = 1:length(fnms)
        nm{i,j} = postpro([bd filesep fnms{j} filesep 'Age_mediation_b_coeff.nii.gz'],[bd filesep fnms{j} filesep 'Age_mediation_b_pval.nii.gz'],allp(i),true,'pos','D:\Science\Matlab\GitHub\VBM\MemoryFin\meanSubMap_GM.nii.gz',0.1,'pos','D:\Science\Matlab\GitHub\VBM\rjhu.nii',[]); % -0.0000001 -0.015
        close all
        atlasVol2Surf(nm{i,j},'D:\Science\Matlab\GitHub\VBM\rjhu.nii','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_LH.nii.gz','D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_RH.nii.gz');
    end
end
