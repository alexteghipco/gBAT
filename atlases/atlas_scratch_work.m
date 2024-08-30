%% THIS CODE SHOWS HOW WE CREATED THE ATLASES REDISTRIBUTED WITH GBAT
%% HO atlas
% load in cortical and subcortical versions
hocort = load_nifti('D:\Science\Matlab\GitHub\VBM\atlases\HO\hocorticalOnly.nii.gz');
hosubcort = load_nifti('D:\Science\Matlab\GitHub\VBM\atlases\HO\hosubcorticalOnly.nii.gz');
hocortROI = readtable('D:\Science\Matlab\GitHub\VBM\atlases\HO\hocorticalOnly.xml');
hosubcortROI = readtable('D:\Science\Matlab\GitHub\VBM\atlases\HO\hosubcorticalOnly.xml');

% convert character arrays for convenience and adjust indexing to match
% matlab's
hocortROI.label = cellstr(hocortROI.label);
hocortROI.indexAttribute = hocortROI.indexAttribute+1;
hosubcortROI.label = cellstr(hosubcortROI.label);
hosubcortROI.indexAttribute = hosubcortROI.indexAttribute+1;

% create a table that splits left and right cortical labels
hocortROI.label = arrayfun(@(x) ['Left ' x{1}], hocortROI.label, 'UniformOutput', false);
hocortROI2 = hocortROI;
hocortROI2.label = extractAfter(hocortROI2.label,'Left ');
hocortROI2.label = arrayfun(@(x) ['Right ' x{1}], hocortROI2.label, 'UniformOutput', false);
hocortROI2.indexAttribute = hocortROI2.indexAttribute+max(hocortROI.indexAttribute);

% the subcortical atlas gives us the left/right splits--get voxels
id = find(contains(hosubcortROI.label,'Right'));
rid = find(ismember(hosubcort.vol, hosubcortROI.indexAttribute(id)));
id = find(contains(hosubcortROI.label,'Left'));
lid = find(ismember(hosubcort.vol, hosubcortROI.indexAttribute(id)));

% now split left and right hemisphere regions from the cortical atlas
for i = 1:height(hocortROI)
    id = find(hocort.vol == hocortROI.indexAttribute(i)); % find all voxels for region in bilateral atlas
    [c,ia,ib] = intersect(id,rid);
    hocort.vol(c) = hocortROI2.indexAttribute(i);
end
save_nifti(hocort,'D:\Science\Matlab\GitHub\VBM\atlases\HO\hocort_split.nii.gz')

% now combined tables for left and right cortical 
hocortROI = vertcat(hocortROI,hocortROI2);

% create output table for combined atlas
hoROI = table();
hoROI.Index = hocortROI.indexAttribute;
hoROI.Label = hocortROI.label;

% adjust the values of the subcortical atlas based on the number of ROIs in
% the cortical atlas
mx = max(hocort.vol(:));
hosubcortROI.indexAttribute = hosubcortROI.indexAttribute+mx;
id = find(hosubcort.vol ~= 0);
hosubcort.vol(id) = hosubcort.vol(id)+mx;
% mn = min(hosubcort.vol(id)); % check that this is 1 (indexing isn't unusual 

% combine the tables for both atlases
h = height(hoROI)+1;
hoROI.Index(h:h-1+length(hosubcortROI.indexAttribute)) = hosubcortROI.indexAttribute;
hoROI.Label(h:h-1+length(hosubcortROI.label)) = hosubcortROI.label;

% define areas to fill
id = find(hocort.vol == 0 & hosubcort.vol ~= 0);
hocort.vol(id) = hosubcort.vol(id);
save_nifti(hocort,'D:\Science\Matlab\GitHub\VBM\atlases\HO\ho.nii.gz')

% classify tissue types
tis = load_nifti('D:\Science\Matlab\GitHub\VBM\tissue\MNI152_T1_2mm_brain_pveseg.nii.gz');
un = setdiff(unique(hocort.vol),0); % find all unique integers that are not zero
for i = 1:length(un)
    id = find(hocort.vol == un(i));
    idcsf = length(find(tis.vol(id) == 1));
    if isempty(idcsf)
        idcsf = 0;
    end
    %idcsf= 0;
    idgm = length(find(tis.vol(id) == 2));
    if isempty(idgm)
        idgm = 0;
    end
    idwm = length(find(tis.vol(id) == 3));
    if isempty(idwm)
        idwm = 0;
    end
    [~,hoROI.Tissue(i)] = max([idcsf idgm idwm]);
end

% write out txt files with regions based on our formatting
for i = 1:size(hoROI,1)
    txt{i,1} = [num2str(hoROI.Index(i)) '|' '|' hoROI.Label{i} '|' num2str(hoROI.Tissue(i))];
end

fileID = fopen(['D:\Science\Matlab\GitHub\VBM\atlases\HO\ho.txt'], 'w');
for i = 1:length(txt)
    fprintf(fileID, '%s\n', txt{i});
end
fclose(fileID);

txt1 = txt(1:height(hocortROI));
fileID = fopen(['D:\Science\Matlab\GitHub\VBM\atlases\HO\hocort_split.txt'], 'w');
for i = 1:length(txt1)
    fprintf(fileID, '%s\n', txt1{i});
end
fclose(fileID);

txt1 = txt(height(hocortROI)+1:end);
fileID = fopen(['D:\Science\Matlab\GitHub\VBM\atlases\HO\hosubcorticalOnly.txt'], 'w');
for i = 1:length(txt1)
    fprintf(fileID, '%s\n', txt1{i});
end
fclose(fileID);

tmp = load_nifti('D:\Science\Matlab\GitHub\VBM\atlases\HO\HarvardOxford-cort-maxprob-thr0-2mm_CombinedClusters_FSSpace_RH.nii.gz');
id = find(tmp.vol ~= 0);
tmp.vol(id) = tmp.vol(id)+height(hocortROI2);
save_nifti(tmp,['D:\Science\Matlab\GitHub\VBM\atlases\HO\HarvardOxford-cort-maxprob-thr0-2mm_CombinedClusters_FSSpace_RH.nii.gz'])

% We need to adjust our table so that 1 is GM, 2 is WM, everything else is
% 3
fileID = fopen('D:\Science\Matlab\GitHub\VBM\atlases\HO\ho.txt', 'r');
txt1Read = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
txt1Read = txt1Read{1};
for i = 1:size(txt1Read,1)
    id = strfind(txt1Read{i},'|');
    if str2num(txt1Read{i}(id(end)+1:end)) == 1
        txt1Read{i} = [txt1Read{i}(1:id(end)) '3'];
    elseif str2num(txt1Read{i}(id(end)+1:end)) == 2
        txt1Read{i} = [txt1Read{i}(1:id(end)) '1'];
    elseif str2num(txt1Read{i}(id(end)+1:end)) == 3
        txt1Read{i} = [txt1Read{i}(1:id(end)) '2'];
    end
end
fileID = fopen(['D:\Science\Matlab\GitHub\VBM\atlases\HO\ho2.txt'], 'w');
for i = 1:length(txt1Read)
    fprintf(fileID, '%s\n', txt1Read{i});
end
fclose(fileID);

fileID = fopen('D:\Science\Matlab\GitHub\VBM\atlases\HO\hosubcorticalOnly.txt', 'r');
txt1Read = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
txt1Read = txt1Read{1};
for i = 1:size(txt1Read,1)
    id = strfind(txt1Read{i},'|');
    if str2num(txt1Read{i}(id(end)+1:end)) == 1
        txt1Read{i} = [txt1Read{i}(1:id(end)) '3'];
    elseif str2num(txt1Read{i}(id(end)+1:end)) == 2
        txt1Read{i} = [txt1Read{i}(1:id(end)) '1'];
    elseif str2num(txt1Read{i}(id(end)+1:end)) == 3
        txt1Read{i} = [txt1Read{i}(1:id(end)) '2'];
    end
end
fileID = fopen(['D:\Science\Matlab\GitHub\VBM\atlases\HO\hosubcorticalOnly2.txt'], 'w');
for i = 1:length(txt1Read)
    fprintf(fileID, '%s\n', txt1Read{i});
end
fclose(fileID);

fileID = fopen('D:\Science\Matlab\GitHub\VBM\atlases\HO\hocort_split.txt', 'r');
txt1Read = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
txt1Read = txt1Read{1};
for i = 1:size(txt1Read,1)
    id = strfind(txt1Read{i},'|');
    if str2num(txt1Read{i}(id(end)+1:end)) == 1
        txt1Read{i} = [txt1Read{i}(1:id(end)) '3'];
    elseif str2num(txt1Read{i}(id(end)+1:end)) == 2
        txt1Read{i} = [txt1Read{i}(1:id(end)) '1'];
    elseif str2num(txt1Read{i}(id(end)+1:end)) == 3
        txt1Read{i} = [txt1Read{i}(1:id(end)) '2'];
    end
end
fileID = fopen(['D:\Science\Matlab\GitHub\VBM\atlases\HO\hocort_split2.txt'], 'w');
for i = 1:length(txt1Read)
    fprintf(fileID, '%s\n', txt1Read{i});
end
fclose(fileID);

%% XTRACT
roi = readtable('D:\Science\Matlab\GitHub\VBM\atlases\XTRACT\XTRACT.xml');
roi.Tissue = ones(height(roi),1)*2;
clear txt
for i = 1:size(roi,1)
    txt{i,1} = [num2str(roi.indexAttribute(i)+1) '|' '|' char(roi.label(i)) '|' num2str(roi.Tissue(i))];
end
fileID = fopen(['D:\Science\Matlab\GitHub\VBM\atlases\XTRACT\XTRACT.txt'], 'w');
for i = 1:length(txt)
    fprintf(fileID, '%s\n', txt{i});
end
fclose(fileID);

%% Julich
latlas = load_nifti('D:\Science\Matlab\GitHub\VBM\atlases\Julich\JulichBrainAtlas_3.1_207areas_MPM_lh_MNI152_resampled_2mm.nii.gz');
ratlas = load_nifti('D:\Science\Matlab\GitHub\VBM\atlases\Julich\JulichBrainAtlas_3.1_207areas_MPM_rh_MNI152_resampled2mm.nii.gz');
lroi = readtable('D:\Science\Matlab\GitHub\VBM\atlases\Julich\JulichBrainAtlas_3.1_207areas_MPM_lh_MNI152.xml');
rroi = readtable('D:\Science\Matlab\GitHub\VBM\atlases\Julich\JulichBrainAtlas_3.1_207areas_MPM_rh_MNI152.xml');
rroi.grayvalueAttribute = rroi.grayvalueAttribute+max(latlas.vol(:));
lroi.Structure = arrayfun(@(x) ['Left ' x{1}], lroi.Structure, 'UniformOutput', false);
rroi.Structure = arrayfun(@(x) ['Right ' x{1}], rroi.Structure, 'UniformOutput', false);
roi = vertcat(lroi,rroi);
id = find(ratlas.vol ~= 0);
ratlas.vol(id) = ratlas.vol(id)+max(latlas.vol(:));
id = find(ratlas.vol ~= 0 & latlas.vol ~= 0); % check that there is no overlap...there is not
ratlas.vol = ratlas.vol+latlas.vol;
save_nifti(ratlas,'D:\Science\Matlab\GitHub\VBM\atlases\Julich\JulichBrainAtlas_3.1_207areas_MPM_MNI152_resampled_2mm.nii.gz')

roi.Tissue = ones(height(roi),1);
clear txt
for i = 1:size(roi,1)
    txt{i,1} = [num2str(roi.grayvalueAttribute(i)) '|' '|' char(roi.Structure(i)) '|' num2str(roi.Tissue(i))];
end
fileID = fopen(['C:\Users\alext\Downloads\Julich\Julich\Julich.txt'], 'w');
for i = 1:length(txt)
    fprintf(fileID, '%s\n', txt{i});
end
fclose(fileID);
g = gifti('C:\Users\alext\Downloads\Julich\Julich\rh.JulichBrainAtlas_3.1.label.gii');
g.cdata = int32(g.cdata+max(latlas.vol(:)));
g2 = gifti('D:\Science\Matlab\GitHub\brainSurfer\atlases\jhu_resampled_to_FSL_fixed_noWM_CombinedClusters_FSSpace_LH.nii_to_fs_LR.label.gii');
g2.cdata = g.cdata;
save(g2,'C:\Users\alext\Downloads\Julich\Julich\rh.JulichBrainAtlas_3.1.label_fixed2.gii','ExternalFileBinary')

addpath('D:\Science\Matlab\GitHub\brainSurfer\scripts\HCP\cifti\ft_cifti\@gifti');

%% AICHA AND JHU ATLASES ARE EDITED VERSIONS OF THOSE SUPPLIED BY NIISTAT. THE MULTIRES ATLAS IS FROM LQT (see manual)
