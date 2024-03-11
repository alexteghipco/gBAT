function oNm = postpro(inCoeff,inP,pthresh,fdr,tail,montageUnderlay,gmThresh,anTail,atPth,inLim)
% requires load_nifti.m
% tail does NOT halve p-values at the moment, only removes positive or
% negative values
cin = load_nifti(inCoeff);
pin = load_nifti(inP);

id = find(cin.vol ~= 0);

p = zeros(size(pin.vol));
if fdr
    ptmp = mafdr(pin.vol(id),'BHFDR',true);
else
    ptmp = pin.vol(id);
end
p(id) = ptmp;

if fdr
    pin.vol = p;
    [p1,p2,p3] = fileparts(inP);
    if strcmpi(p3,'.gz')
        p2 = p2(1:end-4);
    end
    save_nifti(pin,[p1 filesep p2 '_corrected_FDR.nii.gz']);
end

id = find(p > pthresh);
cin.vol(id) = 0;
[p1,p2,p3] = fileparts(inCoeff);
if strcmpi(p3,'.gz')
    p2 = p2(1:end-4);
end

if strcmpi(tail,'neg')
    id = find(cin.vol > 0);
    cin.vol(id) = 0;
elseif strcmpi(tail,'pos')
    id = find(cin.vol < 0);
    cin.vol(id) = 0;
end
save_nifti(cin,[p1 filesep p2 '_corrected_FDR_' num2str(pthresh) '_tail_' tail '.nii.gz']);
oNm = [p1 filesep p2 '_corrected_FDR_' num2str(pthresh) '_tail_' tail '.nii.gz'];

if ~isempty(montageUnderlay)
    [tmp1,tmp2,~,~,olim] = brainMontagerWithAtlas([p1 filesep p2 '_corrected_FDR_' num2str(pthresh) '_tail_' tail '.nii.gz'],montageUnderlay,atPth,round(linspace(1,size(cin.vol,3),20)),1000,inLim,'viridis',[],[1 1 1],false,1000,[],'gray',[],0.8,'ud',true,false,true,gmThresh,anTail);
    %title([p1 filesep p2 '_corrected_FDR_' num2str(pthresh) '_tail_' tail '.nii.gz'])
    saveas(gcf,[p1 filesep p2 '_corrected_FDR_' num2str(pthresh) '_tail_' tail '_montage_lims_' num2str(min(olim)) '_' num2str(max(olim))  '.fig']);
end


