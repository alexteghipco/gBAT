function atlasVol2Surf(inVol,inVolAtlas,inSurfAtlasLH,inSurfAtlasRH)

volin = load_nifti(inVol);
atin = load_nifti(inVolAtlas);
lhin = load_nifti(inSurfAtlasLH);
rhin = load_nifti(inSurfAtlasRH);

un = unique(atin.vol);
id = find(un == 0);
un(id) = [];
lho = lhin;
rho = rhin;
lho.vol = zeros(size(lho.vol));
rho.vol = zeros(size(rho.vol));

for i = 1:length(un)
    id1 = find(atin.vol == un(i));
    id2 = find(lhin.vol == un(i));
    id3 = find(rhin.vol == un(i));

    if isempty(id3) && ~isempty(id2)
        lho.vol(id2) = unique(volin.vol(id1));
    elseif isempty(id2) && ~isempty(id3)
        rho.vol(id3) = unique(volin.vol(id1));
    end
end
[p1,p2,p3] = fileparts(inSurfAtlasLH);
if strcmpi(p3,'.gz')
    p2 = p2(1:end-4);
end
[p1v,p2v,p3v] = fileparts(inVol);

if strcmpi(p3v,'.gz')
    p2v = p2v(1:end-4);
end

save_nifti(lho,[p1v filesep p2v p2])
[p1,p2,p3] = fileparts(inSurfAtlasRH);
if strcmpi(p3,'.gz')
    p2 = p2(1:end-4);
end
save_nifti(rho,[p1v filesep p2v p2])