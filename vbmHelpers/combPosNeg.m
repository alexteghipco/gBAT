function s = combPosNeg(sNeg,sPos)
% This function combines sNeg and sPos vectors into one by assuming that
% there is no overlap among them...
id1 = find(sNeg ~= 0);
id2 = find(sPos ~= 0);

if ~isempty(intersect(id1,id2))
    error('Your maps have overlap across rows...')
end
s = sNeg + sPos;