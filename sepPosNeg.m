function [sNeg,sPos,pNeg,pPos] = sepPosNeg(s,p)
% This function findgs negative and positive values in s to generate
% separate positive and negative matrices/vectors for s and p

id = find(s > 0); sPos = zeros(size(s)); sPos(id) = s(id); pPos = zeros(size(p)); pPos(id) = p(id);
id = find(s < 0); sNeg = zeros(size(s)); sNeg(id) = s(id); pNeg = zeros(size(p)); pNeg(id) = p(id);
