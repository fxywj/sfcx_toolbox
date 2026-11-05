function [fcVec, nNodes] = boldToFc(bold, kind, timeAxis)
%BOLDTOFC Upper-triangle FC edge vector from BOLD series.
%   [fcVec, nNodes] = sfcx.boldToFc(bold, 'pearson')
if nargin < 2 || isempty(kind)
    kind = 'pearson';
end
if nargin < 3 || isempty(timeAxis)
    timeAxis = 'auto';
end
mats = sfcx.boldToMatrices(bold, kind, timeAxis);
n = size(mats, 1);
nNodes = size(mats, 2);
[iu, iv] = upperTriIndex(nNodes);
X = reshape(mats, n, nNodes * nNodes);
fcVec = X(:, sub2ind([nNodes nNodes], iu, iv));
end
