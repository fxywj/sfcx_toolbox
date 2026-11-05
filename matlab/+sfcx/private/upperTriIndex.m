function [iu, iv] = upperTriIndex(nNodes)
%UPPERTRIINDEX Edge endpoints in Python/C order (row-major upper triangle).
%   1-based indices, so they can index MATLAB arrays directly.
nEdges = nNodes * (nNodes - 1) / 2;
iu = zeros(nEdges, 1);
iv = zeros(nEdges, 1);
k = 0;
for i = 1:(nNodes - 1)
    m = nNodes - i;
    iu(k + (1:m)) = i;
    iv(k + (1:m)) = (i + 1):nNodes;
    k = k + m;
end
end
