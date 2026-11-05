function [M, counts] = incidence(nNodes)
%INCIDENCE Sparse region-edge incidence matrix, cached by node count.
persistent lastN Mlast CLast
if ~isempty(lastN) && lastN == nNodes
    M = Mlast;
    counts = CLast;
    return
end
nEdges = nNodes * (nNodes - 1) / 2;
[iu, iv] = upperTriIndex(nNodes);
cols = (1:nEdges)';
M = sparse([iu; iv], [cols; cols], 1, nNodes, nEdges);
counts = full(sum(M, 2));
counts(counts < 1) = 1;
lastN = nNodes;
Mlast = M;
CLast = counts;
end
