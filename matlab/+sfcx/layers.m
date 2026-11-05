function res = layers(edgeV, coh, nPerm, seed)
%LAYERS Within- versus between-network dispersion of an edge-level map.
if nargin < 3 || isempty(nPerm)
    nPerm = 1000;
end
if nargin < 4 || isempty(seed)
    seed = 0;
end
v = double(edgeV(:));
net = coh.networks(:);
same = net(coh.iu) == net(coh.iv);
w = nanMeanVec(v(same), 1);
b = nanMeanVec(v(~same), 1);
delta = w - b;

stream = RandStream('mt19937ar', 'Seed', seed);
null = nan(nPerm, 1);
for i = 1:nPerm
    lab = net(randperm(stream, numel(net)));
    s = lab(coh.iu) == lab(coh.iv);
    null(i) = nanMeanVec(v(s), 1) - nanMeanVec(v(~s), 1);
end
res = struct('within', w, 'between', b, 'delta', delta, ...
             'p', twoSidedP(delta, null), 'nEdges', numel(v), ...
             'sameNetEdges', sum(same), 'null', 'network-label permutation', ...
             'nPerm', nPerm);
end
