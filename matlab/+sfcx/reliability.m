function res = reliability(edgeField, coh, seed)
%RELIABILITY Split-half reliability of the edge-level dispersion map.
if nargin < 3 || isempty(seed)
    seed = 0;
end
[a, b] = coh.splitHalf(seed);
va = nanStd(double(edgeField(a, :)), 1);
vb = nanStd(double(edgeField(b, :)), 1);
ok = isfinite(va) & isfinite(vb);
res = struct('splitHalfR', spearmanR(va(ok), vb(ok)), ...
             'nEdges', sum(ok), 'nHalf', numel(a));
end
