function res = axis(v, hierarchy, engine)
%AXIS Spearman r between a dispersion map and the cortical hierarchy.
%   res = sfcx.axis(v, hierarchy, engine)
%   Regional maps are compared with the regional hierarchy; edge maps with
%   the endpoint-mean hierarchy, and the null re-expands a rotated node map
%   so its spatial structure is inherited from the parcel layout.
v = double(v(:));
h = double(hierarchy(:));
if numel(h) ~= numel(v)
    error('sfcx:axis', 'hierarchy length %d does not match map length %d', numel(h), numel(v));
end
ok = isfinite(v) & isfinite(h);
r = spearmanR(v(ok), h(ok));
res = struct('r', r, 'p', NaN, 'null', 'none', 'nPerm', 0, 'n', sum(ok), ...
             'nullMean', NaN, 'nullSd', NaN, 'nNullValid', 0);

if nargin < 3 || isempty(engine)
    return
end
isEdge = numel(v) ~= engine.coh.nNodes;
nullR = nan(engine.nPerm, 1);
if ~isEdge
    vf = v;
    vf(~ok) = mean(v(ok));
    S = engine.surrogates(vf);
    for i = 1:engine.nPerm
        nullR(i) = spearmanR(S(i, ok), h(ok));
    end
else
    S = engine.surrogates(engine.coh.hierarchy);
    for i = 1:engine.nPerm
        eh = 0.5 * (S(i, engine.coh.iu) + S(i, engine.coh.iv));
        nullR(i) = spearmanR(v(ok), eh(ok));
    end
end
nullR = nullR(isfinite(nullR));
if ~isempty(nullR)
    res.p = twoSidedP(r, nullR);
    res.nullMean = mean(nullR);
    res.nullSd = std(nullR, 1);
    res.nNullValid = numel(nullR);
end
res.null = engine.kind;
res.nPerm = engine.nPerm;
end
