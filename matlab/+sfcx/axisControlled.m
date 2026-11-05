function res = axisControlled(v, hierarchy, control, engine)
%AXISCONTROLLED Partial Spearman of v and hierarchy given a control map.
v = double(v(:));
h = double(hierarchy(:));
c = double(control(:));
ok = isfinite(v) & isfinite(h) & isfinite(c);
rv = rankVector(v(ok));
rh = rankVector(h(ok));
rc = rankVector(c(ok));
ev = residualise(rv, rc);
eh = residualise(rh, rc);
r = corrCoeff(ev, eh);
res = struct('r', r, 'p', NaN, 'null', 'none', 'nPerm', 0, 'n', sum(ok));

if nargin < 4 || isempty(engine)
    return
end
vf = v;
vf(~ok) = mean(v(ok));
S = engine.surrogates(vf);
nullR = nan(engine.nPerm, 1);
for i = 1:engine.nPerm
    evn = residualise(rankVector(S(i, ok)'), rc);
    nullR(i) = corrCoeff(evn, eh);
end
nullR = nullR(isfinite(nullR));
if ~isempty(nullR)
    res.p = twoSidedP(r, nullR);
end
res.null = engine.kind;
res.nPerm = engine.nPerm;
end

function r = corrCoeff(a, b)
a = a(:);
b = b(:);
sa = std(a, 1);
sb = std(b, 1);
if sa < 1e-12 || sb < 1e-12
    r = NaN;
    return
end
r = mean((a - mean(a)) .* (b - mean(b))) / (sa * sb);
end
