function [kind, field] = evaluateCustom(fn, coh)
%EVALUATECUSTOM Run a user function and coerce its output to a coupling field.
%   Accepted signatures:  f(sc, fc, iu, iv, nNodes)  or  f(coh).
n = coh.n;
e = coh.nEdges;
p = coh.nNodes;
if nargin(fn) == 5
    out = fn(double(coh.sc), double(coh.fc), coh.iu, coh.iv, coh.nNodes);
else
    out = fn(coh);
end
out = double(out);
if isscalar(out)
    kind = 'edge';
    field = repmat(out, n, e);
elseif isequal(size(out), [n, e])
    kind = 'edge';
    field = out;
elseif numel(out) == e && (size(out, 1) == 1 || size(out, 2) == 1)
    kind = 'edge';
    field = repmat(out(:)', n, 1);
elseif isequal(size(out), [n, p])
    kind = 'regional';
    field = out;
elseif numel(out) == p && (size(out, 1) == 1 || size(out, 2) == 1)
    kind = 'regional';
    field = repmat(out(:)', n, 1);
elseif size(out, 1) == n && size(out, 2) == p && size(out, 3) == p
    kind = 'edge';
    field = reshape(out(:, sub2ind([p p], coh.iu, coh.iv)), n, e);
elseif isequal(size(out), [p, p])
    kind = 'edge';
    field = repmat(out(sub2ind([p p], coh.iu, coh.iv))', n, 1);
else
    error('sfcx:custom', ...
        ['return size [%s] is not usable; expected [%d %d] (subjects x edges), ' ...
         '[%d 1] (edges), [%d 1] (nodes) or [%d %d]'], ...
        num2str(size(out)), n, e, e, p, p, p);
end
if all(isnan(field(:)))
    error('sfcx:custom', 'custom function returned all-NaN values');
end
end
