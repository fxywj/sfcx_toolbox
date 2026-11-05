function m = nanMeanVec(x, dim)
%NANMEANVEC NaN-aware mean that keeps NaN for all-NaN slices.
%   Row/column vectors are reduced to a scalar; matrices are reduced along
%   DIM (default 1).
if nargin < 2 || isempty(dim)
    dim = 1;
end
if isempty(x)
    m = NaN;
    return
end
if isvector(x)
    dim = find(size(x) > 1, 1);
    if isempty(dim)
        dim = 1;
    end
end
ok = isfinite(x);
n = sum(ok, dim);
xs = x;
xs(~ok) = 0;
m = sum(xs, dim) ./ n;
m(n == 0) = NaN;
end
