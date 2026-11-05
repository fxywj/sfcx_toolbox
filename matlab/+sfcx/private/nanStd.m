function s = nanStd(x, dim)
%NANSTD Population SD along dim, NaN-aware; all-NaN slices return NaN.
if nargin < 2
    dim = 1;
end
ok = isfinite(x);
n = sum(ok, dim);
xs = x;
xs(~ok) = 0;
mu = sum(xs, dim) ./ max(n, 1);
dev = x;
dev(~ok) = 0;
dev = (dev - mu) .^ 2;
ss = sum(dev, dim);
s = sqrt(ss ./ n);
s(n < 2) = NaN;
s = squeeze(s);
end
