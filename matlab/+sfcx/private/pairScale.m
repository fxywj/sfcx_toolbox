function out = pairScale(m)
%PAIRSCALE Log-compress and rescale a non-negative pair matrix.
m = double(m);
m(~isfinite(m)) = 0;
m = log1p(max(m, 0));
sd = std(m(:), 1);
if sd > 1e-12
    out = m / sd;
else
    out = m;
end
end
