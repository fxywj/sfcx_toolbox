function res = partialSpearman(x, y, cov)
%PARTIALSPEARMAN Rank correlation, optionally partialling out covariates.
x = double(x(:));
y = double(y(:));
if numel(x) ~= numel(y)
    error('sfcx:partialSpearman', 'x and y must have the same length');
end
ok = isfinite(x) & isfinite(y);
if nargin >= 3 && ~isempty(cov)
    cov = double(cov);
    if isvector(cov)
        cov = cov(:);
    end
    ok = ok & all(isfinite(cov), 2);
end
rx = rankVector(x(ok));
ry = rankVector(y(ok));
if nargin >= 3 && ~isempty(cov)
    c = zeros(sum(ok), size(cov, 2));
    for k = 1:size(cov, 2)
        c(:, k) = rankVector(cov(ok, k));
    end
    rx = residualise(rx, c);
    ry = residualise(ry, c);
end
n = numel(rx);
sa = std(rx, 1);
sb = std(ry, 1);
if n < 3 || sa < 1e-12 || sb < 1e-12
    r = NaN;
    p = NaN;
else
    r = mean((rx - mean(rx)) .* (ry - mean(ry))) / (sa * sb);
    if abs(r) >= 1
        p = 0;
    else
        p = betainc(1 - r .^ 2, (n - 2) / 2, 0.5);
    end
end
res = struct('r', r, 'p', p, 'n', n);
end
