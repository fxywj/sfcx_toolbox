function mats = boldToMatrices(bold, kind, timeAxis)
%BOLDMATRICES Per-subject FC matrices from BOLD time series.
%   mats = sfcx.boldToMatrices(bold, 'pearson')
%   bold : (subjects, nodes, time) | (nodes, time) | (time, nodes)
%   kind : 'pearson' (full correlation) | 'precision' (Ledoit-Wolf partial)
if nargin < 2 || isempty(kind)
    kind = 'pearson';
end
if nargin < 3 || isempty(timeAxis)
    timeAxis = 'auto';
end
x = asBold3d(bold, timeAxis);
x = interpolateMissing(x);
[n, p, t] = size(x);
mats = zeros(n, p, p);
switch lower(kind)
    case 'pearson'
        z = standardise(x);
        for s = 1:n
            zs = reshape(z(s, :, :), p, t);
            mats(s, :, :) = reshape((zs * zs') / t, 1, p, p);
        end
    case {'precision', 'partial'}
        for s = 1:n
            xs = reshape(x(s, :, :), p, t)';
            mats(s, :, :) = reshape(precisionFromCov(ledoitWolfCov(xs)), 1, p, p);
        end
    otherwise
        error('sfcx:bold', 'unknown functional pipeline ''%s''', kind);
end
end

function x = asBold3d(bold, timeAxis)
x = double(bold);
if ndims(x) == 2
    if strcmpi(timeAxis, 'time')
        x = x';
    elseif strcmpi(timeAxis, 'nodes')
        % already nodes x time
    elseif size(x, 2) < size(x, 1)
        x = x';
    end
    x = reshape(x, 1, size(x, 1), size(x, 2));
elseif ndims(x) == 3
    if strcmpi(timeAxis, 'time')
        x = permute(x, [1 3 2]);
    end
else
    error('sfcx:bold', 'BOLD must be 2-D or 3-D');
end
if size(x, 3) < 3
    error('sfcx:bold', 'BOLD needs at least 3 time points, got %d', size(x, 3));
end
end

function x = interpolateMissing(x)
if all(isfinite(x(:)))
    return
end
[n, p, t] = size(x);
tt = (1:t)';
for s = 1:n
    for i = 1:p
        row = reshape(x(s, i, :), 1, t);
        bad = ~isfinite(row);
        if ~any(bad)
            continue;
        end
        good = ~bad;
        if sum(good) < 2
            row(bad) = 0;
        else
            row(bad) = interp1(tt(good), row(good), tt(bad), 'linear');
        end
        x(s, i, :) = reshape(row, 1, 1, t);
    end
end
end

function z = standardise(x)
x = x - mean(x, 3);
sd = sqrt(sum(x .^ 2, 3));
flat = sd < 1e-12;
sd(flat) = 1;
z = x ./ sd;
z(repmat(flat, [1, 1, size(x, 3)])) = 0;
end

function s = ledoitWolfCov(x)
[t, p] = size(x);
xc = x - mean(x, 1);
s = (xc' * xc) / t;
mu = trace(s) / p;
target = mu * eye(p);
d2 = sum((s(:) - target(:)) .^ 2);
if d2 < 1e-18
    s = target;
    return
end
b2 = 0;
for k = 1:t
    outer = xc(k, :)' * xc(k, :);
    b2 = b2 + sum((outer(:) - s(:)) .^ 2);
end
b2 = b2 / (t * t);
lam = min(b2 / d2, 1);
s = lam * target + (1 - lam) * s;
end

function pc = precisionFromCov(cov)
p = size(cov, 1);
d = sqrt(max(diag(cov), 1e-12));
cor = cov ./ (d * d');
cor(1:p + 1:end) = 1;
if rcond(cor) < 1e-15
    invCor = pinv(cor);
else
    invCor = inv(cor);
end
dd = sqrt(max(abs(diag(invCor)), 1e-12));
pc = -invCor ./ (dd * dd');
pc(1:p + 1:end) = 0;
end
