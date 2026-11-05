function out = relianceDerivNonlinear(coh, alpha)
%RELIANCEDERIVNONLINEAR Nonlinear derivative reliance: corr(a*(tanh(ts*W)-ts), dts/dt).
if isempty(coh.bold)
    error('sfcx:relianceDerivNonlinear', 'requires BOLD time series in the cohort');
end
if nargin < 2 || isempty(alpha)
    alpha = 1.0;
end
n = coh.n;
p = coh.nNodes;
W = colNormalize(coh.groupMatrix('sc'));
out = zeros(n, p);
for s = 1:n
    ts = squeeze(coh.bold(s, :, :))';  % (T, nNodes)
    dy = ts(2:end, :) - ts(1:end-1, :);
    yCurr = ts(1:end-1, :);
    drive = yCurr * W;
    dyPred = alpha * (tanh(drive) - yCurr);
    for i = 1:p
        out(s, i) = robustPearsonr(dyPred(:, i), dy(:, i));
    end
end
end
