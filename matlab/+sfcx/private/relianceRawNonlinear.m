function out = relianceRawNonlinear(coh, alpha)
%RELIANCERAWNONLINEAR Nonlinear reliance (leak + tanh): corr((1-a)ts + a*tanh(ts*W), ts_t).
if isempty(coh.bold)
    error('sfcx:relianceRawNonlinear', 'requires BOLD time series in the cohort');
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
    drive = ts(1:end-1, :) * W;
    yPred = (1 - alpha) * ts(1:end-1, :) + alpha * tanh(drive);
    yActual = ts(2:end, :);
    for i = 1:p
        out(s, i) = robustPearsonr(yPred(:, i), yActual(:, i));
    end
end
end
