function out = relianceDeriv(coh)
%RELIANCEDERIV Derivative reliance: corr(-ts + ts*W, dts/dt) per region.
if isempty(coh.bold)
    error('sfcx:relianceDeriv', 'requires BOLD time series in the cohort');
end
n = coh.n;
p = coh.nNodes;
W = colNormalize(coh.groupMatrix('sc'));
out = zeros(n, p);
for s = 1:n
    ts = squeeze(coh.bold(s, :, :))';  % (T, nNodes)
    dy = ts(2:end, :) - ts(1:end-1, :);
    yCurr = ts(1:end-1, :);
    dyPred = -yCurr + yCurr * W;
    for i = 1:p
        out(s, i) = robustPearsonr(dyPred(:, i), dy(:, i));
    end
end
end
