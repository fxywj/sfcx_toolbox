function out = relianceRaw(coh)
%RELIANCERAW Classical reliance: corr(ts_{t-1} W, ts_t) per region.
if isempty(coh.bold)
    error('sfcx:relianceRaw', 'requires BOLD time series in the cohort');
end
n = coh.n;
p = coh.nNodes;
W = colNormalize(coh.groupMatrix('sc'));
out = zeros(n, p);
for s = 1:n
    ts = squeeze(coh.bold(s, :, :))';  % (T, nNodes)
    yPred = ts(1:end-1, :) * W;
    yActual = ts(2:end, :);
    for i = 1:p
        out(s, i) = robustPearsonr(yPred(:, i), yActual(:, i));
    end
end
end
