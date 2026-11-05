function out = temporalReachReliance(coh, lags)
%TEMPORALREACHRELIANCETRR: weighted sum of multi-lag reliance (1/k weights).
if isempty(coh.bold)
    error('sfcx:temporalReachReliance', 'requires BOLD time series');
end
if nargin < 2 || isempty(lags)
    lags = 1:5;
end
n = coh.n;
p = coh.nNodes;
W = colNormalize(coh.groupMatrix('sc'));
w = 1.0 ./ double(lags);
w = w / sum(w);
out = zeros(n, p);

for s = 1:n
    ts = squeeze(coh.bold(s, :, :))';  % (T, nNodes)
    P = ts * W;
    R = nan(numel(lags), p);
    for row = 1:numel(lags)
        k = lags(row);
        if size(ts, 1) - k < 10
            continue;
        end
        for i = 1:p
            R(row, i) = robustPearsonr(P(1:end-k, i), ts(k+1:end, i));
        end
    end
    out(s, :) = nansum(w' .* R, 1);
end
end
