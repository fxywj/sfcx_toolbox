function out = relianceNullCoupling(coh)
%RELIANCENULLCOUPLING RNC: reliance minus weight-shuffle null.
if isempty(coh.bold)
    error('sfcx:relianceNullCoupling', 'requires BOLD time series');
end
n = coh.n;
p = coh.nNodes;
a = coh.groupMatrix('sc');
aLog = log1p(max(a, 0));
W = colNormalize(aLog);
out = zeros(n, p);

for s = 1:n
    ts = squeeze(coh.bold(s, :, :))';  % (T, nNodes)
    yPred = ts(1:end-1, :) * W;
    yActual = ts(2:end, :);
    rObs = zeros(1, p);
    for i = 1:p
        rObs(i) = robustPearsonr(yPred(:, i), yActual(:, i));
    end
    seed = mod(sum(double(char(coh.state))) + s, 2^31);
    aNull = weightShuffleNull(aLog, seed);
    WNull = colNormalize(aNull);
    yPredNull = ts(1:end-1, :) * WNull;
    rNull = zeros(1, p);
    for i = 1:p
        rNull(i) = robustPearsonr(yPredNull(:, i), yActual(:, i));
    end
    out(s, :) = rObs - rNull;
end
end
