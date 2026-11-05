function out = grangerCausality(coh, maxLag)
%GRANGERCAUSALITY Granger causality F-statistic: SC-weighted past predicts BOLD future.
if isempty(coh.bold)
    error('sfcx:grangerCausality', 'requires BOLD time series');
end
if nargin < 2 || isempty(maxLag)
    maxLag = 3;
end
n = coh.n;
p = coh.nNodes;
W = colNormalize(coh.groupMatrix('sc'));
out = zeros(n, p);

for s = 1:n
    ts = squeeze(coh.bold(s, :, :))';  % (T, nNodes)
    T = size(ts, 1);
    if T <= maxLag + 2
        continue;
    end
    scWeightedFull = ts(1:end-1, :) * W;
    for i = 1:p
        y = ts(maxLag+1:end, i);
        XRestricted = zeros(T - maxLag, maxLag);
        for lag = 1:maxLag
            XRestricted(:, lag) = ts(maxLag - lag + 1:T - lag, i);
        end
        scCols = cell(maxLag, 1);
        for lag = 1:maxLag
            startIdx = maxLag - lag + 1;
            endIdx = T - lag - 1;
            if startIdx < endIdx
                scCols{lag} = scWeightedFull(startIdx:endIdx, i);
            end
        end
        valid = ~cellfun(@isempty, scCols);
        if ~any(valid)
            continue;
        end
        scCols = scCols(valid);
        minLen = min(cellfun(@numel, scCols));
        scCols = cellfun(@(c) c(end-minLen+1:end), scCols, 'UniformOutput', false);
        XRestrictedAdj = XRestricted(end-minLen+1:end, :);
        yAdj = y(end-minLen+1:end);
        XFull = [XRestrictedAdj, horzcat(scCols{:})];
        try
            betaR = XRestrictedAdj \ yAdj;
            betaF = XFull \ yAdj;
            rssR = sum((yAdj - XRestrictedAdj * betaR) .^ 2);
            rssF = sum((yAdj - XFull * betaF) .^ 2);
            nObs = numel(yAdj);
            q = size(XFull, 2);
            if rssF > 1e-10
                F = ((rssR - rssF) / maxLag) / (rssF / (nObs - q - 1));
                out(s, i) = max(F, 0);
            end
        catch
            out(s, i) = 0;
        end
    end
end
end
