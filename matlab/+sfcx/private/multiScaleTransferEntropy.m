function out = multiScaleTransferEntropy(coh, scales, nCondBins)
%MULTISCALETRANSFERENTROPY Multi-scale transfer entropy from SC to BOLD.
if isempty(coh.bold)
    error('sfcx:multiScaleTransferEntropy', 'requires BOLD time series');
end
if nargin < 2 || isempty(scales)
    scales = [1, 2, 4, 8];
end
if nargin < 3 || isempty(nCondBins)
    nCondBins = 5;
end
n = coh.n;
p = coh.nNodes;
W = colNormalize(coh.groupMatrix('sc'));
out = zeros(n, p);

for s = 1:n
    ts = squeeze(coh.bold(s, :, :))';  % (T, nNodes)
    teScores = [];
    for scale = scales
        tsDown = ts(1:scale:end, :);
        T = size(tsDown, 1);
        if T < 20
            teScores(end+1, :) = zeros(1, p); %#ok<AGROW>
            continue;
        end
        scPred = tsDown(1:end-1, :) * W;
        te = zeros(1, p);
        for i = 1:p
            scVals = scPred(:, i);
            fcFuture = tsDown(2:end, i);
            fcPast = tsDown(1:end-1, i);
            nGroups = min(nCondBins, floor(numel(fcPast) / 20));
            if nGroups < 2
                te(i) = 0;
                continue;
            end
            quantiles = prctile(fcPast, linspace(0, 100, nGroups + 1));
            weightedMi = 0;
            totalWeight = 0;
            for g = 1:nGroups
                if g == nGroups
                    mask = (fcPast >= quantiles(g)) & (fcPast <= quantiles(g + 1));
                else
                    mask = (fcPast >= quantiles(g)) & (fcPast < quantiles(g + 1));
                end
                if sum(mask) < 10
                    continue;
                end
                scG = scVals(mask);
                fcFG = fcFuture(mask);
                scBin = discretize(scG, linspace(min(scG), max(scG), 15));
                fcFBin = discretize(fcFG, linspace(min(fcFG), max(fcFG), 15));
                miG = mutualInfo(scBin, fcFBin);
                weightedMi = weightedMi + miG * sum(mask);
                totalWeight = totalWeight + sum(mask);
            end
            te(i) = weightedMi / max(totalWeight, 1);
        end
        teScores(end+1, :) = te; %#ok<AGROW>
    end
    out(s, :) = mean(teScores, 1);
end
end

function mi = mutualInfo(x, y)
%MUTUALINFO Discrete mutual information.
x = x(:);
y = y(:);
joint = accumarray([x, y], 1, [15, 15]);
px = sum(joint, 2) / sum(joint(:));
py = sum(joint, 1) / sum(joint(:));
pxy = joint / sum(joint(:));
mi = 0;
for i = 1:15
    for j = 1:15
        if pxy(i, j) > 0 && px(i) > 0 && py(j) > 0
            mi = mi + pxy(i, j) * log(pxy(i, j) / (px(i) * py(j)));
        end
    end
end
end
