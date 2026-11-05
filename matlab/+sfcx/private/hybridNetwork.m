function out = hybridNetwork(coh, nClusters, nSplits)
%HYBRIDNETWORK Network-constrained hybrid with K-means clusters.
if nargin < 2 || isempty(nClusters)
    nClusters = 7;
end
if nargin < 3 || isempty(nSplits)
    nSplits = 3;
end
M = hybridMatrices(coh);
nMethods = numel(M);
n = coh.n;
p = coh.nNodes;
out = zeros(n, p);
lb = zeros(nMethods, 1);
ub = ones(nMethods, 1);
x0 = ones(nMethods, 1) / nMethods;
options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 500, ...
    'OptimalityTolerance', 1e-10, 'Algorithm', 'sqp');

rng(42, 'twister');

for s = 1:n
    fc = vecToMat(coh.fc(s, :), coh.iu, coh.iv, coh.nNodes);
    regionIdx = randperm(p);
    foldSize = floor(p / nSplits);
    for fold = 1:nSplits
        testIdx = regionIdx((fold - 1) * foldSize + 1:min(fold * foldSize, p));
        trainIdx = setdiff(regionIdx, testIdx);
        fcTrain = fc(trainIdx, trainIdx);
        kmeansIdx = kmeans(fcTrain, nClusters, 'Replicates', 10, 'Display', 'off');
        uniqueLabels = unique(kmeansIdx);
        labelMap = containers.Map(uniqueLabels, 1:numel(uniqueLabels));
        trainLabels = arrayfun(@(x) labelMap(x), kmeansIdx);
        % Predict test labels by nearest centroid
        centroids = zeros(numel(uniqueLabels), size(fcTrain, 2));
        for c = 1:numel(uniqueLabels)
            centroids(c, :) = mean(fcTrain(trainLabels == c, :), 1);
        end
        testLabels = zeros(size(testIdx));
        for t = 1:numel(testIdx)
            dists = sum((fc(testIdx(t), trainIdx) - centroids) .^ 2, 2);
            [~, testLabels(t)] = min(dists);
        end
        for net = 1:numel(uniqueLabels)
            netTrain = trainIdx(trainLabels == net);
            netTest = testIdx(testLabels == net);
            if numel(netTrain) < 2 || isempty(netTest)
                continue;
            end
            obj = @(w) -netCoupling(w, M, fc, p, netTrain);
            w = fmincon(obj, x0, [], [], ones(1, nMethods), 1, lb, ub, [], options);
            w = max(w, 0);
            w = w / sum(w);
            combined = zeros(p, p);
            for k = 1:nMethods
                combined = combined + w(k) * M{k};
            end
            for i = netTest
                mask = true(p, 1);
                mask(i) = false;
                out(s, i) = robustPearsonr(combined(mask, i), fc(mask, i));
            end
        end
    end
end
end

function r = netCoupling(w, M, fc, p, netTrain)
w = max(w(:), 0);
w = w / sum(w);
combined = zeros(p, p);
for k = 1:numel(M)
    combined = combined + w(k) * M{k};
end
r = 0;
for i = netTrain
    mask = true(p, 1);
    mask(i) = false;
    mask(~ismember(1:p, netTrain)) = false;
    r = r + robustPearsonr(combined(mask, i), fc(mask, i));
end
r = r / max(numel(netTrain), 1);
end
