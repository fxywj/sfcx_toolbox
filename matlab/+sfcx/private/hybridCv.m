function out = hybridCv(coh, nSplits)
%HYBRIDCV Cross-validated global hybrid.
if nargin < 2 || isempty(nSplits)
    nSplits = 5;
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
regionIdx = randperm(p);
foldSize = floor(p / nSplits);

for s = 1:n
    fc = vecToMat(coh.fc(s, :), coh.iu, coh.iv, coh.nNodes);
    for fold = 1:nSplits
        testIdx = regionIdx((fold - 1) * foldSize + 1:min(fold * foldSize, p));
        trainIdx = setdiff(regionIdx, testIdx);
        obj = @(w) -foldCoupling(w, M, fc, p, trainIdx);
        w = fmincon(obj, x0, [], [], ones(1, nMethods), 1, lb, ub, [], options);
        w = max(w, 0);
        w = w / sum(w);
        combined = zeros(p, p);
        for k = 1:nMethods
            combined = combined + w(k) * M{k};
        end
        for i = testIdx
            mask = true(p, 1);
            mask(i) = false;
            out(s, i) = robustPearsonr(combined(mask, i), fc(mask, i));
        end
    end
end
end

function r = foldCoupling(w, M, fc, p, trainIdx)
w = max(w(:), 0);
w = w / sum(w);
combined = zeros(p, p);
for k = 1:numel(M)
    combined = combined + w(k) * M{k};
end
r = 0;
for i = trainIdx
    mask = true(p, 1);
    mask(i) = false;
    mask(~ismember(1:p, trainIdx)) = false;
    r = r + robustPearsonr(combined(mask, i), fc(mask, i));
end
r = r / max(numel(trainIdx), 1);
end
