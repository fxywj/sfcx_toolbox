function out = hybridPerRegion(coh)
%HYBRIDPERREGION Per-region hybrid: one weight vector per region.
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

for s = 1:n
    fc = vecToMat(coh.fc(s, :), coh.iu, coh.iv, coh.nNodes);
    for i = 1:p
        obj = @(w) -regionCoupling(w, M, fc, p, i);
        w = fmincon(obj, x0, [], [], ones(1, nMethods), 1, lb, ub, [], options);
        w = max(w, 0);
        w = w / sum(w);
        combinedCol = zeros(p, 1);
        for k = 1:nMethods
            combinedCol = combinedCol + w(k) * M{k}(:, i);
        end
        mask = true(p, 1);
        mask(i) = false;
        out(s, i) = robustPearsonr(combinedCol(mask), fc(mask, i));
    end
end
end

function r = regionCoupling(w, M, fc, p, idx)
w = max(w(:), 0);
w = w / sum(w);
combinedCol = zeros(p, 1);
for k = 1:numel(M)
    combinedCol = combinedCol + w(k) * M{k}(:, idx);
end
mask = true(p, 1);
mask(idx) = false;
r = robustPearsonr(combinedCol(mask), fc(mask, idx));
end
