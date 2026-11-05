function out = hybridGlobal(coh)
%HYBRIDGLOBAL Global hybrid: one set of weights across all regions.
M = hybridMatrices(coh);
nMethods = numel(M);
n = coh.n;
p = coh.nNodes;
out = zeros(n, p);
lb = zeros(nMethods, 1);
ub = ones(nMethods, 1);
x0 = ones(nMethods, 1) / nMethods;
options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 1000, ...
    'OptimalityTolerance', 1e-12, 'Algorithm', 'sqp');

for s = 1:n
    fc = vecToMat(coh.fc(s, :), coh.iu, coh.iv, coh.nNodes);
    obj = @(w) -meanCoupling(w, M, fc, p);
    w = fmincon(obj, x0, [], [], ones(1, nMethods), 1, lb, ub, [], options);
    w = max(w, 0);
    w = w / sum(w);
    combined = zeros(p, p);
    for k = 1:nMethods
        combined = combined + w(k) * M{k};
    end
    for i = 1:p
        mask = true(p, 1);
        mask(i) = false;
        out(s, i) = robustPearsonr(combined(mask, i), fc(mask, i));
    end
end
end

function r = meanCoupling(w, M, fc, p)
w = max(w(:), 0);
w = w / sum(w);
combined = zeros(p, p);
for k = 1:numel(M)
    combined = combined + w(k) * M{k};
end
r = 0;
for i = 1:p
    mask = true(p, 1);
    mask(i) = false;
    r = r + robustPearsonr(combined(mask, i), fc(mask, i));
end
r = r / p;
end
