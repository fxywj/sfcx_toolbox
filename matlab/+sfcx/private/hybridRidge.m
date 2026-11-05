function out = hybridRidge(coh)
%HYBRIDRIDGE Ridge regression hybrid.
M = hybridMatrices(coh);
n = coh.n;
p = coh.nNodes;
out = zeros(n, p);
alphas = logspace(-3, 3, 20);

for s = 1:n
    fc = vecToMat(coh.fc(s, :), coh.iu, coh.iv, coh.nNodes);
    for i = 1:p
        mask = true(p, 1);
        mask(i) = false;
        X = zeros(p - 1, numel(M));
        for k = 1:numel(M)
            X(:, k) = M{k}(mask, i);
        end
        y = fc(mask, i);
        % Ridge with CV over alphas
        bestR = -inf;
        for alpha = alphas
            beta = (X' * X + alpha * eye(size(X, 2))) \ (X' * y);
            pred = X * beta;
            r = robustPearsonr(pred, y);
            if r > bestR
                bestR = r;
            end
        end
        out(s, i) = bestR;
    end
end
end
