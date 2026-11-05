function out = criticalityIndex(coh)
%CRITICALITYINDEX Distance of SC eigenvalue spectrum from lambda=1.
n = coh.n;
p = coh.nNodes;
out = zeros(n, p);
for s = 1:n
    sc = vecToMat(coh.sc(s, :), coh.iu, coh.iv, coh.nNodes);
    fc = vecToMat(coh.fc(s, :), coh.iu, coh.iv, coh.nNodes);
    scCoupling = sc ./ (sum(sc, 1) + 1e-10);
    eigvals = eig(scCoupling);
    lambdaMax = max(abs(eigvals));
    ci = 1.0 / (1.0 + abs(lambdaMax - 1.0));
    for i = 1:p
        mask = true(p, 1);
        mask(i) = false;
        scRow = scCoupling(i, mask);
        fcCol = fc(mask, i);
        r = robustPearsonr(scRow(:), fcCol);
        out(s, i) = r * ci;
    end
end
end
