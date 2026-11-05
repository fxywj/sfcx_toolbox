function out = informationGeometry(coh)
%INFORMATIONGEOMETRY Log-Euclidean geodesic coupling on the SPD manifold.
n = coh.n;
p = coh.nNodes;
out = zeros(n, p);
for s = 1:n
    sc = vecToMat(coh.sc(s, :), coh.iu, coh.iv, coh.nNodes);
    fc = vecToMat(coh.fc(s, :), coh.iu, coh.iv, coh.nNodes);
    minEigSc = min(eig(sc));
    if minEigSc <= 0
        sc = sc + (abs(minEigSc) + 1e-6) * eye(p);
    end
    minEigFc = min(eig(fc));
    if minEigFc <= 0
        fc = fc + (abs(minEigFc) + 1e-6) * eye(p);
    end
    logSc = logm(sc);
    logFc = logm(fc);
    for i = 1:p
        mask = true(p, 1);
        mask(i) = false;
        scCol = real(logSc(mask, i));
        fcCol = real(logFc(mask, i));
        out(s, i) = robustPearsonr(scCol, fcCol);
    end
end
end
