function mats = hybridMatrices(coh)
%HYBRIDMATRICES Z-scored group-mean SC + communication matrices for hybrid models.
a = coh.groupMatrix('sc');
aLog = log1p(max(a, 0));
mats = cell(8, 1);
mats{1} = zscoreCols(aLog);
aBin = double(aLog > 1e-9);
mats{2} = zscoreCols(cosineSimilarity(aBin));
mats{3} = zscoreCols(matchingIndex(aBin));
mats{4} = zscoreCols(communicability(aBin));
mats{5} = zscoreCols(meanFirstPassage(aBin));
mats{6} = zscoreCols(pathEase(aBin));
mats{7} = zscoreCols(searchInfo(aBin));
mats{8} = zscoreCols(flowTwoStep(aBin));
end
