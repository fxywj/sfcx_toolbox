function a = structural(coh, binarise)
%STRUCTURAL Group-mean structural matrix, optionally binarised.
if nargin < 2
    binarise = true;
end
meanSc = mean(coh.sc, 1, 'omitnan');
a = vecToMat(meanSc, coh.iu, coh.iv, coh.nNodes);
if binarise
    a = double(a > 1e-9);
end
end
