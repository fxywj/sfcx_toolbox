function out = coupleVector(measureVec, fc, iu, nNodes)
%COUPLEVECTOR Edge-wise product of two standardised vectors, per subject.
mg = zscoreByGroup(reshape(measureVec, 1, []), iu, nNodes);
fg = zscoreByGroup(double(fc), iu, nNodes);
out = fg .* mg;
end
