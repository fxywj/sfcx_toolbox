function out = residualise(y, cov)
%RESIDUALISE Remove a linear covariate model from each column of y.
if isempty(cov)
    out = y;
    return;
end
if isvector(cov)
    cov = cov(:);
end
x = [ones(size(cov, 1), 1) cov];
beta = x \ y;
out = y - x * beta;
end
