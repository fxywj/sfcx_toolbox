function v = regionalVariability(regional, covariates)
%REGIONALVARIABILITY SD across subjects of regional coupling, optionally adjusted.
x = double(regional);
if nargin >= 2 && ~isempty(covariates)
    x = residualise(x, covariates);
end
v = nanStd(x, 1);
v = v(:)';
end
