function v = variability(edgeField, coh, estimator)
%VARIABILITY Dispersion map for an edge-level coupling field.
%   'edge' : SD across subjects of every edge
%   'node' : mean of each region's edge SDs (node projection)
if nargin < 3 || isempty(estimator)
    estimator = 'edge';
end
edgeField = double(edgeField);
switch lower(estimator)
    case 'edge'
        v = nanStd(edgeField, 1);
        v = v(:)';
    case 'node'
        es = nanStd(edgeField, 1);
        es = es(:);
        v = nan(1, coh.nNodes);
        for i = 1:coh.nNodes
            idx = find((coh.iu == i) | (coh.iv == i));
            if ~isempty(idx)
                v(i) = nanMeanVec(es(idx), 1);
            end
        end
    otherwise
        error('sfcx:variability', 'unknown estimator ''%s''', estimator);
end
end
