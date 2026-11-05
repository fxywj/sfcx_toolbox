function out = coupling(coh, metric)
%COUPLING One coupling definition, layer 1 of the coordinate system.
%   out = sfcx.coupling(coh, 'communicability')
%   out.kind   'edge' (nSubjects x nEdges) or 'regional' (nSubjects x nNodes)
%   out.field  the coupling field itself
if nargin < 2 || isempty(metric)
    metric = 'communicability';
end
metric = char(metric);

e = registry('get', metric);
if ~isempty(e)
    [kind, field] = evaluateCustom(e.fn, coh);
    out = struct('kind', kind, 'field', field);
    return
end

switch metric
    case 'direct'
        sc = zscoreByGroup(double(coh.sc), coh.iu, coh.nNodes);
        fc = zscoreByGroup(double(coh.fc), coh.iu, coh.nNodes);
        out = struct('kind', 'edge', 'field', sc .* fc);
    case 'rank_direct'
        sc = zeros(size(coh.sc));
        fc = zeros(size(coh.fc));
        for s = 1:coh.n
            sc(s, :) = rankVector(coh.sc(s, :));
            fc(s, :) = rankVector(coh.fc(s, :));
        end
        sc = zscoreByGroup(sc, coh.iu, coh.nNodes);
        fc = zscoreByGroup(fc, coh.iu, coh.nNodes);
        out = struct('kind', 'edge', 'field', sc .* fc);
    case 'reliance'
        out = struct('kind', 'edge', 'field', relianceField(coh));
    case {'cosine', 'matching', 'communicability', 'degree_product', ...
          'mfpt', 'path_ease', 'diffusion_efficiency', 'search_info', 'flow'}
        a = structural(coh, true);
        switch metric
            case 'cosine'
                m = cosineSimilarity(a);
            case 'matching'
                m = matchingIndex(a);
            case 'communicability'
                m = communicability(a);
            case 'degree_product'
                m = degreeProduct(a);
            case 'mfpt'
                m = meanFirstPassage(a);
            case 'path_ease'
                m = pathEase(a);
            case 'diffusion_efficiency'
                m = diffusionEfficiency(a);
            case 'search_info'
                m = searchInfo(a);
            case 'flow'
                m = flowTwoStep(a);
        end
        m = pairScale(m);
        v = matToVec(m, coh.iu, coh.iv);
        out = struct('kind', 'edge', ...
                     'field', coupleVector(v, coh.fc, coh.iu, coh.nNodes));
    case {'reliance_raw', 'reliance_deriv', 'reliance_raw_nonlinear', ...
          'reliance_deriv_nonlinear', 'information_geometry', ...
          'multi_scale_transfer_entropy', 'criticality', 'granger', ...
          'rnc', 'trr', 'hybrid_global', 'hybrid_per_region', ...
          'hybrid_cv', 'hybrid_ridge', 'hybrid_network'}
        out = struct('kind', 'regional', 'field', regionalCoupling(coh, metric));
    otherwise
        error('sfcx:coupling', 'unknown metric ''%s''; run sfcx.catalog() to list', metric);
end
end
