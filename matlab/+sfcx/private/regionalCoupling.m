function out = regionalCoupling(coh, metric)
%REGIONALCOUPLING Per-subject regional coupling (nSubjects x nNodes).
%   Dispatches to the appropriate private implementation.
switch metric
    case 'reliance_raw'
        out = relianceRaw(coh);
    case 'reliance_deriv'
        out = relianceDeriv(coh);
    case 'reliance_raw_nonlinear'
        out = relianceRawNonlinear(coh);
    case 'reliance_deriv_nonlinear'
        out = relianceDerivNonlinear(coh);
    case 'information_geometry'
        out = informationGeometry(coh);
    case 'multi_scale_transfer_entropy'
        out = multiScaleTransferEntropy(coh);
    case 'criticality'
        out = criticalityIndex(coh);
    case 'granger'
        out = grangerCausality(coh);
    case 'rnc'
        out = relianceNullCoupling(coh);
    case 'trr'
        out = temporalReachReliance(coh);
    case 'hybrid_global'
        out = hybridGlobal(coh);
    case 'hybrid_per_region'
        out = hybridPerRegion(coh);
    case 'hybrid_cv'
        out = hybridCv(coh);
    case 'hybrid_ridge'
        out = hybridRidge(coh);
    case 'hybrid_network'
        out = hybridNetwork(coh);
    otherwise
        error('sfcx:regionalCoupling', 'unknown regional metric ''%s''', metric);
end
end
