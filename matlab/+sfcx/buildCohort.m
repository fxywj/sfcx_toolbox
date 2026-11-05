function [coh, rep] = buildCohort(cfg)
%BUILDCOHORT Assemble a Cohort from any mix of connectivity / BOLD inputs.
%
%   cfg.sc, cfg.fc   : (subjects, nodes, nodes) | (subjects, edges) | (nodes, nodes)
%   cfg.bold         : (subjects, nodes, time) | (nodes, time)
%   cfg.networks, cfg.hierarchy, cfg.coords, cfg.behaviour, cfg.covariates
%   cfg.fcKind       : 'pearson' (default) | 'precision'
%
%   Missing annotation is defaulted, and every default is reported so the
%   caller can surface honest warnings.
rep = struct('inferred', {{}}, 'warnings', {{}});
fields = {'sc', 'fc', 'bold', 'networks', 'hierarchy', 'coords', ...
          'behaviour', 'covariates', 'state', 'subjectIds', 'fcKind', 'timeAxis'};
for i = 1:numel(fields)
    if ~isfield(cfg, fields{i})
        cfg.(fields{i}) = [];
    end
end
if isempty(cfg.sc) && isempty(cfg.bold)
    error('sfcx:build', 'need at least a structural matrix (sc) or BOLD series (bold)');
end

nNodes = [];
scV = [];
fcV = [];
if ~isempty(cfg.sc)
    [scV, nNodes] = toEdgeVectors(cfg.sc, 'sc');
    rep.inferred{end + 1} = sprintf('SC: %d subjects x %d edges (P=%d)', ...
        size(scV, 1), size(scV, 2), nNodes);
end
if ~isempty(cfg.bold)
    if isempty(cfg.fcKind)
        fcKind = 'pearson';
    else
        fcKind = cfg.fcKind;
    end
    if isempty(cfg.timeAxis)
        timeAxis = 'auto';
    else
        timeAxis = cfg.timeAxis;
    end
    fcV = sfcx.boldToFc(cfg.bold, fcKind, timeAxis);
    pBold = round((1 + sqrt(1 + 8 * size(fcV, 2))) / 2);
    if ~isempty(nNodes) && pBold ~= nNodes
        error('sfcx:build', 'SC has %d nodes but BOLD has %d', nNodes, pBold);
    end
    nNodes = pBold;
    rep.inferred{end + 1} = sprintf('BOLD -> FC (%s): %d subjects', fcKind, size(fcV, 1));
end
if ~isempty(cfg.fc)
    [fcV, pFc] = toEdgeVectors(cfg.fc, 'fc');
    if ~isempty(nNodes) && pFc ~= nNodes
        error('sfcx:build', 'SC has %d nodes but FC has %d', nNodes, pFc);
    end
    nNodes = pFc;
elseif isempty(fcV)
    error('sfcx:build', 'need an FC matrix or a BOLD series');
end
if size(scV, 1) ~= size(fcV, 1)
    error('sfcx:build', 'SC has %d subjects but FC has %d', size(scV, 1), size(fcV, 1));
end

if isempty(cfg.networks)
    networks = zeros(1, nNodes);
    rep.warnings{end + 1} = ...
        'no network labels supplied: one community generated; the within/between contrast is meaningless';
else
    networks = double(cfg.networks(:)');
    if numel(networks) ~= nNodes
        error('sfcx:build', 'networks must have %d entries', nNodes);
    end
end

if isempty(cfg.hierarchy)
    [iu, iv] = upperTriIndex(nNodes);
    absmean = mean(abs(double(scV)), 1);
    strength = accumarray(iu, absmean(:), [nNodes 1]) + accumarray(iv, absmean(:), [nNodes 1]);
    strength = strength';
    if all(abs(strength - strength(1)) < 1e-12)
        hierarchy = linspace(1, 7, nNodes);
    else
        order = rankVector(strength) - 1;
        hierarchy = 1 + 6 * order / max(nNodes - 1, 1);
        hierarchy = hierarchy(:)';
    end
    rep.warnings{end + 1} = ...
        'no cortical hierarchy supplied: node structural strength used as a proxy';
else
    hierarchy = double(cfg.hierarchy(:)');
    if numel(hierarchy) ~= nNodes
        error('sfcx:build', 'hierarchy must have %d entries', nNodes);
    end
end

if ~isempty(cfg.coords) && ~isequal(size(cfg.coords), [nNodes, 3])
    error('sfcx:build', 'coords must be (nNodes, 3)');
end
if ~isempty(cfg.behaviour) && numel(cfg.behaviour) ~= size(scV, 1)
    error('sfcx:build', 'behaviour must have %d entries', size(scV, 1));
end
if ~isempty(cfg.covariates) && size(cfg.covariates, 1) ~= size(scV, 1)
    error('sfcx:build', 'covariates must have %d rows', size(scV, 1));
end

opts = struct();
opts.coords = cfg.coords;
opts.behaviour = cfg.behaviour;
opts.covariates = cfg.covariates;
if ~isempty(cfg.state)
    opts.state = cfg.state;
end
if ~isempty(cfg.subjectIds)
    opts.subjectIds = cfg.subjectIds;
end
if ~isempty(cfg.bold)
    opts.bold = cfg.bold;
end
coh = sfcx.Cohort(scV, fcV, nNodes, networks, hierarchy, opts);
coh.meta = struct('report', rep);
end

function [v, p] = toEdgeVectors(x, name)
x = double(x);
if ndims(x) == 3
    n = size(x, 1);
    p = size(x, 2);
    if size(x, 3) ~= p
        error('sfcx:build', '%s: last two dimensions must be square', name);
    end
    [iu, iv] = upperTriIndex(p);
    X = reshape(x, n, p * p);
    v = X(:, sub2ind([p p], iu, iv));
    return
end
if ndims(x) ~= 2
    error('sfcx:build', '%s must be 2-D or 3-D', name);
end
if size(x, 1) == size(x, 2) && size(x, 1) > 1
    p = size(x, 1);
    [iu, iv] = upperTriIndex(p);
    v = reshape(x(sub2ind([p p], iu, iv)), 1, []);
    return
end
e = size(x, 2);
p = round((1 + sqrt(1 + 8 * e)) / 2);
if p * (p - 1) / 2 ~= e
    error('sfcx:build', '%s: %d columns is not an upper triangle', name, e);
end
v = x;
end
