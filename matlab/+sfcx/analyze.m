function bundle = analyze(coh, metricsList, varargin)
%ANALYZE Run the full analysis battery and return a serialisable bundle.
%
%   bundle = sfcx.analyze(coh, {'communicability','direct'}, ...
%                         'estimators', {'regional','node','edge'}, ...
%                         'nPerm', 500, 'null', 'auto', 'seed', 0)
%
%   bundle.meta      dataset facts and the null model that was used
%   bundle.metrics   one struct per definition: coupling summary, dispersion
%                         per estimator, axis, layers, reliability, epicentres
%   bundle.sweep     definition x estimator matrix of axis r and p
%   bundle.tables    CSV text for every table
if nargin < 2 || isempty(metricsList)
    metricsList = {'communicability'};
end
if ischar(metricsList)
    metricsList = {metricsList};
end
metricsList = cellfun(@char, metricsList, 'UniformOutput', false);

opts = struct('estimators', {{'regional', 'node', 'edge'}}, ...
              'nPerm', 500, 'null', 'auto', 'seed', 0, ...
              'doSweep', true, 'progress', [], 'state', '');
k = 1;
while k <= numel(varargin)
    name = varargin{k};
    if k + 1 > numel(varargin)
        error('sfcx:analyze', 'missing value for option %s', name);
    end
    if ~isfield(opts, name)
        error('sfcx:analyze', 'unknown option %s', name);
    end
    opts.(name) = varargin{k + 1};
    k = k + 2;
end
if ischar(opts.estimators)
    opts.estimators = {opts.estimators};
end
validEst = {'regional', 'node', 'edge'};
bad = setdiff(opts.estimators, validEst);
if ~isempty(bad)
    error('sfcx:analyze', 'unknown estimator ''%s''', bad{1});
end
if coh.n < 2
    error('sfcx:analyze:tooFewSubjects', ['数据只包含 %d 个被试，无法计算变异度类统计量', ...
        '（SD 在数学上未定义）。单个被试请改用 sfcx.couplingValues 计算耦合值；', ...
        '多被试数据请以 subjects x P x P 三维数组、subjects x E 边向量或 ', ...
        'nSubjects x P x T BOLD 提供。'], coh.n);
end
if coh.n < 8
    warning('sfcx:analyze', ...
        '被试数 %d < 8：变异度、可靠性与个体信号估计噪声较大，结果仅作演示参考。', coh.n);
end

t0 = tic;
if isempty(opts.state)
    state = coh.state;
else
    state = char(opts.state);
end

nullKind = opts.null;
if strcmpi(nullKind, 'auto')
    if isempty(coh.coords)
        nullKind = 'msr';
    else
        nullKind = 'spin';
    end
end
engine = sfcx.nullEngine(coh, nullKind, opts.nPerm, opts.seed);

meta = struct('nSubjects', coh.n, 'nNodes', coh.nNodes, 'nEdges', coh.nEdges, ...
              'nNetworks', coh.nNetworks(), 'state', state, ...
              'hasBehaviour', ~isempty(coh.behaviour), ...
              'hasCovariates', ~isempty(coh.covariates), ...
              'hasCoords', ~isempty(coh.coords), ...
              'nullModel', engine.kind, 'nPerm', opts.nPerm);
settings = struct('metrics', {metricsList}, 'estimators', {opts.estimators}, ...
                  'nPerm', opts.nPerm, 'null', engine.kind, 'seed', opts.seed);

nM = numel(metricsList);
nE = numel(opts.estimators);
sweepR = nan(nM, nE);
sweepP = nan(nM, nE);
blocks = struct([]);

primaryEst = opts.estimators{1};
if any(strcmp(opts.estimators, 'regional'))
    primaryEst = 'regional';
end

for i = 1:nM
    metric = metricsList{i};
    report(opts.progress, (i - 1) / nM, sprintf('coupling: %s', metric));
    c = sfcx.coupling(coh, metric);
    if strcmp(c.kind, 'regional')
        regional = c.field;
        edgeField = sfcx.edgeProjection(regional, coh);
    else
        regional = sfcx.regionalReadout(coh, metric, c.field);
        edgeField = c.field;
    end

    info = metricMeta(metric);
    block = struct();
    block.name = metric;
    block.family = info.family;
    block.label = info.label;
    block.description = info.description;
    block.fieldKind = c.kind;
    block.couplingSummary = struct( ...
        'mean', nanMeanVec(c.field(:), 1), ...
        'sd', nanStd(c.field(:), 1), ...
        'min', min(c.field(:), [], 'omitnan'), ...
        'max', max(c.field(:), [], 'omitnan'), ...
        'nanFraction', mean(~isfinite(c.field(:))));

    estItems = struct([]);
    for j = 1:nE
        est = opts.estimators{j};
        [v, level] = estVariability(coh, est, regional, edgeField);
        control = estMagnitude(coh, est, regional, edgeField);
        ax = axisBlock(coh, v, engine, control);
        sweepR(i, j) = ax.r;
        sweepP(i, j) = ax.p;
        item = struct();
        item.name = est;
        item.level = level;
        item.n = numel(v);
        item.min = min(v, [], 'omitnan');
        item.max = max(v, [], 'omitnan');
        item.median = median(v(isfinite(v)));
        item.mean = nanMeanVec(v(:), 1);
        item.axis = ax;
        item.values = v(:)';
        estItems = appendStruct(estItems, item);
    end
    block.estimators = estItems;
    block.primaryEstimator = primaryEst;

    prim = estItems(strcmp({estItems.name}, primaryEst));
    primV = prim.values;

    block.epicentres = emptyEpicentres();
    if strcmp(prim.level, 'regional')
        block.epicentres = sfcx.epicentres(primV, coh, 2, opts.seed, min(opts.nPerm, 1000));
    end

    block.layers = sfcx.layers(sfcx.variability(edgeField, coh, 'edge'), coh, ...
                               min(opts.nPerm, 1000), opts.seed);
    block.reliability = sfcx.reliability(edgeField, coh, opts.seed);

    dev = sfcx.deviationSpectrum(edgeField, coh, ...
                                 {'regional', 'network', 'edge', 'gap'}, 20000, opts.seed);
    scales = fieldnames(dev);
    indiv = struct();
    bestScale = '';
    bestP = Inf;
    for s = 1:numel(scales)
        sc = scales{s};
        entry = struct('mean', nanMeanVec(dev.(sc), 1), ...
                       'sd', nanStd(dev.(sc), 1));
        if ~isempty(coh.behaviour)
            ps = sfcx.partialSpearman(dev.(sc), coh.behaviour, coh.covariates);
            entry.r = ps.r;
            entry.p = ps.p;
            entry.n = ps.n;
            if isfinite(ps.p) && ps.p < bestP
                bestP = ps.p;
                bestScale = sc;
            end
        end
        indiv.(sc) = entry;
    end
    if ~isempty(bestScale)
        indiv.bestScale = bestScale;
    end
    block.individual = indiv;

    blocks = appendStruct(blocks, block);
    report(opts.progress, i / nM, sprintf('done: %s', metric));
end

sweep = struct();
sweep.rows = metricsList;
sweep.cols = opts.estimators;
sweep.r = sweepR;
sweep.p = sweepP;
sweep.done = logical(opts.doSweep);

tables = struct();
for i = 1:nM
    name = metricsList{i};
    v = blocks(i).estimators(strcmp({blocks(i).estimators.name}, primaryEst)).values;
    tables.(sprintf('%s_variability', name)) = ...
        sprintf('dispersion\n%s', sprintf('%.10g\n', v));
end
tables.axis_r = matrixCsv(sweepR, metricsList, opts.estimators);
tables.axis_p = matrixCsv(sweepP, metricsList, opts.estimators);
tables.axis_summary = axisSummaryCsv(sweepR, sweepP, metricsList, opts.estimators);

bundle = struct('meta', meta, 'settings', settings, 'metrics', blocks, ...
                'sweep', sweep, 'tables', tables, 'runtimeS', toc(t0));
end

% ---------------------------------------------------------------- helpers --

function [v, level] = estVariability(coh, est, regional, edgeField)
switch est
    case 'regional'
        v = sfcx.regionalVariability(regional, coh.covariates);
        level = 'regional';
    case 'node'
        v = sfcx.variability(edgeField, coh, 'node');
        level = 'regional';
    otherwise
        v = sfcx.variability(edgeField, coh, 'edge');
        level = 'edge';
end
end

function m = estMagnitude(coh, est, regional, edgeField)
switch est
    case 'regional'
        m = nanMeanVec(regional, 1);
    case 'node'
        me = nanMeanVec(edgeField, 1);
        m = nan(1, coh.nNodes);
        for i = 1:coh.nNodes
            idx = find((coh.iu == i) | (coh.iv == i));
            if ~isempty(idx)
                m(i) = nanMeanVec(me(idx), 1);
            end
        end
    otherwise
        m = nanMeanVec(edgeField, 1);
end
end

function ax = axisBlock(coh, v, engine, control)
if numel(v) == coh.nNodes
    h = coh.hierarchy;
else
    h = 0.5 * (coh.hierarchy(coh.iu) + coh.hierarchy(coh.iv));
end
ax = sfcx.axis(v, h, engine);
ax.magnitude = nanMeanVec(v(:), 1);
ax.controlledForMagnitude = struct('r', NaN, 'p', NaN, 'null', 'none', ...
                                   'nPerm', 0, 'n', 0);
if numel(v) == coh.nNodes && ~isempty(control)
    ax.controlledForMagnitude = sfcx.axisControlled(v, coh.hierarchy, control, engine);
end
end

function s = appendStruct(s, item)
if isempty(s)
    s = item;
else
    s(end + 1) = item;
end
end

function e = emptyEpicentres()
e = struct('regions', [], 'n', 0, 'internalEdges', 0, 'z', NaN, 'p', NaN, ...
           'nullMean', NaN, 'clusterLabels', [], 'networkCounts', []);
end

function m = metricMeta(metric)
t = metricCatalog();
idx = find(strcmp({t.id}, metric), 1);
if ~isempty(idx)
    m = struct('family', t(idx).family, 'label', t(idx).label, ...
               'description', t(idx).description);
    return
end
e = registry('get', metric);
if isempty(e)
    m = struct('family', 'custom', 'label', metric, 'description', '');
else
    m = struct('family', e.family, 'label', e.label, 'description', e.description);
end
end

function report(fn, frac, msg)
if ~isempty(fn)
    fn(min(max(frac, 0), 1), msg);
end
end

function txt = matrixCsv(m, rowNames, colNames)
lines = {strjoin(['metric', colNames], ',')};
for i = 1:numel(rowNames)
    parts = rowNames{i};
    for j = 1:numel(colNames)
        if isfinite(m(i, j))
            parts = [parts ',' sprintf('%.6f', m(i, j))];
        else
            parts = [parts ','];
        end
    end
    lines{end + 1} = parts;
end
txt = strjoin(lines, newline);
end

function txt = axisSummaryCsv(sweepR, sweepP, rowNames, colNames)
lines = {'metric,estimator,axis_r,axis_p'};
for i = 1:numel(rowNames)
    for j = 1:numel(colNames)
        r = sweepR(i, j);
        p = sweepP(i, j);
        if isfinite(r)
            rs = sprintf('%.6f', r);
        else
            rs = '';
        end
        if isfinite(p)
            ps = sprintf('%.6g', p);
        else
            ps = '';
        end
        lines{end + 1} = [rowNames{i} ',' colNames{j} ',' rs ',' ps];
    end
end
txt = strjoin(lines, newline);
end
