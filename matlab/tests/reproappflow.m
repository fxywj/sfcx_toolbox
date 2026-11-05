% Simulate the app flow exactly: load demo5_group.mat -> buildCohort -> analyze.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);

mat = fullfile(fileparts(root), 'examples', 'data', 'demo5_group.mat');
fprintf('data file: %s (exists=%d)\n', mat, exist(mat, 'file'));

s = load(mat);
disp(fieldnames(s));

% what pickField would return for each slot
S = struct();
S.sc = s.sc; S.fc = s.fc; S.bold = [];
S.networks = []; S.hierarchy = []; S.coords = [];
S.behaviour = []; S.covariates = [];
if isfield(s, 'networks'), S.networks = s.networks; end
if isfield(s, 'hierarchy'), S.hierarchy = s.hierarchy; end
if isfield(s, 'coords'), S.coords = s.coords; end
if isfield(s, 'behaviour'), S.behaviour = s.behaviour; end

fprintf('sc size: %s\n', mat2str(size(S.sc)));
fprintf('fc size: %s\n', mat2str(size(S.fc)));

cfg = struct('sc', S.sc, 'fc', S.fc, 'bold', S.bold, ...
             'networks', S.networks, 'hierarchy', S.hierarchy, ...
             'coords', S.coords, 'behaviour', S.behaviour, ...
             'covariates', S.covariates, ...
             'fcKind', 'pearson', 'state', 'rest');
try
    coh = sfcx.buildCohort(cfg);
    fprintf('cohort: n=%d nodes=%d edges=%d\n', coh.n, coh.nNodes, coh.nEdges);
catch err
    fprintf('BUILDCOHORT FAILED: %s\n', err.message);
    exit(1);
end

try
    t0 = tic;
    bundle = sfcx.analyze(coh, {'communicability'}, ...
                          'estimators', {'regional', 'node', 'edge'}, ...
                          'nPerm', 100, 'null', 'auto', 'seed', 0, ...
                          'state', 'rest');
    fprintf('analyze ok in %.1fs; metrics=%d\n', toc(t0), numel(bundle.metrics));
    e = bundle.metrics(1).estimators(1);
    fprintf('first axis: r=%.3f p=%.3g null=%s\n', e.axis.r, e.axis.p, e.axis.null);
catch err
    fprintf('ANALYZE FAILED: %s\n', err.message);
    disp(err.getReport('extended'));
    exit(1);
end

% refreshViews logic on the bundle (table rows)
b = bundle;
rows = cell(0, 6);
for i = 1:numel(b.metrics)
    for j = 1:numel(b.metrics(i).estimators)
        e = b.metrics(i).estimators(j);
        rows(end + 1, :) = {b.metrics(i).name, e.name, ...
                            e.axis.r, e.axis.p, e.axis.n, e.axis.null}; %#ok<AGROW>
    end
end
fprintf('axis rows: %d\n', size(rows, 1));
disp(rows);

% epicentres branch
i0 = 1;
e0 = b.metrics(i0).estimators(strcmp({b.metrics(i0).estimators.name}, ...
                                     b.metrics(i0).primaryEstimator));
fprintf('primaryEstimator=%s level values=%d epicentres=%d\n', ...
    b.metrics(i0).primaryEstimator, numel(e0.values), b.metrics(i0).epicentres.n);
fprintf('APP FLOW OK\n');
