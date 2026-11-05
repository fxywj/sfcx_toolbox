%RUN_TESTS Headless verification of the MATLAB toolbox.
%   matlab -batch "run('...\matlab\tests\run_tests.m')"
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
fprintf('toolbox root: %s\n', root);

checks = {};
checks{end + 1} = @checkCohort;
checks{end + 1} = @checkCatalog;
checks{end + 1} = @checkCoupling;
checks{end + 1} = @checkAxis;
checks{end + 1} = @checkLayers;
checks{end + 1} = @checkAnalyze;
checks{end + 1} = @checkCustom;
checks{end + 1} = @checkBold;
checks{end + 1} = @checkReference;

passed = 0;
for i = 1:numel(checks)
    name = func2str(checks{i});
    try
        checks{i}();
        fprintf('PASS  %s\n', name);
        passed = passed + 1;
    catch err
        fprintf('FAIL  %s  ->  %s\n', name, err.message);
    end
end
fprintf('\n%d/%d passed\n', passed, numel(checks));
if passed < numel(checks)
    exit(1);
end

function checkCohort()
coh = sfcx.demo(60, 40, 0);
assert(coh.n == 60 && coh.nNodes == 40 && coh.nEdges == 780);
rc = coh.regionalCoupling();
assert(isequal(size(rc), [60, 40]));
assert(isfinite(nanmean(rc(:))));
[a, b] = coh.splitHalf(0);
assert(numel(a) == 30 && numel(b) == 30 && isempty(intersect(a, b)));
end

function checkCatalog()
c = sfcx.catalog();
assert(numel(c.metrics) >= 27);
assert(any(strcmp({c.metrics.id}, 'communicability')));
assert(any(strcmp({c.metrics.id}, 'diffusion_efficiency')));
assert(any(strcmp({c.metrics.id}, 'rnc')));
assert(any(strcmp({c.metrics.id}, 'trr')));
assert(any(strcmp({c.metrics.id}, 'hybrid_global')));
assert(any(strcmp({c.estimators.id}, 'regional')));
assert(~isempty(sfcx.describe()));
end

function checkCoupling()
coh = sfcx.demo(30, 30, 0, false);
edgeMetrics = {'direct', 'rank_direct', 'cosine', 'matching', 'communicability', ...
               'degree_product', 'mfpt', 'path_ease', 'diffusion_efficiency', ...
               'search_info', 'flow', 'reliance'};
for m = edgeMetrics
    out = sfcx.coupling(coh, m{1});
    assert(strcmp(out.kind, 'edge'));
    assert(isequal(size(out.field), [30, 435]));
    assert(any(isfinite(out.field(:))), 'metric %s has no finite values', m{1});
end
regionalMetrics = {'reliance_raw', 'reliance_deriv', 'reliance_raw_nonlinear', ...
                   'reliance_deriv_nonlinear', 'information_geometry', ...
                   'multi_scale_transfer_entropy', 'criticality', 'granger', ...
                   'rnc', 'trr'};
for m = regionalMetrics
    out = sfcx.coupling(coh, m{1});
    assert(strcmp(out.kind, 'regional'));
    assert(isequal(size(out.field), [30, 30]));
    assert(any(isfinite(out.field(:))), 'metric %s has no finite values', m{1});
end
% Hybrid models are slow; test one to verify the plumbing
out = sfcx.coupling(coh, 'hybrid_ridge');
assert(strcmp(out.kind, 'regional'));
assert(isequal(size(out.field), [30, 30]));
assert(any(isfinite(out.field(:))));
end

function checkAxis()
coh = sfcx.demo(60, 40, 0, false);
c = sfcx.coupling(coh, 'direct');
v = sfcx.variability(c.field, coh, 'edge');
h = 0.5 * (coh.hierarchy(coh.iu) + coh.hierarchy(coh.iv));
res = sfcx.axis(v, h, []);
assert(isfinite(res.r));

vr = sfcx.regionalVariability(coh.regionalCoupling(), []);
engine = sfcx.nullEngine(coh, 'msr', 60, 0);
res = sfcx.axis(vr, coh.hierarchy, engine);
assert(isfinite(res.r) && res.p >= 0 && res.p <= 1);
assert(strcmp(res.null, 'msr') && res.nPerm == 60);

engine = sfcx.nullEngine(coh, 'spin', 40, 0);
res = sfcx.axis(vr, coh.hierarchy, engine);
assert(isfinite(res.p) && strcmp(res.null, 'spin'));
end

function checkLayers()
coh = sfcx.demo(80, 60, 0, false);
c = sfcx.coupling(coh, 'direct');
v = sfcx.variability(c.field, coh, 'edge');
lay = sfcx.layers(v, coh, 100, 0);
assert(isfinite(lay.within) && isfinite(lay.between));
assert(lay.p >= 0 && lay.p <= 1);
end

function checkAnalyze()
coh = sfcx.demo(60, 40, 0);
b = sfcx.analyze(coh, {'communicability', 'direct'}, ...
                 'estimators', {'regional', 'node', 'edge'}, ...
                 'nPerm', 80, 'seed', 0);
assert(numel(b.metrics) == 2);
assert(isequal(size(b.sweep.r), [2, 3]));
assert(isfinite(b.metrics(1).estimators(1).axis.r));
assert(isfinite(b.metrics(1).layers.delta));
assert(b.metrics(1).reliability.splitHalfR > -1);
assert(isfield(b.metrics(1).individual, 'edge'));
assert(isfield(b.tables, 'axis_r'));
assert(b.runtimeS >= 0);
fprintf('      sweep r = [%s]\n', num2str(b.sweep.r, '%.3f  '));
end

function checkCustom()
src = sprintf(['function out = couple(sc, fc, iu, iv, nNodes)\n' ...
               '    out = sc .* fc;\n' ...
               'end\n']);
sfcx.registerCustom('my_coupling', src, 'elementwise product');
items = sfcx.listCustom();
assert(any(strcmp({items.name}, 'my_coupling')));
coh = sfcx.demo(20, 24, 0, false);
out = sfcx.coupling(coh, 'my_coupling');
assert(isequal(size(out.field), [20, 276]));
sfcx.unregisterCustom('my_coupling');
assert(isempty(sfcx.listCustom()) || ~any(strcmp({sfcx.listCustom().name}, 'my_coupling')));
end

function checkBold()
stream = RandStream('mt19937ar', 'Seed', 1);
bold = randn(stream, 10, 20, 50);
[fc, p] = sfcx.boldToFc(bold, 'pearson');
assert(p == 20 && isequal(size(fc), [10, 190]));
fc2 = sfcx.boldToFc(bold, 'precision');
assert(isequal(size(fc2), [10, 190]));
assert(all(abs(diag(corrcoef(fc)) - 1) < 1e-6));
coh = sfcx.buildCohort(struct('sc', ones(10, 190), 'bold', bold));
assert(coh.nNodes == 20);
end

function checkReference()
here = fileparts(mfilename('fullpath'));
file = fullfile(here, 'reference.mat');
if ~exist(file, 'file')
    fprintf('      reference.mat absent, parity check skipped\n');
    return
end
s = load(file);
% Metric list, strict mask and tolerances all come from the reference file
% itself (written by make_reference.py), so the two sides cannot drift apart.
if isfield(s, 'metric_names')
    if iscell(s.metric_names)
        metrics = s.metric_names(:);
    else
        metrics = cellstr(squeeze(s.metric_names));   % char array -> cell
    end
else
    metrics = {'direct', 'rank_direct', 'cosine', 'matching', 'communicability', ...
               'degree_product', 'mfpt', 'path_ease', 'diffusion_efficiency', ...
               'search_info', 'flow', 'reliance'};
end
coh = sfcx.Cohort(s.sc, s.fc, [], s.networks, s.hierarchy, ...
                  struct('coords', s.coords, 'covariates', s.covariates));
b = sfcx.analyze(coh, metrics, ...
                 'estimators', {'regional', 'node', 'edge'}, ...
                 'nPerm', 50, 'seed', 0, 'doSweep', false);
assert(isequal(size(b.sweep.r), size(s.expected_r)), ...
       sprintf('sweep %dx%d vs reference %dx%d', ...
               size(b.sweep.r, 1), size(b.sweep.r, 2), ...
               size(s.expected_r, 1), size(s.expected_r, 2)));
d = abs(b.sweep.r - s.expected_r);
if isfield(s, 'strict')
    strict = logical(s.strict);        % false = not comparable across runtimes
    tol = s.tol;
else
    strict = ones(size(d), 'logical');
    tol = 1e-5 * ones(size(d));
end
fail = strict & (d > tol);             % only strict cells can fail
bad = find(fail);
assert(isempty(bad), sprintf('parity failed in %d cells, worst |dr| = %.3g', ...
        numel(bad), max([d(fail); 0])));
n_checked = sum(strict(:));
fprintf('      parity: %d/%d strict cells within tolerance, worst |dr| = %.2g\n', ...
        n_checked - numel(bad), n_checked, max(d(:)));
end
