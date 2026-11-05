function check_example_data()
%CHECK_EXAMPLE_DATA Load the shipped 5-subject dataset and compare with the
%   reference numbers in examples/data/demo5_expected.json.
%
%   Run from anywhere:
%       addpath('...\sfcx_toolbox\matlab\tests'); check_example_data
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'matlab'));

S = load(fullfile(root, 'examples', 'data', 'demo5_group.mat'));
cfg = struct();
cfg.sc = S.sc;
cfg.fc = S.fc;
cfg.networks = S.networks;
cfg.hierarchy = S.hierarchy;
cfg.coords = S.coords;
cfg.behaviour = S.behaviour;
cfg.covariates = S.covariates;
coh = sfcx.buildCohort(cfg);
fprintf('cohort     : %d subjects x %d nodes x %d edges\n', ...
        coh.n, coh.nNodes, numel(coh.iu));

b = sfcx.analyze(coh, {'communicability', 'direct', 'reliance'}, ...
                 'estimators', {'regional', 'node', 'edge'}, ...
                 'nPerm', 500, 'seed', 0);
R = b.sweep.r;

ref = jsondecode(fileread(fullfile(root, 'examples', 'data', 'demo5_expected.json')));
names = {'communicability', 'direct', 'reliance'};
ests = {'regional', 'node', 'edge'};
fprintf('\n%-16s %-9s %9s %9s\n', 'metric', 'estimator', 'axis r', 'reference');
fprintf('%s\n', repmat('-', 1, 48));
worst = 0;
for i = 1:numel(names)
    for j = 1:numel(ests)
        r = R(i, j);
        e = ref.axis.(names{i}).(ests{j}).r;
        worst = max(worst, abs(r - e));
        fprintf('%-16s %-9s %+9.4f %+9.4f\n', names{i}, ests{j}, r, e);
    end
end
assert(worst < 1e-3, sprintf('axis r deviates from the Python reference by %.3g', worst));
fprintf('\nmax |dr| vs Python reference = %.2g  ->  PASS\n', worst);
end
