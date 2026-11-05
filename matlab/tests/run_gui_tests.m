function run_gui_tests()
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
cd(fullfile(here, '..'));

fprintf('== 1. batchStack on temp folder ==\n');
tmp = fullfile(tempdir, sprintf('sfcx_batch_%d', randi(1e6)));
if ~exist(tmp, 'dir'), mkdir(tmp); end
rng(7);
for s = 1:5
    m = rand(12) * (rand(12) < 0.3);
    m = triu(m, 1); m = m + m';
    writematrix(m, fullfile(tmp, sprintf('sub%02d.csv', s)));
end
[S, ids] = sfcx.batchStack(tmp, 'sc');
assert(isequal(size(S), [5 12 12]), 'batchStack size');
assert(numel(ids) == 5, 'batchStack ids');
fprintf('   PASS  5 subjects stacked, size %s\n', mat2str(size(S)));

fprintf('== 2. buildCohort with batch 3-D ==\n');
bold = rand(5, 12, 60);
cfg = struct('sc', S, 'bold', bold);
[coh, ~] = sfcx.buildCohort(cfg);
assert(coh.n == 5 && coh.nNodes == 12, 'batch cohort shape');
fprintf('   PASS  cohort %d subjects x %d nodes\n', coh.n, coh.nNodes);

fprintf('== 3. couplingValues multi + single ==\n');
vals = sfcx.couplingValues(coh, {'direct', 'communicability'});
assert(numel(vals) == 2, 'two metrics');
assert(isequal(size(vals(1).field), [5, coh.nEdges]), 'edge field size');
assert(isequal(size(vals(1).regional), [5, 12]), 'regional size');
assert(~isempty(vals(1).meanMatrix) && isequal(size(vals(1).meanMatrix), [12 12]), 'mean matrix');
coh1 = sfcx.buildCohort(struct('sc', S(1, :, :), 'bold', reshape(bold(1, :, :), 1, 12, 60)));
vals1 = sfcx.couplingValues(coh1, {'direct'});
assert(size(vals1(1).field, 1) == 1, 'single-subject field');
fprintf('   PASS  multi (n=5) and single (n=1) coupling values\n');

fprintf('== 4. exportCouplingValues ==\n');
out = fullfile(tempdir, sprintf('sfcx_export_%d', randi(1e6)));
sfcx.exportCouplingValues(vals, out, {'csv', 'mat'});
f1 = fullfile(out, 'direct_coupling_edges.csv');
f2 = fullfile(out, 'communicability_coupling_regional.csv');
f3 = fullfile(out, 'direct_coupling.mat');
assert(exist(f1, 'file') == 2 && exist(f2, 'file') == 2 && exist(f3, 'file') == 2, 'export files');
T = readtable(f1);
assert(height(T) == 5, 'exported rows');
fprintf('   PASS  csv + mat exported to %s\n', out);

fprintf('== 5. analyze rejects single subject ==\n');
threw = false;
try
    sfcx.analyze(coh1, {'direct'});
catch err
    threw = contains(err.identifier, 'tooFewSubjects');
end
assert(threw, 'expected tooFewSubjects error');
fprintf('   PASS  single subject blocked with readable error\n');

fprintf('== 6. GUI creates headless ==\n');
figH = sfcx.app('zh');
assert(isa(figH, 'matlab.ui.Figure'), 'figure type');
kids = figH.Children;
assert(~isempty(kids), 'figure has content');
% tab titles present?
tbg = findall(figH, 'Type', 'uitabgroup');
assert(numel(tbg) >= 1, 'tab groups exist');
tabTitles = findall(figH, 'Type', 'uitab');
assert(numel(tabTitles) >= 8, 'expected >= 8 tabs (4 main + result subtabs)');
close(figH);
fprintf('   PASS  GUI built and closed\n');

cleanup = {tmp, out};
for i = 1:numel(cleanup)
    if exist(cleanup{i}, 'dir')
        rmdir(cleanup{i}, 's');
    end
end
fprintf('\nALL GUI/EXPORT TESTS PASSED\n');
end
