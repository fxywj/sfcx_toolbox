function outdir = exportResults(bundle, outdir, coh)
%EXPORTRESULTS Write CSVs, a .mat bundle and a standalone HTML report.
%
%   sfcx.exportResults(bundle, 'my_run')
%   sfcx.exportResults(bundle, 'my_run', coh)
if nargin < 2 || isempty(outdir)
    outdir = fullfile(pwd, 'sfcx_out');
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
tabdir = fullfile(outdir, 'tables');
if ~exist(tabdir, 'dir')
    mkdir(tabdir);
end

names = fieldnames(bundle.tables);
for i = 1:numel(names)
    fid = fopen(fullfile(tabdir, [names{i} '.csv']), 'w');
    fprintf(fid, '%s', bundle.tables.(names{i}));
    fclose(fid);
end

payload = struct('bundle', bundle);
if nargin >= 3 && ~isempty(coh)
    payload.sc = coh.sc;
    payload.fc = coh.fc;
    payload.iu = coh.iu;
    payload.iv = coh.iv;
    payload.networks = coh.networks;
    payload.hierarchy = coh.hierarchy;
    payload.coords = coh.coords;
    payload.behaviour = coh.behaviour;
    payload.covariates = coh.covariates;
    payload.state = coh.state;
end
save(fullfile(outdir, 'sfcx_bundle.mat'), '-struct', 'payload');

writeHtml(bundle, fullfile(outdir, 'report.html'));
fprintf('results written to %s\n', outdir);
end

function writeHtml(bundle, path)
fid = fopen(path, 'w');
fprintf(fid, '<!doctype html><html><head><meta charset="utf-8">');
fprintf(fid, '<title>SFCX report</title>');
fprintf(fid, '<style>body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;');
fprintf(fid, 'max-width:960px;margin:40px auto;color:#243038;line-height:1.6}');
fprintf(fid, 'h1{font-size:22px;font-weight:600}table{border-collapse:collapse;width:100%%}');
fprintf(fid, 'td,th{border-bottom:1px solid #e6ebef;padding:6px 8px;font-size:13px;text-align:left}');
fprintf(fid, 'th{color:#6b7c88;font-weight:600}code{background:#f4f7f9;padding:1px 4px;border-radius:3px}');
fprintf(fid, '</style></head><body>');
fprintf(fid, '<h1>Structure-function coupling report</h1>');
m = bundle.meta;
fprintf(fid, '<p>%d subjects &middot; %d regions &middot; %d edges &middot; state <code>%s</code><br>', ...
    m.nSubjects, m.nNodes, m.nEdges, m.state);
fprintf(fid, 'null model <code>%s</code> with %d draws &middot; runtime %.1f s</p>', ...
    m.nullModel, m.nPerm, bundle.runtimeS);
fprintf(fid, '<h2>Variability axis across the coordinate system</h2><table>');
fprintf(fid, '<tr><th>definition</th>');
for j = 1:numel(bundle.sweep.cols)
    fprintf(fid, '<th>%s</th>', bundle.sweep.cols{j});
end
fprintf(fid, '</tr>');
for i = 1:numel(bundle.sweep.rows)
    fprintf(fid, '<tr><td>%s</td>', bundle.sweep.rows{i});
    for j = 1:numel(bundle.sweep.cols)
        r = bundle.sweep.r(i, j);
        p = bundle.sweep.p(i, j);
        if isfinite(r)
            fprintf(fid, '<td>r = %+.3f<br><span style="color:#6b7c88">p = %.3g</span></td>', r, p);
        else
            fprintf(fid, '<td>&mdash;</td>');
        end
    end
    fprintf(fid, '</tr>');
end
fprintf(fid, '</table>');
fprintf(fid, '<h2>Per-definition detail</h2><table>');
fprintf(fid, '<tr><th>definition</th><th>primary estimator</th><th>layers</th>');
fprintf(fid, '<th>split-half r</th><th>epicentres</th></tr>');
for i = 1:numel(bundle.metrics)
    b = bundle.metrics(i);
    fprintf(fid, '<tr><td>%s</td><td>%s</td>', b.name, b.primaryEstimator);
    fprintf(fid, '<td>within %.3f / between %.3f</td>', b.layers.within, b.layers.between);
    fprintf(fid, '<td>%.3f</td>', b.reliability.splitHalfR);
    fprintf(fid, '<td>%d (z = %.2f)</td>', b.epicentres.n, b.epicentres.z);
    fprintf(fid, '</tr>');
end
fprintf(fid, '</table>');
fprintf(fid, '<p style="color:#6b7c88;font-size:12px;margin-top:32px">');
fprintf(fid, 'p values come from the stated null model; reporting p without the null is not reporting p.</p>');
fprintf(fid, '</body></html>');
fclose(fid);
end
