function out = couplingValues(coh, metrics, subjectIds)
%COUPLINGVALUES Per-subject coupling magnitudes, ready for export.
%
%   out = sfcx.couplingValues(coh, {'direct','communicability'})
%
%   Works for a single subject (n = 1) and for cohorts of any size; it needs
%   no dispersion statistics, only the coupling definitions themselves.
%   out(i) has, for metrics{i}:
%       .name        metric id
%       .kind        'edge' (nSubjects x nEdges) or 'regional' (nSubjects x P)
%       .field       the coupling values themselves, one row per subject
%       .regional    nSubjects x P regional read-out (edge metrics: mean of
%                    the edges touching each region)
%       .meanVec     group-mean edge vector (nEdges x 1)  [edge metrics]
%       .meanMatrix  P x P group-mean matrix               [edge metrics]
%       .matrixNote  '' or an explanation why no matrix exists
%       .subjectIds  cell array of subject labels
if nargin < 3 || isempty(subjectIds)
    subjectIds = arrayfun(@(s) sprintf('sub%03d', s), 1:coh.n, 'UniformOutput', false);
end
if ischar(metrics) || isstring(metrics)
    metrics = cellstr(metrics);
end
if iscell(metrics) && all(cellfun(@isstring, metrics))
    metrics = cellfun(@char, metrics, 'UniformOutput', false);
end
out = struct([]);
for i = 1:numel(metrics)
    m = metrics{i};
    c = sfcx.coupling(coh, m);
    o = struct('name', m, 'kind', c.kind, 'field', double(c.field), ...
               'subjectIds', {subjectIds});
    if strcmp(c.kind, 'edge')
        o.regional = sfcx.regionalReadout(coh, m, c.field);
        o.meanVec = nanMeanVec(double(c.field), 1)';
        [iu, iv] = upperTriIndex(coh.nNodes);
        M = zeros(coh.nNodes);
        M(sub2ind([coh.nNodes coh.nNodes], iu, iv)) = o.meanVec;
        M = M + M';
        o.meanMatrix = M;
        o.matrixNote = '';
    else
        o.regional = double(c.field);
        o.meanVec = nanMeanVec(double(c.field), 1)';
        o.meanMatrix = [];
        o.matrixNote = sprintf(['%s produces a regional (not edge-wise) coupling ', ...
            'value, so there is no P x P matrix; see .regional'], m);
    end
    if isempty(out)
        out = o;
    else
        out(i) = o; %#ok<AGROW>
    end
end
end
