function [coh, rep] = loadCohort(path, varargin)
%LOADCOHORT Load a cohort from a .mat bundle or a numeric CSV.
%
%   [coh, rep] = sfcx.loadCohort('cohort.mat')
%   Recognised .mat fields: sc, fc, bold, networks, hierarchy, coords,
%   behaviour, covariates (plus a few aliases).
opts = struct('fcKind', 'pearson', 'timeAxis', 'auto', 'state', 'rest');
k = 1;
while k <= numel(varargin)
    if ~isfield(opts, varargin{k})
        error('sfcx:load', 'unknown option %s', varargin{k});
    end
    opts.(varargin{k}) = varargin{k + 1};
    k = k + 2;
end

[~, ~, ext] = fileparts(path);
if any(strcmpi(ext, {'.mat', '.mat73'}))
    s = load(path);
else
    m = readmatrix(path);
    s = struct('sc', m, 'fc', m);
end

cfg = struct();
cfg.sc = pick(s, {'sc', 'SC', 'sc_matrix', 'struct'});
cfg.fc = pick(s, {'fc', 'FC', 'fc_matrix', 'func'});
cfg.bold = pick(s, {'bold', 'BOLD', 'timeseries', 'ts'});
cfg.networks = pick(s, {'networks', 'network', 'yeo', 'communities'});
cfg.hierarchy = pick(s, {'hierarchy', 'hier', 'level'});
cfg.coords = pick(s, {'coords', 'xyz', 'sphere'});
cfg.behaviour = pick(s, {'behaviour', 'behavior', 'behav', 'y'});
cfg.covariates = pick(s, {'covariates', 'cov'});
cfg.fcKind = opts.fcKind;
cfg.timeAxis = opts.timeAxis;
cfg.state = opts.state;

if isempty(cfg.sc) && isempty(cfg.bold)
    error('sfcx:load', 'could not find sc/fc/bold fields in %s', path);
end
[coh, rep] = sfcx.buildCohort(cfg);
end

function v = pick(s, names)
v = [];
for i = 1:numel(names)
    if isfield(s, names{i}) && ~isempty(s.(names{i}))
        v = s.(names{i});
        return
    end
end
end
