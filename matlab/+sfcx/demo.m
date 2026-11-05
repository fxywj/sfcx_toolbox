function coh = demo(nSubjects, nNodes, seed, withBehaviour)
%DEMO A small generative cohort with the same logic as the Python demo.
%
%   coh = sfcx.demo()                % 200 subjects, 120 regions
%   coh = sfcx.demo(80, 60, 0)
%
%   Structure is modular, function is a noisy individually weighted copy of
%   structure, and the noise load rises towards association cortex.
if nargin < 1 || isempty(nSubjects)
    nSubjects = 200;
end
if nargin < 2 || isempty(nNodes)
    nNodes = 120;
end
if nargin < 3 || isempty(seed)
    seed = 0;
end
if nargin < 4 || isempty(withBehaviour)
    withBehaviour = true;
end

stream = RandStream('mt19937ar', 'Seed', seed);
nNetworks = 6;

centres = randn(stream, nNetworks, 3);
centres = centres ./ sqrt(sum(centres .^ 2, 2));

perNet = floor(nNodes / nNetworks);
net = repelem((0:nNetworks - 1)', perNet, 1);
if numel(net) < nNodes
    net = [net; repmat(nNetworks - 1, nNodes - numel(net), 1)];
end

xyz = zeros(nNodes, 3);
for g = 0:(nNetworks - 1)
    idx = find(net == g);
    dirs = repmat(centres(g + 1, :), numel(idx), 1) + 0.55 * randn(stream, numel(idx), 3);
    dirs = dirs ./ sqrt(sum(dirs .^ 2, 2));
    xyz(idx, :) = dirs;
end

hierarchy = net / (nNetworks - 1) + 0.08 * randn(stream, nNodes, 1);
hierarchy = min(max(hierarchy, 0), 1);
hierarchyOrder = 1 + 6 * hierarchy;

d = sqrt(sum((reshape(xyz, nNodes, 1, 3) - reshape(xyz, 1, nNodes, 3)) .^ 2, 3));
same = net == net';
sc = exp(-d / 0.55) .* (1 + 0.9 * double(same));
knn = 14;
sorted = sort(sc, 2, 'descend');
thr = sorted(:, min(knn, nNodes));
sc(sc < thr) = 0;
sc = max(sc, sc');
sc(1:nNodes + 1:end) = 0;

[iu, iv] = upperTriIndex(nNodes);
scE = reshape(sc(sub2ind([nNodes nNodes], iu, iv)), 1, []);
base = scE / (max(abs(scE)) + 1e-9);

jitter = randn(stream, nNodes, 1);
adj = double(sc > 0);
deg = max(sum(adj, 2), 1);
for r = 1:2
    jitter = 0.35 * jitter + 0.65 * (adj * jitter) ./ deg;
end
jitter = (jitter - mean(jitter)) / (std(jitter, 1) + 1e-12);
nHot = max(6, round(nNodes / 12));
hot = randperm(stream, nNodes);
hot = hot(1:nHot);
for c = hot
    jitter = jitter + 1.4 * exp(-((1:nNodes)' - c) .^ 2 / (2 * 3 ^ 2));
end

noiseLoad = (0.26 + 1.35 * hierarchy) .* exp(0.7 * jitter);
sdMap = 0.05 + 0.30 * (noiseLoad / (max(noiseLoad) + 1e-12));
subjScale = exp(randn(stream, nSubjects, nNodes) .* sdMap');

sameE = double(same(sub2ind([nNodes nNodes], iu, iv)));
ampE = 0.5 * (noiseLoad(iu) + noiseLoad(iv));
ampE = ampE .* (0.15 + 0.85 * abs(base'));

fc = zeros(nSubjects, numel(scE));
for s = 1:nSubjects
    gain = 1 + 0.2 * randn(stream, 1);
    eps = randn(stream, numel(scE), 1);
    fc(s, :) = gain * base + 0.9 * 0.22 * (ampE(:) .* eps)';
end
fc = min(max(fc, -1.4), 1.4);

covariates = [randn(stream, nSubjects, 1), randn(stream, nSubjects, 1)];

% Generate synthetic BOLD time series (nSubjects, nNodes, nTime)
nTime = 200;
bold = zeros(nSubjects, nNodes, nTime);
W = sc ./ (sum(sc, 1) + 1e-10);
for s = 1:nSubjects
    ts = randn(stream, nNodes, nTime);
    % Smooth along time to create autocorrelation
    for t = 2:nTime
        ts(:, t) = 0.7 * ts(:, t - 1) + 0.3 * ts(:, t) + 0.2 * W * ts(:, t - 1);
    end
    bold(s, :, :) = ts;
end

opts = struct();
opts.coords = xyz;
opts.covariates = covariates;
opts.state = 'synthetic';
opts.bold = bold;
coh = sfcx.Cohort(repmat(scE, nSubjects, 1), fc, nNodes, net', hierarchyOrder', opts);

if withBehaviour
    c = sfcx.coupling(coh, 'communicability');
    dev = c.field - mean(c.field, 1);
    order = zeros(size(dev));
    for s = 1:nSubjects
        [~, o] = sort(dev(s, :), 'descend');
        order(s, :) = o;
    end
    ranked = zeros(size(dev));
    for s = 1:nSubjects
        ranked(s, :) = dev(s, order(s, :));
    end
    m = size(ranked, 2);
    r = (1:m)';
    r = (r - mean(r)) / std(r, 1);
    slope = (ranked * r) / (r' * r);
    z = (slope - mean(slope)) / (std(slope, 1) + 1e-12);
    coh.behaviour = 0.35 * z + randn(stream, nSubjects, 1);
end
end
