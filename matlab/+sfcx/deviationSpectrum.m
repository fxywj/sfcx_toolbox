function out = deviationSpectrum(edgeField, coh, scales, maxEdges, seed)
%DEVIATIONSPECTRUM Per-subject deviation-spectrum slope at several scales.
if nargin < 3 || isempty(scales)
    scales = {'regional', 'network', 'edge', 'gap'};
end
if ischar(scales)
    scales = {scales};
end
if nargin < 4 || isempty(maxEdges)
    maxEdges = 20000;
end
if nargin < 5 || isempty(seed)
    seed = 0;
end
ef = double(edgeField);
out = struct();

if any(strcmp(scales, 'regional'))
    out.regional = profileSlope(coh.regionalCoupling());
end
if any(strcmp(scales, 'network'))
    reg = coh.regionalCoupling();
    net = nan(size(reg, 1), coh.nNetworks());
    for g = 0:(coh.nNetworks() - 1)
        m = coh.networks == g;
        if any(m)
            net(:, g + 1) = nanMeanVec(reg(:, m), 2);
        end
    end
    out.network = profileSlope(net);
end
if any(strcmp(scales, 'edge'))
    e = ef;
    if size(e, 2) > maxEdges
        stream = RandStream('mt19937ar', 'Seed', seed);
        idx = randperm(stream, size(e, 2));
        e = e(:, idx(1:maxEdges));
    end
    out.edge = profileSlope(e);
end
if any(strcmp(scales, 'gap'))
    net = coh.networks(:);
    same = net(coh.iu) == net(coh.iv);
    w = nanMeanVec(ef(:, same), 2);
    b = nanMeanVec(ef(:, ~same), 2);
    out.gap = w - b;
end
end
