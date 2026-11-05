function res = epicentres(v, coh, k, seed, nPerm)
%EPICENTRES Focal maxima with supporting neighbourhoods, plus separation z.
if nargin < 3 || isempty(k)
    k = 2;
end
if nargin < 4 || isempty(seed)
    seed = 0;
end
if nargin < 5 || isempty(nPerm)
    nPerm = 1000;
end
v = double(v(:));
a = coh.adjacency() > 0;
n = coh.nNodes;
med = median(v(isfinite(v)));

focal = [];
for i = 1:n
    nb = find(a(i, :));
    if isempty(nb)
        continue;
    end
    if v(i) >= max(v(nb), [], 'omitnan') && mean(v(nb), 'omitnan') > med
        focal(end+1) = i; %#ok<AGROW>
    end
end

sep = 0;
for x = 1:numel(focal)
    for y = (x + 1):numel(focal)
        if a(focal(x), focal(y))
            sep = sep + 1;
        end
    end
end

stream = RandStream('mt19937ar', 'Seed', seed);
null = zeros(nPerm, 1);
for i = 1:nPerm
    idx = randperm(stream, n);
    idx = idx(1:max(numel(focal), 1));
    c = 0;
    for x = 1:numel(idx)
        for y = (x + 1):numel(idx)
            if a(idx(x), idx(y))
                c = c + 1;
            end
        end
    end
    null(i) = c;
end
z = (sep - mean(null)) / (std(null, 1) + 1e-12);

prof = v;
prof(~isfinite(prof)) = med;
prof = (prof - mean(prof)) / (std(prof, 1) + 1e-12);
labels = nan(1, n);
if exist('linkage', 'file') == 2
    tree = linkage(prof, 'ward');
    labels = cluster(tree, 'maxclust', k)';
end

netCounts = zeros(1, coh.nNetworks);
for i = focal
    g = coh.networks(i);
    netCounts(g + 1) = netCounts(g + 1) + 1;
end

res = struct();
res.regions = focal;
res.n = numel(focal);
res.internalEdges = sep;
res.z = z;
res.p = (1 + sum(null <= sep)) / (1 + nPerm);
res.nullMean = mean(null);
res.clusterLabels = labels;
res.networkCounts = netCounts;
end
