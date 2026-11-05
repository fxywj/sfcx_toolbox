function out = zscoreByGroup(x, groups, nGroups)
%ZSCOREBYGROUP Standardise each row of x within each group of columns.
%   groups are 1-based labels in 1..nGroups (MATLAB-native; the Python core
%   uses 0-based labels, the arithmetic is identical).
out = zeros(size(x), 'like', double(x));
for g = 1:nGroups
    idx = find(groups == g);
    if isempty(idx)
        continue;
    end
    if numel(idx) < 2
        out(:, idx) = 0;
        continue;
    end
    block = double(x(:, idx));
    mu = mean(block, 2);
    sd = std(block, 1, 2);
    sd(sd < 1e-12) = 1;
    out(:, idx) = (block - mu) ./ sd;
end
end
