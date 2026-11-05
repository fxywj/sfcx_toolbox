function r = rankVector(x)
%RANKVECTOR Average ranks, ties shared (drop-in for tiedrank, core MATLAB only).
x = x(:);
[~, idx] = sort(x, 'ascend');
r = zeros(numel(x), 1);
r(idx) = (1:numel(x))';
[~, ~, ic] = unique(x, 'stable');
for k = 1:max(ic)
    sel = ic == k;
    if sum(sel) > 1
        r(sel) = mean(r(sel));
    end
end
end
