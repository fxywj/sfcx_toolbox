function out = relianceField(coh)
%RELIANCEFIELD Per-subject slope of FC on SC, re-expressed at edge level.
[n, e] = size(coh.sc);
sc = double(coh.sc);
fc = double(coh.fc);
out = zeros(n, e);
for i = 1:coh.nNodes
    idx = find(coh.iu == i);
    if numel(idx) < 3
        continue;
    end
    s = sc(:, idx);
    f = fc(:, idx);
    sm = s - mean(s, 2);
    fm = f - mean(f, 2);
    v = sum(sm .^ 2, 2);
    v(v < 1e-12) = Inf;
    b = sum(sm .* fm, 2) ./ v;
    out(:, idx) = sm .* (b .* sm);
end
out = zscoreByGroup(out, coh.iu, coh.nNodes);
end
