function r = spearmanR(a, b)
%SPEARMANR Rank correlation with pairwise-complete observations.
a = a(:);
b = b(:);
ok = isfinite(a) & isfinite(b);
if sum(ok) < 3
    r = NaN;
    return
end
ra = rankVector(a(ok));
rb = rankVector(b(ok));
sa = std(ra, 1);
sb = std(rb, 1);
if sa < 1e-12 || sb < 1e-12
    r = NaN;
    return
end
r = mean((ra - mean(ra)) .* (rb - mean(rb))) / (sa * sb);
end
