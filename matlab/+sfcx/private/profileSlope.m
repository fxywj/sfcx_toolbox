function s = profileSlope(x)
%PROFILESLOPE Slope of each row's deviation profile against its rank.
x = double(x);
ok = isfinite(x);
c = sum(ok, 2);
meanv = sum(x .* ok, 2) ./ max(c, 1);
dev = x - meanv;
dev(~ok) = 0;
ranked = sort(dev, 2, 'descend');
m = size(ranked, 2);
if m < 2
    s = zeros(size(ranked, 1), 1);
    return
end
r = (1:m)';
r = (r - mean(r)) / std(r, 1);
s = (ranked * r) / (r' * r);
end
