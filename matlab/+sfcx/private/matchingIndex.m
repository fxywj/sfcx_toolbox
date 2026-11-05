function m = matchingIndex(a)
%MATCHINGINDEX Shared-neighbour ratio between two regions.
a2 = a * a;
k = sum(a, 2);
den = k + k' - 2 * a;
m = zeros(size(a));
ok = den > 0;
m(ok) = a2(ok) ./ den(ok);
m(~isfinite(m)) = 0;
end
