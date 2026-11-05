function m = pathEase(a)
%PATHEASE Inverse shortest-path length.
w = zeros(size(a));
ok = a > 0;
w(ok) = 1 ./ max(a(ok), 1e-12);
g = graph(w);
d = distances(g);
ease = 1 ./ (1 + d);
ease(~isfinite(d)) = 0;
m = ease;
end
