function m = diffusionEfficiency(a)
%DIFFUSIONEFFICIENCY 1 / L_ij on weighted shortest paths (0 when unreachable).
w = zeros(size(a));
ok = a > 0;
w(ok) = 1 ./ max(a(ok), 1e-12);
g = graph(w);
d = distances(g);
m = zeros(size(a));
valid = isfinite(d) & d > 0;
m(valid) = 1 ./ d(valid);
end
