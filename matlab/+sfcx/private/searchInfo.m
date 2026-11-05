function m = searchInfo(a)
%SEARCHINFO Cost of the shortest path under the random-walk prior.
n = size(a, 1);
w = zeros(size(a));
ok = a > 0;
w(ok) = 1 ./ max(a(ok), 1e-12);
g = graph(w);
k = sum(a, 2);
k(k < 1e-12) = 1;
logp = -log2(1 ./ k);
m = zeros(n);
for i = 1:n
    [~, pred] = shortestpathtree(g, i);
    for j = 1:n
        if i == j || isnan(pred(j))
            continue;
        end
        s = 0;
        node = j;
        guard = 0;
        while node ~= i && node > 0 && guard < n
            s = s + logp(node);
            node = pred(node);
            guard = guard + 1;
        end
        m(i, j) = s;
    end
end
end
