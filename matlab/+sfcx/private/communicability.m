function m = communicability(a)
%COMMUNICABILITY Weighted walks of every length, expm(A).
[V, D] = eig((a + a') / 2);
V = real(V);
w = max(real(diag(D)), 0);
m = (V .* exp(w)') * V';
m = (m + m') / 2;
end
