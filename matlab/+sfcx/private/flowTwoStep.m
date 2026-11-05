function m = flowTwoStep(a)
%FLOWTWOSTEP Two-step random-walk reachability, (D^-1 A)^2, symmetrised.
k = sum(a, 2);
k(k < 1e-12) = 1;
p = a ./ k;
m = p * p;
m = (m + m') / 2;
end
