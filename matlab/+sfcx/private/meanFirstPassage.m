function m = meanFirstPassage(a)
%MEANFIRSTPASSAGE Random-walk hitting times, inverse-normalised.
n = size(a, 1);
k = sum(a, 2);
k(k < 1e-12) = 1;
P = a ./ k;
pi = k / sum(k);
Z = inv(eye(n) - P + ones(n, 1) * pi');
d = diag(Z);
m = (d' - Z) ./ max(pi', 1e-12);
m = max(m, 0);
end
