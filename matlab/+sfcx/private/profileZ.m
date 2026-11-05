function [zi, zj] = profileZ(x, coh)
%PROFILEZ Standardise every edge inside the connection profile of each endpoint.
%   Returns the two endpoint-specific standardised versions of the same edge
%   vector, so that a coupling field can be symmetrised by averaging them.
x = double(x);
if isvector(x)
    x = reshape(x, 1, []);
end
[M, counts] = incidence(coh.nNodes);
s1 = x * M';
s2 = (x .^ 2) * M';
mu = s1 ./ counts';
va = max(s2 ./ counts' - mu .^ 2, 0);
sd = sqrt(va);
sd(sd < 1e-12) = 1;
zi = (x - mu(:, coh.iu)) ./ sd(:, coh.iu);
zj = (x - mu(:, coh.iv)) ./ sd(:, coh.iv);
end
