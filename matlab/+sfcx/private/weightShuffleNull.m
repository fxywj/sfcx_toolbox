function aNull = weightShuffleNull(aLog, seed)
%WEIGHTSHUFFLENULL Weight-multiset permutation null model.
%   Preserves connection density and weight distribution, destroys topology.
n = size(aLog, 1);
[iu, iv] = upperTriIndex(n);
w = aLog(sub2ind([n, n], iu, iv));
wv = w(w > 0);
if isempty(wv)
    aNull = zeros(size(aLog));
    return;
end
rng = xorshift32(seed);
perm = zeros(size(wv));
for i = 1:numel(wv)
    perm(i) = rng();
end
[~, order] = sort(perm);
wvShuffled = wv(order);
tmp = zeros(size(w));
tmp(w > 0) = wvShuffled;
aNull = zeros(size(aLog));
aNull(sub2ind([n, n], iu, iv)) = tmp;
aNull = aNull + aNull';
end
