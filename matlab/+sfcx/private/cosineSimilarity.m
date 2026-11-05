function m = cosineSimilarity(a)
%COSINESIMILARITY Cosine between the structural profiles of two regions.
nrm = sqrt(sum(a .^ 2, 2));
nrm(nrm < 1e-12) = 1;
an = a ./ nrm;
m = an * an';
end
