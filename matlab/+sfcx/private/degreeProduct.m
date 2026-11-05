function m = degreeProduct(a)
%DEGREEPRODUCT k_i * k_j on the binarised structural matrix.
k = sum(a, 2);
m = k * k';
end
