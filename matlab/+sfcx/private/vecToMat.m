function m = vecToMat(v, iu, iv, nNodes)
%VECTOMAT Scatter an edge vector into a symmetric node matrix.
m = zeros(nNodes, nNodes);
v = v(:);
m(sub2ind([nNodes nNodes], iu, iv)) = v;
m(sub2ind([nNodes nNodes], iv, iu)) = v;
end
