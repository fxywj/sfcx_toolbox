function W = colNormalize(A)
%COLNORMALIZE Column-normalise a weight matrix (reliance / hybrid / frontier).
W = double(A);
W = W ./ (sum(W, 1) + 1e-10);
W(1:size(W,1)+1:end) = 0;
W = W ./ (sum(W, 1) + 1e-10);
end
