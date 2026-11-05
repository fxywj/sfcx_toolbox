function v = matToVec(m, iu, iv)
%MATTOVEC Read the upper-triangle edge vector out of a square matrix.
v = m(sub2ind(size(m), iu, iv));
v = v(:);
end
