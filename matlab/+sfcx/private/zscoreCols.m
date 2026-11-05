function M = zscoreCols(M)
%ZSCORECOLS Z-score each column of a matrix.
M = double(M);
for i = 1:size(M, 2)
    col = M(:, i);
    s = std(col, 1);
    if s > 1e-10
        M(:, i) = (col - mean(col)) / s;
    end
end
end
