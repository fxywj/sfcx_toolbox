function [r, p] = robustPearsonr(x, y)
%ROBUSTPEARSONR Pearson correlation with NaN masking.
mask = ~(isnan(x) | isnan(y));
if sum(mask) < 3
    r = 0;
    p = 1;
    return;
end
[r, p] = corr(x(mask), y(mask));
end
