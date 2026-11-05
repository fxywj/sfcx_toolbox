function p = twoSidedP(obs, null)
%TWOSIDEDP Empirical two-sided p from a surrogate distribution.
null = null(isfinite(null));
if isempty(null)
    p = NaN;
    return;
end
ge = sum(null >= obs);
le = sum(null <= obs);
p = (1 + min(ge, le)) / (1 + numel(null));
end
