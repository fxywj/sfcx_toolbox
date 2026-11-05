function reg = regionalReadout(coh, metric, field)
%REGIONALREADOUT Regional summary of a coupling field.
if strcmpi(metric, 'direct')
    reg = coh.regionalCoupling();
    return
end
arr = double(field);
if isequal(size(arr), [coh.n, coh.nNodes])
    reg = arr;
    return
end
reg = nan(size(arr, 1), coh.nNodes);
for i = 1:coh.nNodes
    idx = find((coh.iu == i) | (coh.iv == i));
    if isempty(idx)
        continue;
    end
    reg(:, i) = nanMeanVec(arr(:, idx), 2);
end
end
