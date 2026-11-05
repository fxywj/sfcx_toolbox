function f = edgeProjection(reg, coh)
%EDGEPROJECTION Endpoint mean of a regional map, giving an edge field.
reg = double(reg);
f = 0.5 * (reg(:, coh.iu) + reg(:, coh.iv));
end
