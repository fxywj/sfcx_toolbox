function txt = describe()
%DESCRIBE Human-readable description of the four-layer coordinate system.
c = sfcx.catalog();
lines = {'SFC readout coordinate system', ''};
lines{end+1} = 'Layer 1 - definition';
for i = 1:numel(c.metrics)
    m = c.metrics(i);
    lines{end+1} = sprintf('  %-16s %-24s %s', m.id, m.family, m.description);
end
lines{end+1} = '';
lines{end+1} = 'Layer 2 - functional pipeline';
for i = 1:numel(c.pipelines)
    lines{end+1} = sprintf('  %-16s %s', c.pipelines(i).id, c.pipelines(i).description);
end
lines{end+1} = '';
lines{end+1} = 'Layer 3 - dispersion estimator';
for i = 1:numel(c.estimators)
    lines{end+1} = sprintf('  %-16s %s', c.estimators(i).id, c.estimators(i).description);
end
lines{end+1} = '';
lines{end+1} = 'Layer 4 - observation scale';
for i = 1:numel(c.scales)
    lines{end+1} = sprintf('  %-16s %s', c.scales(i).id, c.scales(i).description);
end
txt = strjoin(lines, newline);
end
