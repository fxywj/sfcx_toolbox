function items = listCustom()
%LISTCUSTOM Registered user-defined coupling definitions.
s = registry('list');
items = struct('name', {}, 'family', {}, 'label', {}, 'description', {});
for i = 1:numel(s)
    items(end + 1) = struct('name', s(i).name, 'family', s(i).family, ...
                            'label', s(i).label, 'description', s(i).description);
end
end
