function out = registry(action, varargin)
%REGISTRY Process-wide registry of user-supplied coupling definitions.
switch lower(action)
    case 'register'
        name = char(varargin{1});
        fn = varargin{2};
        meta = varargin{3};
        store = getStore();
        idx = find(strcmp({store.name}, name), 1);
        entry = struct('name', name, ...
                       'family', getField(meta, 'family', 'custom'), ...
                       'label', getField(meta, 'label', name), ...
                       'description', getField(meta, 'description', ''), ...
                       'fn', fn);
        if isempty(idx)
            store(end + 1) = entry;
        else
            store(idx) = entry;
        end
        setStore(store);
        out = name;
    case 'unregister'
        name = char(varargin{1});
        store = getStore();
        if ~isempty(store)
            store(strcmp({store.name}, name)) = [];
        end
        setStore(store);
        out = name;
    case 'list'
        out = getStore();
    case 'get'
        name = char(varargin{1});
        store = getStore();
        idx = find(strcmp({store.name}, name), 1);
        if isempty(idx)
            out = [];
        else
            out = store(idx);
        end
    case 'clear'
        setStore(struct('name', {}, 'family', {}, 'label', {}, 'description', {}, 'fn', {}));
        out = true;
    otherwise
        error('sfcx:registry', 'unknown action %s', action);
end
end

function s = getStore()
s = getappdata(0, 'sfcx_registry');
if isempty(s)
    s = struct('name', {}, 'family', {}, 'label', {}, 'description', {}, 'fn', {});
end
end

function setStore(s)
setappdata(0, 'sfcx_registry', s);
end

function v = getField(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
