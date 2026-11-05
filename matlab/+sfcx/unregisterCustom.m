function name = unregisterCustom(name)
%UNREGISTERCUSTOM Remove a user-defined coupling definition.
name = char(name);
registry('unregister', name);
end
