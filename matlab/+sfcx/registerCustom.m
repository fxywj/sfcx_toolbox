function name = registerCustom(name, fnOrSource, description)
%REGISTERCUSTOM Make a user-defined coupling definition a first-class metric.
%
%   sfcx.registerCustom('my_coupling', @myCoupling)
%   sfcx.registerCustom('my_coupling', sourceText)   % writes a temp .m file
%
%   Contract, both signatures accepted:
%       out = couple(sc, fc, iu, iv, nNodes)
%       out = couple(coh)
%   Return (subjects x edges), (1 x edges), (subjects x nodes), (1 x nodes),
%   (nodes x nodes) or a scalar.
name = strtrim(char(name));
if isempty(name) || ~all(isletter(name) | name == '_' | (name >= '0' & name <= '9'))
    error('sfcx:custom', 'name may contain letters, digits and underscores only');
end
builtin = {'direct', 'reliance', 'cosine', 'matching', 'communicability', ...
           'mfpt', 'path_ease', 'search_info'};
if any(strcmp(builtin, name))
    error('sfcx:custom', '''%s'' is a built-in definition; pick another name', name);
end
if nargin < 3 || isempty(description)
    description = '';
end

if ischar(fnOrSource) || isstring(fnOrSource)
    fn = compileSource(char(fnOrSource), name);
elseif isa(fnOrSource, 'function_handle')
    fn = fnOrSource;
else
    error('sfcx:custom', 'second argument must be a function handle or source text');
end

coh = sfcx.demo(8, 24, 0, false);
[kind, field] = evaluateCustom(fn, coh);
if ~any(isfinite(field(:)))
    error('sfcx:custom', 'custom function returned no finite values');
end

registry('register', name, fn, struct('family', 'custom', 'label', name, ...
                                      'description', description));
end

function fn = compileSource(source, name)
source = strtrim(source);
if isempty(source)
    error('sfcx:custom', 'empty source');
end
tok = regexp(source, '^\s*function\s+[^\n(]*?(\w+)\s*\(', 'tokens', 'once');
if isempty(tok)
    error('sfcx:custom', 'source must start with a function declaration');
end
source = regexprep(source, ...
    '^\s*function\s+([^\n(]*?)(\w+)\s*\(', ['function $1' name '('], 'once');

dirpath = fullfile(tempdir, 'sfcx_custom');
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
file = fullfile(dirpath, [name '.m']);
fid = fopen(file, 'w');
if fid < 0
    error('sfcx:custom', 'cannot write %s', file);
end
fprintf(fid, '%s\n', source);
fclose(fid);
addpath(dirpath);
% rehash rescans the whole MATLAB path, which flags unrelated user toolboxes
% that shadow builtins; keep that noise away from the app's console.
w = warning('off', 'MATLAB:dispatcher:nameConflict');
rehash path;
warning(w);
clear(name);
fn = str2func(name);
end
