function [stack, ids] = batchStack(path, kind)
%BATCHSTACK Stack a folder of per-subject files into one 3-D array.
%
%   [S, ids] = sfcx.batchStack(folder, 'sc')     % each file: one P x P matrix
%   [S, ids] = sfcx.batchStack(folder, 'bold')   % each file: one P x T matrix
%
%   Accepts .csv / .txt / .tsv / .mat (a .mat must contain exactly one
%   numeric array, or one whose name matches kind).  Files are sorted by
%   name; ids are the file names without extension.
%   Returns S: subjects x P x P (sc/fc) or subjects x P x T (bold).
if nargin < 2
    kind = 'sc';
end
files = [dir(fullfile(path, '*.csv')); dir(fullfile(path, '*.txt')); ...
         dir(fullfile(path, '*.tsv')); dir(fullfile(path, '*.mat'))];
if isempty(files)
    error('sfcx:batch', 'no .csv/.txt/.tsv/.mat files in %s', path);
end
[~, order] = sort(lower({files.name}));
files = files(order);
n = numel(files);
cells = cell(1, n);
ids = cell(1, n);
for k = 1:n
    full = fullfile(files(k).folder, files(k).name);
    [~, base, ext] = fileparts(files(k).name);
    ids{k} = base;
    if strcmpi(ext, '.mat')
        s = load(full);
        fn = fieldnames(s);
        v = [];
        if numel(fn) == 1
            v = s.(fn{1});
        else
            hit = find(strcmpi(fn, kind), 1);
            if isempty(hit)
                hit = find(cellfun(@(f) isnumeric(s.(f)) && ndims(s.(f)) >= 2, fn), 1);
            end
            if ~isempty(hit)
                v = s.(fn{hit});
            end
        end
        if isempty(v)
            error('sfcx:batch', '%s: no numeric array to read', files(k).name);
        end
    else
        v = readmatrix(full);
    end
    v = double(v);
    if strcmpi(kind, 'bold') && size(v, 1) > size(v, 2)
        v = v';  % time x nodes -> nodes x time
    end
    cells{k} = v;
end
ref = size(cells{1});
for k = 2:n
    if ~isequal(size(cells{k}), ref)
        error('sfcx:batch', '%s has size %s but %s has %s', ...
            ids{k}, mat2str(size(cells{k})), ids{1}, mat2str(ref));
    end
end
stack = zeros([n, ref], 'like', cells{1});
for k = 1:n
    stack(k, :, :, :) = cells{k};
end
end
