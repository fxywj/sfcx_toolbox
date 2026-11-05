function out = exportCouplingValues(vals, dirpath, fmt)
%EXPORTCOUPLINGVALUES Write coupling magnitudes to CSV and/or MAT files.
%
%   sfcx.exportCouplingValues(vals, dir)         % both csv and mat
%   sfcx.exportCouplingValues(vals, dir, 'csv')  % csv only
%
%   For every metric in vals this writes
%       <metric>_coupling_edges.csv      nSubjects x (subjectId + nEdges)
%       <metric>_coupling_regional.csv   nSubjects x (subjectId + P)
%       <metric>_mean_matrix.mat         struct with meanMatrix (P x P)
%   and a README.txt describing every column.
if nargin < 3 || isempty(fmt)
    fmt = {'csv', 'mat'};
end
if ischar(fmt)
    fmt = {fmt};
end
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
lines = {'SFCX coupling-value export', ...
         sprintf('generated %s', datestr(now)), ...
         'Each <metric>_coupling_edges.csv: column 1 = subject id, columns 2..end = upper-triangle edge coupling (same ordering as sfcx.upperTriIndex).', ...
         'Each <metric>_coupling_regional.csv: column 1 = subject id, columns 2..end = regional coupling per parcel.', ...
         'Each <metric>_mean_matrix.mat: struct meanMatrix (P x P), group mean of the per-subject matrices.'};
for i = 1:numel(vals)
    v = vals(i);
    sid = v.subjectIds(:);
    if any(strcmp(fmt, 'csv'))
        if strcmp(v.kind, 'edge')
            T = [table(sid), array2table(double(v.field))];
            writetable(T, fullfile(dirpath, sprintf('%s_coupling_edges.csv', v.name)));
        end
        T2 = [table(sid), array2table(v.regional)];
        writetable(T2, fullfile(dirpath, sprintf('%s_coupling_regional.csv', v.name)));
    end
    if any(strcmp(fmt, 'mat'))
        S = struct('metric', v.name, 'kind', v.kind, ...
                   'field', double(v.field), 'regional', v.regional, ...
                   'subjectIds', {v.subjectIds});
        if ~isempty(v.meanMatrix)
            S.meanMatrix = v.meanMatrix;
        end
        save(fullfile(dirpath, sprintf('%s_coupling.mat', v.name)), '-struct', 'S');
    end
    lines{end + 1} = sprintf('%s: kind=%s, %d subjects', v.name, v.kind, size(v.field, 1)); %#ok<AGROW>
end
fid = fopen(fullfile(dirpath, 'README.txt'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\r\n', lines{:});
fclose(fid);
out = dirpath;
end
