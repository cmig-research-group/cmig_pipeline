function cmig_write_tabdat_postqc(all_qcstruct, fname_tsv, forceflag)
%function cmig_write_tabdat_postqc(all_qcstruct, fname_tsv, forceflag)
%
% Purpose: write post-processing QC info struct array to a TSV file
%
% Required Input:
%   all_qcstruct: cell array of structs, one per TSV row
%   fname_tsv: output TSV file path
%   forceflag: [0|1] overwrite existing output file
%
% Created:  06/05/2026 by Don Hagler
% Prev Mod: 06/10/2026 by Don Hagler
% Prev Mod: 06/11/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end
if isempty(all_qcstruct), return; end
if exist(fname_tsv,'file') && ~forceflag, return; end

% collect all unique field names across rows (preserving insertion order)
all_fields = {};
for i = 1:length(all_qcstruct)
  flds = fieldnames(all_qcstruct{i});
  for k = 1:length(flds)
    if ~ismember(flds{k}, all_fields)
      all_fields{end+1} = flds{k};
    end
  end
end

% build data cell matrix (nrows x ncols)
nrows = length(all_qcstruct);
ncols = length(all_fields);
data = cell(nrows, ncols);
for i = 1:nrows
  row = all_qcstruct{i};
  for k = 1:ncols
    fn = all_fields{k};
    if isfield(row, fn)
      val = row.(fn);
      if isnumeric(val) && numel(val) == 1 && isnan(val)
        data{i,k} = 'n/a';
      elseif ischar(val) && strcmp(val,'_NaN_')
        data{i,k} = 'n/a';
      else
        data{i,k} = val;
      end
    else
      data{i,k} = 'n/a';
    end
  end
end

mmil_write_csv(fname_tsv, cat(1, all_fields, data), 'separator', '\t', 'quote_flag', 0);

return;
