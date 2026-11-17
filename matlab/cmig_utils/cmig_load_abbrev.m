function [names, abbrevs, types] = cmig_load_abbrev(varargin)
%function [names, abbrevs, types] = cmig_load_abbrev([options])
%
% Purpose: load ROI or network abbreviation lookup table from CSV file
%
% Optional Input:
%   'indir': directory containing lookup table CSV files
%     {default = [cmig_getenv('CMIG_PIPELINE_DIR') '/parms/tabulated']}
%   'instem': input file stem (without .csv extension)
%     {default = 'roi_abbrev'}
%     other common values: 'gpnet_abbrev', 'aseg_abbrev'
%
% Output:
%   names: cell array of full names {1 x n}
%   abbrevs: cell array of abbreviations {1 x n}
%   types: cell array of type labels {1 x n}; {} if CSV has no type column
%
% Created:  05/20/2026 by Don Hagler
% Prev Mod: 05/20/2026 by Don Hagler
% Last Mod: 10/26/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {}; abbrevs = {}; types = {};

parms = mmil_args2parms(varargin,{...
  'indir',[],[],...
  'instem','roi_abbrev',[],...
});

if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

if isempty(parms.indir)
  parms.indir = [cmig_getenv('CMIG_PIPELINE_DIR') '/parms/tabulated'];
end

fname = sprintf('%s/%s.csv', parms.indir, parms.instem);
if ~exist(fname, 'file')
  fprintf('%s: WARNING: lookup table not found: %s\n', mfilename, fname);
  return;
end

s = mmil_csv2struct(fname);
names   = {s.name};
abbrevs = {s.abbrev};
if isfield(s, 'type')
  types = {s.type};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

