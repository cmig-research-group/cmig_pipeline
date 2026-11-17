function cmig_run_bids_module(module_name,varargin)
%function cmig_run_bids_module(module_name,[args])
%
% Purpose: run one cmig bids module, selected by name, so that a single
%   compiled executable can serve every module
%
% Required Input:
%   module_name: name of the module function (e.g. 'cmig_proc_bids_func')
%
% Optional Input:
%   args: passed unchanged to the module, in the module's own order
%     (positional arguments, then 'key',value pairs)
%
% Environment:
%   CMIG_MATLAB_NTHREADS: number of computational threads (maxNumCompThreads);
%     if not set, MATLAB's default (the cores available)
%
% Created:  10/06/2026 by Don Hagler
% Prev Mod: 10/06/2026 by Don Hagler
% Last Mod: 11/03/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the compiler does not follow feval of a name held in a string, so the
%   modules are listed here to be included in the build
%#function cmig_proc_bids_anat cmig_proc_bids_func
%#function cmig_analyze_bids_anat cmig_analyze_bids_func cmig_analyze_bids_dwi
%#function cmig_rawqc_bids
%#function cmig_postqc_bids_anat cmig_postqc_bids_func cmig_postqc_bids_dwi
%#function cmig_postqc_bids_fsurf cmig_postqc_bids_fmap
%#function cmig_concat_bids_sessions

module_names = {...
  'cmig_proc_bids_anat','cmig_proc_bids_func',...
  'cmig_analyze_bids_anat','cmig_analyze_bids_func','cmig_analyze_bids_dwi',...
  'cmig_rawqc_bids',...
  'cmig_postqc_bids_anat','cmig_postqc_bids_func','cmig_postqc_bids_dwi',...
  'cmig_postqc_bids_fsurf','cmig_postqc_bids_fmap',...
  'cmig_concat_bids_sessions'};

if ~mmil_check_nargs(nargin,1), return; end

if ~ismember(module_name,module_names)
  error('%s: unknown module %s; expected one of: %s',...
    mfilename,module_name,strjoin(module_names,', '));
end

% a fixed thread count, since MATLAB's results depend on its thread count,
%   which otherwise follows the host
nthreads = str2double(getenv('CMIG_MATLAB_NTHREADS'));
if ~isnan(nthreads) && nthreads > 0
  maxNumCompThreads(nthreads);
end
fprintf('%s: %d computational threads\n',mfilename,maxNumCompThreads);

feval(module_name,varargin{:});

return;
