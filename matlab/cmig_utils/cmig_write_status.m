function cmig_write_status(participant_id,session_id,outdir,module,errcode,errmsg,extra,startdate,procflag)
%function cmig_write_status(participant_id,session_id,outdir,module,errcode,errmsg,[extra],[startdate],[procflag])
%
% Required Input:
%   participant_id: participant identifier string
%   session_id: session identifier string
%   outdir: output directory for the status TSV file
%   module: module identifier string, carrying the desc- entity, since it
%          becomes the {participant_id}_{session_id}_{module}_status.tsv
%          filename that cmig_mkjobs_bids_pipeline.py globs for
%          e.g. 'desc-rawqc', 'desc-proc', 'desc-fsurf', 'desc-postqc',
%               'desc-analyze', sprintf('task-%s_desc-analyze',task)
%   errcode: error code (0=success, 1=error)
%   errmsg: error message string (empty if no error)
%     an empty errmsg or extra value is written as n/a (BIDS), and printed
%     as given
%
% Optional Input:
%   extra: struct whose fields are written as additional columns
%     after errmsg; field values may be numeric, logical, or char
%     {default = []}
%   startdate: the module's entry time, written as a column after errmsg and
%     printed above outdir, so a duration needs no differencing against the
%     module before it
%     {default = []}
%   procflag: [0|1] whether the module did its work. 0 with the status file
%     already present means it skipped, so the file is left as it is: a
%     rewrite would re-stamp procdate for an invocation that computed nothing
%     and discard the version columns the python wrapper appended
%     {default = 1}
%
% Created:  05/19/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 06/03/2026 by Don Hagler
% Prev Mod: 06/05/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~exist('extra','var'), extra = []; end
  if ~exist('startdate','var'), startdate = []; end
  if ~exist('procflag','var') || isempty(procflag), procflag = 1; end
  fname_status = sprintf('%s/%s_%s_%s_status.tsv',...
    outdir,participant_id,session_id,module);
  % the file-exists half matters: outputs can be present with no status file,
  %   and writing nothing there would leave mkjobs never seeing done
  if ~procflag && exist(fname_status,'file')
    [~,fstem,fext] = fileparts(fname_status);
    fprintf('%s: %s skipped, leaving %s\n',mfilename,module,[fstem fext]);
    return;
  end
  if ~ischar(errmsg), errmsg = ''; end
  errmsg = strrep(errmsg, sprintf('\n'), ' ');
  errmsg = strrep(errmsg, sprintf('\r'), ' ');
  procdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');
  fprintf('\n');
  % printed first so the block names itself; anything reading the .out would
  %   otherwise have to reconstruct the module from surrounding context, and
  %   the offsets of the four lines below are relied on downstream
  fprintf('%s: module = %s\n',mfilename,module);
  if ~isempty(startdate)
    fprintf('%s: startdate = %s\n',mfilename,startdate);
  end
  fprintf('%s: outdir = %s\n',mfilename,outdir);
  fprintf('%s: procdate = %s\n',mfilename,procdate);
  fprintf('%s: errcode = %d\n',mfilename,errcode);
  fprintf('%s: errmsg = %s\n',mfilename,errmsg);
  % convert extra columns to strings
  extra_names = {}; extra_vals = {};
  if ~isempty(extra) && isstruct(extra)
    extra_names = fieldnames(extra);
    extra_vals = cell(size(extra_names));
    for i=1:length(extra_names)
      val = extra.(extra_names{i});
      if isnumeric(val) || islogical(val)
        extra_vals{i} = num2str(val);
      elseif ischar(val)
        extra_vals{i} = val;
      else
        extra_vals{i} = '';
      end
      if isempty(extra_vals{i}), extra_vals{i} = 'n/a'; end
    end
  end
  mmil_mkdir(outdir);
  fid = fopen(fname_status,'wt');
  if fid < 0
    error('%s: failed to open %s for writing',mfilename,fname_status);
  end
  fprintf(fid,'participant_id\tsession_id\tprocdate\terrcode\terrmsg');
  if ~isempty(startdate), fprintf(fid,'\tstartdate'); end
  for i=1:length(extra_names)
    fprintf(fid,'\t%s',extra_names{i});
  end
  fprintf(fid,'\n');
  if isempty(errmsg)
    errmsg_tsv = 'n/a';
  else
    errmsg_tsv = errmsg;
  end
  fprintf(fid,'%s\t%s\t%s\t%d\t%s',participant_id,session_id,procdate,errcode,errmsg_tsv);
  if ~isempty(startdate), fprintf(fid,'\t%s',startdate); end
  for i=1:length(extra_vals)
    fprintf(fid,'\t%s',extra_vals{i});
  end
  fprintf(fid,'\n');
  fclose(fid);
return;
