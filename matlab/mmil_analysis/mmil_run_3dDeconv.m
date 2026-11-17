function [output,errcode] = mmil_run_3dDeconv(fnames_in,varargin)
%function [output,errcode] = mmil_run_3dDeconv(fnames_in,varargin)
%
% Purpose: write and execute a single AFNI 3dDeconvolve (or 3dREMLfit) shell
%   script for a single set of input files
%
% Usage:
%   mmil_run_3dDeconv(fnames_in,'key1',value1,...);
%
% Required Parameters:
%   fnames_in: string or column cell array of input NIfTI file names
%     if empty, -nodata is used instead of -input (design_flag=1 only)
%
% Optional Parameters:
%   'outstem': full path or relative output file stem
%     if empty and fnames_in is non-empty, will use file stem of first element
%     if empty and fnames_in is empty, will use 'deconv'
%     if full path, outdir will be set to directory containing outstem
%     {default = []}
%   'outdir': output directory
%     if empty, will use directory of first element of fnames_in
%     {default = []}
%   'deconvscript': full path of shell script to write and execute
%     if empty, will be set based on outstem
%     {default = []}
%   'design_flag': [0|1] only write design matrix, do not fit
%     {default = 0}
%   'nodata_nt': number of time points for -nodata (used when fnames_in is empty)
%     {default = []}
%   'nodata_tr': TR in seconds for -nodata (used when fnames_in is empty)
%     {default = []}
%   'errts_flag': [0|1] output residual time series
%     {default = 0}
%   'stim_fnames': cell array of stimulus file names, one per condition
%     {default = {}}
%   'stim_labels': cell array of stimulus condition labels
%     {default = {}}
%   'stim_times_flag': [0|1|2|3] stimulus timing format
%     0 = .1D (TR-locked), 1 = .txt (event times),
%     2 = _block.txt + _dur.txt, 3 = _block.txt with impulse
%     {default = 3}
%   'stim_times_model': HRF model name for -stim_times
%     {default = 'SPMG'}
%   'stim_times_nbasis': number of basis functions (TENT/CSPLIN only)
%     {default = 10}
%   'fname_motion': full path of motion .1D file
%     {default = []}
%   'fname_deriv': full path of motion derivative .1D file
%     {default = []}
%   'fname_censor': full path of censor .1D file
%     {default = []}
%   'motion_labels': cell array of motion regressor labels
%     {default = {'Roll' 'Pitch' 'Yaw' 'dS' 'dL' 'dP'}}
%   'glt_fnames': cell array of GLT matrix file names
%     {default = {}}
%   'glt_labels': cell array of GLT labels
%     {default = {}}
%   'detrend': polynomial detrending order (0-3)
%     {default = 3}
%   'minlag': minimum lag in TRs (for .1D stim files)
%     {default = 0}
%   'maxlag': maximum lag in TRs (for .1D stim files)
%     {default = 4}
%   'censor_flag': [0|1] apply motion censoring
%     {default = 1}
%   'deriv_flag': [0|1] include motion derivatives as regressors
%     {default = 1}
%   'reml_flag': [0|1] use 3dREMLfit for pre-whitening
%     {default = 0}
%   'baseline_flag': [0|1] add polort-0 baseline output for post-normalization
%     when ~reml_flag: adds -bout to 3dDeconvolve
%     when reml_flag: adds -Rbeta to 3dREMLfit; path returned in output.glm_Rbeta_brik
%     {default = 0}
%   'iresp_flag': [0|1] output impulse response functions
%     {default = 0}
%   'goforit_flag': [0|1] use -GOFORIT option
%     {default = 1}
%   'goforit_val': threshold value for -GOFORIT
%     {default = 10}
%   'njobs': number of 3dDeconvolve sub-processes (-jobs), splitting the
%     voxels; 1 writes no -jobs
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
%
% Output:
%   output: struct containing output file names
%   errcode: [0|1] whether there was an error running 3dDeconvolve
%
% Created:  04/09/26 by Don Hagler
% Prev Mod: 05/11/26 by Don Hagler
% Prev Mod: 06/29/26 by Don Hagler
% Prev Mod: 07/02/26 by Don Hagler
% Last Mod: 10/08/26 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on mmil_3dDeconv:run_deconv
%   Created:  09/01/08 by Don Hagler
%   Prev Mod: 01/14/19 by Dani Cornejo
%   Last Mod: 12/27/25 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = []; errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end

parms = mmil_args2parms(varargin,{ ...
  'outdir',[],[],...
  'outstem',[],[],...
  'deconvscript',[],[],...
  'design_flag',false,[false true],...
  'nodata_nt',[],[],...
  'nodata_tr',[],[],...
  'errts_flag',false,[false true],...
  'stim_fnames',{},[],...
  'stim_labels',{},[],...
  'stim_times_flag',3,0:3,...
  'stim_times_model','SPMG',{'SPMG','SPMG1','SPMG2','GAM','TENT','TENTzero','CSPLIN','CSPLINzero','BLOCK'},...
  'stim_times_nbasis',10,[1,100],...
  'fname_motion',[],[],...
  'fname_deriv',[],[],...
  'fname_censor',[],[],...
  'motion_labels',{'Roll' 'Pitch' 'Yaw' 'dS' 'dL' 'dP'},[],...
  'glt_fnames',{},[],...
  'glt_labels',{},[],...
  'detrend',3,0:3,...
  'minlag',0,[0,10],...
  'maxlag',4,[0,30],...
  'censor_flag',true,[false true],...
  'deriv_flag',true,[false true],...
  'reml_flag',false,[false true],...
  'baseline_flag',false,[false true],...
  'iresp_flag',false,[false true],...
  'goforit_flag',true,[false true],...
  'goforit_val',10,1:20,...
  'njobs',1,1:32,...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.reml_flag && ~ismember(parms.stim_times_model,{'GAM','SPMG','SPMG1','SPMG2'})
  fprintf('%s: WARNING: disabling reml_flag; stim_times_model %s not supported with REML\n',...
    mfilename,parms.stim_times_model);
  parms.reml_flag = false;
end

if isempty(fnames_in)
  fnames_in = {};
elseif ~iscell(fnames_in)
  fnames_in = {fnames_in};
end

% set outstem and outdir if not provided
if isempty(parms.outstem)
  if isempty(fnames_in)
    if isempty(parms.outdir), parms.outdir = pwd; end
    parms.outstem = sprintf('%s/deconv',parms.outdir);
  else
    [tmp_fpath,tmp_fstem,~] = fileparts(fnames_in{1});
    if isempty(parms.outdir), parms.outdir = tmp_fpath; end
    if isempty(parms.outdir), parms.outdir = pwd; end
    parms.outstem = sprintf('%s/%s',parms.outdir,tmp_fstem);
  end
elseif mmil_isrelative(parms.outstem)
  if isempty(parms.outdir), parms.outdir = pwd; end
  parms.outstem = sprintf('%s/%s',parms.outdir,parms.outstem);
end
[parms.outdir,outstem_local] = fileparts(parms.outstem);
mmil_mkdir(parms.outdir);

% set deconvscript if not provided
if isempty(parms.deconvscript)
  parms.deconvscript = sprintf('%s-deconv.sh',parms.outstem);
end

% derive counts from inputs
parms.nconds = length(parms.stim_fnames);
if isempty(parms.stim_labels)
  for c=1:parms.nconds
    parms.stim_labels{c} = sprintf('cond%d',c);
  end
end
if ~isempty(parms.fname_motion)
  parms.nmotion = length(parms.motion_labels);
else
  parms.nmotion = 0;
end
parms.nstims = parms.nconds + parms.nmotion;
if parms.deriv_flag
  parms.nstims = parms.nconds + 2*parms.nmotion;
end
parms.nglt = length(parms.glt_fnames);

% construct output file names
output = [];
output.glm_design = sprintf('%s.xmat.1D',parms.outstem);
if ~parms.design_flag
  output.glm_results_brik = sprintf('%s+orig.BRIK',parms.outstem);
  output.glm_results_head = strrep(output.glm_results_brik,'BRIK','HEAD');
  if parms.errts_flag
    output.glm_errts_brik = strrep(output.glm_results_brik,'3dDeconv','3dDeconv_errts');
    output.glm_errts_head = strrep(output.glm_results_head,'3dDeconv','3dDeconv_errts');
  end
  if parms.iresp_flag && ~parms.design_flag
    for c=1:parms.nconds
      output.glm_iresp_brik{c} = sprintf('%s_iresp_%s+orig.BRIK',...
        parms.outstem,parms.stim_labels{c});
      output.glm_iresp_head{c} = sprintf('%s_iresp_%s+orig.HEAD',...
        parms.outstem,parms.stim_labels{c});
    end
  end
  if parms.reml_flag && parms.baseline_flag
    output.glm_Rbeta_brik = sprintf('%s_Rbeta+orig.BRIK',parms.outstem);
    output.glm_Rbeta_head = strrep(output.glm_Rbeta_brik,'BRIK','HEAD');
  end
end

% check if output files already exist
all_fname_out = struct2files(output);
[exist_flag,missing_file] = check_all_exist(all_fname_out);

if ~exist(parms.deconvscript,'file') || ~exist_flag || parms.forceflag
  % remove old outputs
  if exist(parms.deconvscript,'file'), delete(parms.deconvscript); end
  all_fname_out = struct2files(output);

  for i=1:length(all_fname_out)
    if exist(all_fname_out{i},'file')
      fprintf('%s: deleting %s\n',mfilename,all_fname_out{i});
      delete(all_fname_out{i});
    end
  end

  % write shell script
  fprintf('%s: creating shell script..\n',mfilename);
  fid = fopen(parms.deconvscript,'wt');
  if fid == -1
    error('%s: unable to open file %s for writing',mfilename,parms.deconvscript);
  end
  fprintf(fid,'#!/bin/sh\n');
  fprintf(fid,'cd %s\n',parms.outdir);
  fprintf(fid,'3dDeconvolve \\\n');
  if isempty(fnames_in)
    if ~isempty(parms.nodata_nt) && ~isempty(parms.nodata_tr)
      fprintf(fid,'  -nodata %d %g \\\n',parms.nodata_nt,parms.nodata_tr);
    else
      fprintf(fid,'  -nodata \\\n');
    end
  else
    fprintf(fid,'  -input \\\n');
    for f=1:length(fnames_in)
      fprintf(fid,'    ''%s'' \\\n',fnames_in{f});
    end
    fprintf(fid,'  \\\n');
  end
  fprintf(fid,'  -polort %d \\\n',parms.detrend);
  if parms.njobs > 1 && ~parms.design_flag
    fprintf(fid,'  -jobs %d \\\n',parms.njobs);
  end
  if parms.reml_flag
    fprintf(fid,'  -x1D_stop \\\n');
  end
  if parms.censor_flag && ~parms.design_flag
    fprintf(fid,'  -censor %s \\\n',parms.fname_censor);
  end
  fprintf(fid,'  -num_stimts %d \\\n',parms.nstims);
  fprintf(fid,'  -allzero_OK \\\n');
  fprintf(fid,'  -local_times \\\n');
  s = 1;
  for c=1:parms.nconds
    if parms.stim_times_flag
      switch parms.stim_times_model
        case {'GAM'}
          switch parms.stim_times_flag
            case 1
              stim_model = parms.stim_times_model;
            case 2
              [block_filepath,block_name,~] = fileparts(parms.stim_fnames{c});
              block_dur_path = sprintf('%s/%s_dur.txt',block_filepath,block_name);
              block_dur = load(block_dur_path);
              stim_model = sprintf('%s(8.6,.547,%0.4f)',parms.stim_times_model,block_dur);
            case 3
              stim_model = sprintf('%s(8.6,.547)',parms.stim_times_model);
          end
        case {'SPMG','SPMG1','SPMG2'}
          switch parms.stim_times_flag
            case 1
              stim_model = parms.stim_times_model;
            case 2
              [block_filepath,block_name,~] = fileparts(parms.stim_fnames{c});
              block_dur_path = sprintf('%s/%s_dur.txt',block_filepath,block_name);
              block_dur = load(block_dur_path);
              stim_model = sprintf('%s(%0.4f)',parms.stim_times_model,block_dur);
            case 3
              stim_model = sprintf('%s(0)',parms.stim_times_model);
          end
        case {'TENT','TENTzero','CSPLIN','CSPLINzero'}
          stim_model = sprintf('%s(%d,%d,%d)',...
            parms.stim_times_model,parms.minlag,parms.maxlag,parms.stim_times_nbasis);
        case 'BLOCK'
          switch parms.stim_times_flag
            case 1
              stim_model = 'BLOCK(1,1)';
            case 2
              [block_filepath,block_name,~] = fileparts(parms.stim_fnames{c});
              block_dur_path = sprintf('%s/%s_dur.txt',block_filepath,block_name);
              block_dur = load(block_dur_path);
              stim_model = sprintf('BLOCK(%0.4f,1)',block_dur);
            case 3
              stim_model = 'BLOCK(0.01,1)';
          end
        otherwise
          error('%s: unsupported stim times model: %s',mfilename,parms.stim_times_model);
      end
      fprintf(fid,'  -stim_times %d ''%s'' ''%s'' -stim_label %d %s \\\n',...
        s,parms.stim_fnames{c},stim_model,s,parms.stim_labels{c});
    else
      fprintf(fid,'  -stim_file %d ''%s''  -stim_label %d %s \\\n',...
        s,parms.stim_fnames{c},s,parms.stim_labels{c});
      fprintf(fid,'  -stim_minlag %d %d  -stim_maxlag %d %d \\\n',...
        s,parms.minlag,s,parms.maxlag);
    end
    if parms.iresp_flag
      fprintf(fid,'  -iresp %d %s_iresp_%s \\\n',...
        s,outstem_local,parms.stim_labels{c});
    end
    s = s+1;
  end
  for m=1:parms.nmotion
    fprintf(fid,'  -stim_file %d ''%s[%d]'' -stim_label %d %s -stim_base %d \\\n',...
      s,parms.fname_motion,m,s,parms.motion_labels{m},s);
    s = s+1;
  end
  if parms.deriv_flag
    for m=1:parms.nmotion
      fprintf(fid,'  -stim_file %d ''%s[%d]'' -stim_label %d %s_d -stim_base %d \\\n',...
        s,parms.fname_deriv,m,s,parms.motion_labels{m},s);
      s = s+1;
    end
  end
  if parms.nglt > 0
    fprintf(fid,'  -num_glt %d \\\n',parms.nglt);
    for g=1:parms.nglt
      fprintf(fid,'  -glt 1 ''%s'' -glt_label %d "%s" \\\n',...
        parms.glt_fnames{g},g,parms.glt_labels{g});
    end
  end
  if parms.reml_flag && ~parms.design_flag
    % number of stim beta volumes to skip in REMLfit output
    switch parms.stim_times_model
      case {'GAM','SPMG1'}
        reml_stims = 2*size(parms.stim_labels,2);
      case {'SPMG','SPMG2'}
        reml_stims = 4*size(parms.stim_labels,2);
      otherwise
        error('%s: unsupported stim times model with REML: %s',...
          mfilename,parms.stim_times_model);
    end
    fprintf(fid,'  -xjpeg ''%s'' \\\n',outstem_local);
    fprintf(fid,'  -x1D ''%s.xmat.1D'' \n',outstem_local);
    fprintf(fid,'\n');
    fprintf(fid,'3dREMLfit -matrix ''%s.xmat.1D'' \\\n',outstem_local);
    if ~isempty(fnames_in)
      fprintf(fid,'  -input \\\n');
      for f=1:length(fnames_in)
        fprintf(fid,'    ''%s'' \\\n',fnames_in{f});
      end
    end
    fprintf(fid,'  -tout -noFDR \\\n');
    fprintf(fid,'  -Rbuck ''%s_stims'' \\\n',parms.outstem);
    if parms.goforit_flag
      fprintf(fid,'  -GOFORIT ''%s'' \\\n',num2str(parms.goforit_val));
    end
    if parms.errts_flag
      outstem_errts_local = strrep(outstem_local,'3dDeconv','3dDeconv_errts');
      fprintf(fid,'  -Rerrts ''%s'' \\\n',outstem_errts_local);
    end
    if parms.baseline_flag
      fprintf(fid,'  -Rbeta ''%s_Rbeta'' \\\n',parms.outstem);
    end
    fprintf(fid,'  -verb \n');
    fprintf(fid,'\n');
    fprintf(fid,'3dTcat -prefix ''%s'' \\\n',outstem_local);
    fprintf(fid,'''%s_stims+orig[0,%d..$]'' \n',outstem_local,reml_stims+1);
  else
    if parms.baseline_flag
      % -nocout dropped: -nocout suppresses ALL coefficient output (including
      % baseline), overriding -bout; must keep stim coefficient output to get
      % baseline coefficients out via -bout
      fprintf(fid,'  -tout -float \\\n');
      fprintf(fid,'  -bout \\\n');
    else
      fprintf(fid,'  -nocout -tout -float \\\n');
    end
    fprintf(fid,'  -xjpeg ''%s'' \\\n',outstem_local);
    fprintf(fid,'  -x1D ''%s.xmat.1D'' \\\n',outstem_local);
    if parms.design_flag
      fprintf(fid,'  -x1D_stop \\\n');
    end
    if parms.goforit_flag
      fprintf(fid,'  -GOFORIT ''%s'' \\\n',num2str(parms.goforit_val));
    end
    if parms.errts_flag
      outstem_errts_local = strrep(outstem_local,'3dDeconv','3dDeconv_errts');
      fprintf(fid,'  -errts ''%s'' \\\n',outstem_errts_local);
    end
    fprintf(fid,'  -bucket ''%s''\n',outstem_local);
  end
  fclose(fid);

  % run shell script
  fprintf('%s: running 3dDeconvolve with %s...\n',mfilename,parms.deconvscript);
  [status,result] = mmil_unix(sprintf('source %s',parms.deconvscript));
  disp(result);
  if status
    if ~isempty(regexp(result,'Insufficient data'))
      fprintf('%s: ERROR: insufficient data for GLM\n',mfilename);
      errcode = 1; return;
    elseif ~isempty(regexp(result,'Implies strong collinearity in the matrix columns!'))
      fprintf('%s: ERROR: collinearity error for GLM\n',mfilename);
      errcode = 1; return;
    else
      error('%s: failed to run shell script %s\n%s',mfilename,parms.deconvscript,result);
    end
  end

  % remove REML_cmd file if present
  fname_reml_cmd = sprintf('%s.REML_cmd',regexprep(parms.outstem,'\..+',''));
  if exist(fname_reml_cmd,'file'), delete(fname_reml_cmd); end

end % if script/outputs missing or forceflag

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_fname_out = struct2files(output)
  output_fields = fieldnames(output);
  all_fname_out = {};
  for i=1:length(output_fields)
    fname_out_data = output.(output_fields{i});
    if ~iscell(fname_out_data), fname_out_data = {fname_out_data}; end
    for j=1:length(fname_out_data)
      all_fname_out{end+1} = fname_out_data{j};
    end
  end
  all_fname_out = mmil_colvec(all_fname_out);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [exist_flag,missing_file] = check_all_exist(all_files)
  exist_flag = 1; missing_file = [];
  for i=1:length(all_files)
    if ~exist(all_files{i},'file')
      [~,name,ext] = fileparts(all_files{i});
      missing_file = [name,ext];
      exist_flag = 0; break;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

