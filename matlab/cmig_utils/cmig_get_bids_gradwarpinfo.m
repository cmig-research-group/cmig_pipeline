function [gwtype,unwarpflag,isoctrflag] = cmig_get_bids_gradwarpinfo(scan_info)
%function [gwtype,unwarpflag,isoctrflag] = cmig_get_bids_gradwarpinfo(scan_info)
%
% Required Input:
%   scan_info: struct containing scan information
%     as read from BIDS sidecar json file
%
% Output:
%   gwtype: gradient type number
%      0:  Siemens Sonata
%      1:  Siemens Allegra
%      2:  GE BRM
%      3:  GE CRM
%      4:  Siemens Avanto
%      5:  Siemens AXXess/Espree
%      6:  Siemens Quantum/Symphony
%      7:  GE Twin Speed Whole Body
%      8:  GE Twin Speed Zoom
%      9:  GE MR450 or MR750
%      10: GE MR750W
%      11: Siemens Skyra
%      12: Siemens Connectome Skyra
%      13: Siemens Prisma
%      14: GE Signa Premier
%      15: GE Signa UHP
%      16: Siemens Verio/Biograph
%      NaN: unknown or ambiguous
%   unwarpflag: unwarp flag
%      0: unwarp 3D
%      1: unwarp through plan only
%      2: unwarp inplane only
%      NaN: unwarp not required, or NonlinearGradientCorrection not in sidecar
%   isoctrflag: isocenter flag
%      0: scanning at isocenter
%      1: scanning off-isocenter
%
% Created:  08/06/2025 by Don Hagler
% Prev Mod: 10/06/2025 by Don Hagler
% Last Mod: 08/19/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on mmil_get_gradwarpinfo
% Created:  03/19/19 by Feng Xue
% Last Mod: 12/09/24 by Don Hagler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gwtype = NaN; unwarpflag = 0; isoctrflag = 0;

gwtype = get_gwtype(scan_info);

unwarpflag = get_unwarpflag(scan_info);

isoctrflag = get_isoctrflag(scan_info);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gwtype = get_gwtype(scan_info)
  gwtype = NaN;
  manu = lower(scan_info.Manufacturer);
  model = deblank(lower(scan_info.ManufacturersModelName));
  fieldstr = scan_info.MagneticFieldStrength;
  if ~isempty(regexp(manu,'siemens'))
    if ismember(model,{'sonata','trio','sonatavision'})
      gwtype = 0;
    elseif ismember(model,{'allegra'})
      gwtype = 1;
    elseif ismember(model,{'avanto','triotim'})
      gwtype = 4;
    elseif ismember(model,{'espree','axxess'})
      gwtype = 5;
    elseif ismember(model,{'symphony','symphonyvision','symphonytim'}) % Quantum
      gwtype = 6;
    elseif ismember(model,{'skyra'}) % Skyra
      gwtype = 11;
    elseif ismember(model,{'connectome'}) % HCP Connectome Skyra
      gwtype = 12;
    elseif ~isempty(regexp(model,'prisma'))
      gwtype = 13;
    elseif ismember(model,{'verio','biograph_mmr'}) % Verio
      gwtype = 16;
    elseif ~isempty(regexp(model,'vida'))
      gwtype = 17;
    elseif ~isempty(regexp(model,'cimax'))
      gwtype = 18;
    end
  elseif ~isempty(regexp(manu,'ge'))
    if ismember(model,{'discovery mr450','discovery mr750'})
      gwtype = 9;
    elseif ismember(model,{'discovery mr750w'})
      gwtype = 10;
    elseif ismember(model,{'signa premier'})
      gwtype = 14;
    elseif ismember(model,{'signa uhp'})
      gwtype = 15;
    elseif ismember(model,{'signa creator'}) && fieldstr == 3
      gwtype = 9;
    elseif ismember(model,{'brm'})
      gwtype = 2;
    elseif ismember(model,{'crm'})
      gwtype = 3;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function unwarpflag = get_unwarpflag(scan_info)
  unwarpflag = NaN;
  % NOTE: dcm2niix writes NonlinearGradientCorrection, and the Type field
  %   beside it, only for Siemens and GE; Philips sidecars carry neither, so an
  %   absent field is normal and must not error. Leaving unwarpflag = NaN reads
  %   as no additional unwarp required, which is correct for Philips, where the
  %   correction is applied on the scanner, and is the safe default otherwise
  if ~isfield(scan_info,'NonlinearGradientCorrection'), return; end
  if ~scan_info.NonlinearGradientCorrection
    unwarpflag = 0; % unwarp 3D
  else
    switch mmil_getfield(scan_info,'NonlinearGradientCorrectionType','')
      case '1D'
        unwarpflag = 2; % unwarp in-plane only
      case '2D'
        unwarpflag = 1; % unwarp through-plane only
      case '3D'
        unwarpflag = NaN; % no additional unwarp required
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isoctrflag = get_isoctrflag(scan_info)
  isoctrflag = NaN;
  if all(scan_info.TablePosition==0)    
    isoctrflag = 0;
  else
    isoctrflag = 1;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


