function info = cmig_filter_bids_info(scan_info,added)
%function info = cmig_filter_bids_info(scan_info,[added])
%
% Purpose: assemble the json sidecar struct for a processed derivative
%
%   Builds the output from a list of fields to COPY from the rawdata sidecar,
%   rather than copying everything and removing what is unwanted. A field not
%   on the list is simply not written, so a field added by a future dcm2niix
%   needs a deliberate decision before it appears in a derivative.
%
%   Four behaviors:
%     copy      take the rawdata value unchanged, when present
%     fill      populate a BIDS-named field from a vendor-specific counterpart
%               when the BIDS one is absent
%     rewrite   replace with the value true of the processed data
%     add       write a field computed by the pipeline
%   The last two both arrive through 'added', which is merged last and so
%   overrides anything copied.
%
% Required Input:
%   scan_info: struct loaded from the rawdata json sidecar
%
% Optional Input:
%   added: struct of pipeline-computed fields to merge in, overriding copies
%     used both for rewrites (e.g. SliceThickness at the processed value) and
%     for fields with no rawdata counterpart (e.g. registration_matrix_T1)
%     {default = []}
%
% Output:
%   info: struct to serialize as the derivative json sidecar
%
% Created:  08/07/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Last Mod: 08/26/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = [];
if ~mmil_check_nargs(nargin,1), return; end
if ~exist('added','var'), added = []; end

% fields copied from the rawdata sidecar when present
% NOTE: fields absent from a given sidecar are not written; several are
%   modality- or vendor-specific (TaskName func only, InversionTime anat only,
%   ParallelReductionFactor* on a minority of scans)
taglist = {...
  'AcquisitionDateTime' 'AcquisitionTime' 'B0FieldIdentifier'...
  'B0FieldSource' 'BodyPart' 'ConversionSoftware'...
  'ConversionSoftwareVersion' 'DeviceSerialNumber' 'DwellTime'...
  'EchoTime' 'EffectiveEchoSpacing' 'FlipAngle'...
  'IntendedFor' 'InversionTime' 'MRAcquisitionType'...
  'MagneticFieldStrength' 'Manufacturer' 'ManufacturersModelName'...
  'Modality' 'MultibandAccelerationFactor' 'NonlinearGradientCorrection'...
  'NonlinearGradientCorrectionType' 'ParallelReductionFactorInPlane'...
  'ParallelReductionFactorOutOfPlane'...
  'PartialFourier' 'PartialFourierDirection' 'PatientID'...
  'PatientName'...
  'PatientPosition' 'PhaseEncodingDirection' 'ProtocolName'...
  'PulseSequenceDetails' 'ReceiveCoilActiveElements' 'ReceiveCoilName'...
  'RepetitionTime' 'ScanOptions' 'ScanningSequence'...
  'SequenceName' 'SequenceVariant' 'SeriesDescription'...
  'SeriesInstanceUID' 'SeriesNumber' 'SeriesTime'...
  'SliceThickness' 'SoftwareVersions' 'SpacingBetweenSlices'...
  'SpoilingState' 'SpoilingType'...
  'StudyDate' 'StudyDescription' 'StudyInstanceUID'...
  'TaskName' 'TotalReadoutTime'...
};

% BIDS-named field <- vendor-specific counterpart, used when the BIDS one is
%   absent so that one key is populated for every vendor
% NOTE: dcm2niix writes ReceiveCoilName and SequenceName for Siemens only, and
%   CoilString and PulseSequenceName for GE and Philips only
filllist = {...
  'ReceiveCoilName' 'CoilString';...
  'SequenceName'    'PulseSequenceName';...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = struct();

% copy
for i=1:length(taglist)
  tag = taglist{i};
  if isfield(scan_info,tag)
    info.(tag) = scan_info.(tag);
  end
end

% fill from a vendor counterpart where the BIDS-named field is absent or empty
for i=1:size(filllist,1)
  tag = filllist{i,1}; src = filllist{i,2};
  if ~isfield(info,tag) || isempty(info.(tag))
    val = mmil_getfield(scan_info,src,[]);
    if ~isempty(val), info.(tag) = val; end
  end
end

% rewrites and pipeline-added fields, merged last so they override copies
if ~isempty(added) && isstruct(added)
  tags = fieldnames(added);
  for i=1:length(tags)
    info.(tags{i}) = added.(tags{i});
  end
end

return;
