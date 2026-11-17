#!/bin/bash

# cmig_proc_dwi.bash: dMRI processing of one ABCD or HBCD session: input
#   selection, gradient nonlinearity correction, a synthetic undistorted
#   b=0 from Synb0-DISCO (cmig_synb0_disco.sh), and TORTOISEProcess
#   (denoising, Gibbs ringing, motion, eddy current and susceptibility
#   distortion correction); run by cmig_proc_bids_dwi.py
#
# Created:  06/10/2026 by M. Okan Irfanoglu
# Prev mod: 10/28/2026 by Diliana Pecheva
# Prev mod: 11/02/2026 by Don Hagler
# Prev mod: 11/06/2026 by Don Hagler
# Prev mod: 11/11/2026 by Don Hagler
# Prev mod: 11/12/2026 by Don Hagler
# Last mod: 11/13/2026 by Don Hagler
#
# Renamed from ABCD_HBCD_DWI_processing.bash (abcd_preproc 81210b9)

# Output streams: tool messages on stderr join the log on stdout; stderr (fd 3)
# carries only fail's ERROR line, so it is empty unless the run fails.
exec 3>&2 2>&1


# ---------------------------------------------------------------------------
# BENCHMARKING HELPER
# Times a command and appends one TSV row per stage. Completely inert unless
# the caller exports BENCH_STAGE_TSV (e.g. the jobs written by
# make_dwi_jobs.py), so the script behaves exactly as before when run
# on its own.
# Execution semantics of the command are unchanged: same unquoted expansion,
# same exit status returned to the caller.
# ---------------------------------------------------------------------------
BENCH_STAGE_N=0

bench_TIME ()
{
  local bcmd="$1"
  local rc

  if [ -z "${BENCH_STAGE_TSV}" ] || [ ! -x /usr/bin/time ]
  then
      ${bcmd}
      return $?
  fi

  BENCH_STAGE_N=$((BENCH_STAGE_N + 1))
  local tmp="${BENCH_STAGE_TSV}.$$.tmp"
  local t0=`date +%Y-%m-%dT%H:%M:%S`

  /usr/bin/time -o "${tmp}" -f "%e\t%U\t%S\t%M\t%P\t%F\t%R\t%I\t%O" ${bcmd}
  rc=$?

  local tool=`basename "\`echo "${bcmd}" | awk '{print $1}'\`"`
  local flat=`echo "${bcmd}" | tr '\t\n' '  '`
  if [ -s "${tmp}" ]
  then
      printf "%d\t%s\t%s\t%d\t%s\t%s\n" \
             "${BENCH_STAGE_N}" "${t0}" "${tool}" "${rc}" \
             "`tail -n 1 ${tmp}`" "${flat}" >> "${BENCH_STAGE_TSV}"
  fi
  rm -f "${tmp}"
  return ${rc}
}

fail ()
{
  # report the line in the main script that led here, even via my_RUN
  local n=${#BASH_LINENO[@]}
  local msg="ERROR sub-${subj_id} ses-${ses_id} at `basename ${BASH_SOURCE[0]}` line ${BASH_LINENO[$((n-2))]} (script ${script_md5}): $1"
  echo ""
  echo "${msg}"
  echo "${msg}" >&3
  exit 1
}

my_RUN ()
{
  cmd="$1"
  echo "Running: ${cmd}"
  bench_TIME "${cmd}"
  local rc=$?
  [ ${rc} -eq 0 ] || fail "exit ${rc} from: ${cmd}"
}



echo ""


if [ $# -lt 2 ]
then
    echo "ABCD/HBCD processing script that combines several features from different pipelines. "
    echo "The script assumes the data is in BIDS convention..."
    echo ""
    echo "Usage: bash cmig_proc_dwi.bash  options"
    echo ""
    echo "COMMON options for both studies:"
    echo "REQUIRED......"
    echo "--study=name_of_study_data_to_be_processed (ABCD or HBCD)"
    echo "--project_root_folder=full_path_to_DWI_input_directory  or --prf=full_path_to_DWI_input_directory . The BIDS folder holding sub-*/ses-*; the DWI and reverse-PE data are read from it."
    echo "--sub_id=subject_id    . This is the subject id in BIDS format, without including the sub prefix."
    echo "--ses_id=session_id   . This is the subject id in BIDS format, without including the sess prefix."
    echo "--fs_license=FILE   . FreeSurfer licence file; every user needs their own (free registration at surfer.nmr.mgh.harvard.edu)."
    echo "-o=study_output_folder. This would be the folder in the derivatives containing all subjects' subfolders such as /xxx/derivatives/mmps_mproc or /xxx/derivatives/qsiprep"

    
    echo ""
    echo ""    
    echo "INPUT options (OPTIONAL) ......"
    echo "The script will try to automatically locate the images based on BIDS convention. To overwrite the inputs these optional parameters can be used."
    echo "    --run_id=run_id  . DEFAULT:1 . Without the prefix run."
    echo "    -u=full_path_to_the_main_NIFTI_file   OR --up_data=full_path_to_the_main_NIFTI_file (required)"
    echo "    -ub=main_data_bvals_file  OR --up_bvals=main_data_bvals_file  "
    echo "    -uv=main_data_bvecs_file  OR --up_bvecs=main_data_bvecs_file "    
    echo "    -u2=full_path_to_the_second_part_of main_NIFTI_file   OR --up_data2=full_path_to_the__second_part_of_main_NIFTI_file (required)"
    echo "    -ub2=part2_main_data_bvals_file  OR --up_bvals2=part2_main_data_bvals_file "
    echo "    -uv2=part2_main_data_bvecs_file  OR --up_bvecs2=part2_main_data_bvecs_file "    
    echo "    -d=full_path_to_the_reverse_PE_NIFTI_file   OR --down_data=full_path_to_the_reverse_PE_NIFTI_file"
    echo "    -db=down_data_bvals_file  OR --down_bvals=down_data_bvals_file"
    echo "    -dv=down_data_bvecs_file  OR --down_bvecs=down_data_bvecs_file "    
    echo " "
    echo "    -T1=full_path_to_T1W_structural_image  (OPTIONAL)"
    echo "    -T2=full_path_to_T2W_structural_image  (OPTIONAL)"
    echo "    --smri_dir=root_of_sMRI_processing_output. T1w/T2w are taken from {smri_dir}/sub-*/ses-*/anat (not from the DWI input directory) unless -T1/-T2 are given."
    echo "        Matches sub-X_ses-Y_run-NN_desc-preproc_T1w and sub-X_ses-Y_rec-norm_run-NN_desc-preproc_T1w (any run); more than one match stops the run."
    echo "    --useT1=(0/1). By default ABCD uses T1w and HBCD uses T2w anatomicals. This option, when set to 1, forces the pipeline to also use the T1W  (OPTIONAL)"    
    echo "    --useT2=(0/1). By default ABCD uses T1w and HBCD uses T2w anatomicals. This option, when set to 1, forces the pipeline to also use the T2W  (OPTIONAL)"    

    echo ""
    echo "Gradient nonlinearity options (OPTIONAL) ......"
    echo "--gw_coeffs=full_path_to_gradient_nonlinearity_coefficients_file (OPTIONAL). If provided, gradwarp will be applied and the gradient deviation matrix image L will be output in HCP style"
    echo "    --gw_field_full=full_path_to_3D_gradient_nonlinearity_field in ITK displacement field format. (OPTIONAL)"
    echo "    --gw_field_L=full_path_to_gradient_nonlinearity_field_L_fiel_in_MGH_format (OPTIONAL). All three must be provided if one of them is."
    echo "    --gw_field_P=full_path_to_gradient_nonlinearity_field_P_fiel_in_MGH_format (OPTIONAL). All three must be provided if one of them is."
    echo "    --gw_field_H=full_path_to_gradient_nonlinearity_field_H_fiel_in_MGH_format (OPTIONAL). All three must be provided if one of them is."
    echo "    --gw_mode=2d/3d/1d .  Flag to determine if the gradwarp will be applied in 2D or 3D or 1d (which will only be applied along the Z axis). Default:3d"
    echo "    --isocenter=0 or 1 .  Default:0 . When 1, the image center is assumed to be the isocenter of the magnet. When 0, the image header is used for the z axis as well."
    echo "    --gw_auto=0 or 1 .  Default:1 . When 1, gw_mode, isocenter and the gw_field_L/P/H files are derived from the DWI JSON sidecar by py_scripts/cmig_get_bids_gradwarpinfo.py, unless any of --gw_coeffs, --gw_field_full or --gw_field_L was given explicitly. Set 0 to disable gradwarp entirely unless given explicitly."

    echo ""
    echo "TEMPORARY DATA options (OPTIONAL) ..........:"
    echo "     -p=temporary_processing_folder OR --proc_dir=temporary_processing_folder. (OPTIONAL). If not provided, the script will create dwi/proc_dwi in the session's output folder"
    
    
    echo ""
    echo "OUTPUT options (OPTIONAL)...."
    echo "The output will be by default written in the BIDS derivatives folder. This can be overwritten..."
    echo "    --output_res=final_DWI_and_structural_resolution  Default:native_resolution." 
    echo "    --output_size=final_DWI_and_structural_size  Default:Size computed from structural FOV and desired resolution. Example usage: --output_size=150x160x140 "
    echo "    --output_orientation={LPS,RAI,LAI, etc..}  (Default:LPS , i.e. from Right-to-Left (from left of the image to right), Anterior-to-Posterior (from top of the image to bottom), Inferior-to-superior"
    
    
    echo ""
    echo " PROCESSING options (OPTIONAL).... "
    echo "--gpu=use_cuda_versions  .  0 or 1 .  Default:0"    
    echo "    Threads: OMP_NUM_THREADS if set, otherwise SLURM_CPUS_PER_TASK, NSLOTS, or nproc (never the whole node)."
    echo "    --denoising={0,1} (Default:1)"
    echo "    --gibbs={0,1}     (Default:1)"   
    
    echo " By M. Okan Irfanoglu 06/10/2026"
    echo " "
    
    exit 1            
fi

cwd=`pwd`
script_folder=$( dirname $( realpath $0 ))
script_md5=`md5sum "$0" | cut -c1-8`
echo SCRIPT FOLDER:  ${script_folder}


# TOOLS: FSL, FreeSurfer, ANTs, c3d, TORTOISE and Synb0-DISCO from the
# environment (the CMIG BIDS app image sets PATH, FSLDIR, ANTSPATH,
# FREESURFER_HOME and SYNB0_DISCO_DIR)
export FSLOUTPUTTYPE=NIFTI_GZ
[ -n "${SYNB0_DISCO_DIR}" ] || { echo "ERROR: SYNB0_DISCO_DIR not set" >&3 ; exit 1 ; }
[ -n "${CMIG_PIPELINE_ATLASDIR}" ] || { echo "ERROR: CMIG_PIPELINE_ATLASDIR not set" >&3 ; exit 1 ; }

py_scripts_folder=$( realpath ${script_folder}/../python )


# THREADS: keep what the launcher set; otherwise the scheduler allocation, not the node
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-${NSLOTS:-`nproc`}}}
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS:-${OMP_NUM_THREADS}}
echo "Threads: ${OMP_NUM_THREADS}"




#ALL THE MAIN VARIABLES

study=""
prf=""
subj_id=""
ses_id=""
run_id="1"
gpu=0
fs_license=''

up_nifti=''
up_bvals=''
up_bvecs=''

up_nifti2=''
up_bvals2=''
up_bvecs2=''


down_nifti=''
down_bvals=''
down_bvecs=''

T1W_structural=''
T2W_structural=''

useT1=''
useT2=''

gw_coeffs_file=''
gw_field_full=''
gw_field_L_file=''
gw_field_P_file=''
gw_field_H_file=''
gw_mode='3d'
isocenter=0
gw_auto=1

denoising=1
gibbs=1
output_res=""
output_size=""
output_orientation=""
output_dir=""

proc_dir=""
smri_dir=''


#INPUT PARSING

cmd=""

for i in "$@"
do
   cmd="$cmd $i"

case $i in
    --study=*)
        study="${i#*=}"
        study=`echo "${study,,}"`
        shift 
    ;;
    --project_root_folder=*|--prf=*)
        prf="${i#*=}"
        shift 
    ;;
    --sub_id=*)
        subj_id="${i#*=}"
        shift 
    ;;
    --ses_id=*)
        ses_id="${i#*=}"
        shift 
    ;;
    --run_id=*)
        run_id="${i#*=}"
        shift 
    ;;        
    --gpu=*)
        gpu="${i#*=}"
        shift 
    ;;       
    --fs_license=*)
        fs_license="${i#*=}"
        shift
    ;;
    -u=*|--up_data=*)
        up_nifti="${i#*=}"
        shift 
    ;;
    -ub=*|--up_bvals=*)
        up_bvals="${i#*=}"
        shift 
    ;;   
    -uv=*|--up_bvecs=*)
        up_bvecs="${i#*=}"
        shift 
    ;;   
    -u2=*|--up_data2=*)
        up_nifti2="${i#*=}"
        shift 
    ;;
    -ub2=*|--up_bvals2=*)
        up_bvals2="${i#*=}"
        shift 
    ;;   
    -uv2=*|--up_bvecs2=*)
        up_bvecs2="${i#*=}"
        shift 
    ;;  
    -d=*|--down_data=*)
        down_nifti="${i#*=}"
        shift 
    ;;
    -db=*|--down_bvals=*)
        down_bvals="${i#*=}"
        shift 
    ;;   
    -dv=*|--down_bvecs=*)
        down_bvecs="${i#*=}"
        shift 
    ;;    
    -T1=*)
        T1W_structural="${i#*=}"
        shift 
    ;;          
    -T2=*)
        T2W_structural="${i#*=}"
        shift 
    ;;   
    --useT1=*)
        useT1="${i#*=}"
        shift 
    ;;          
    --useT2=*)
        useT2="${i#*=}"
        shift 
    ;;       
    --gw_coeffs=*)
        gw_coeffs_file="${i#*=}"
        shift 
    ;;  
    --gw_field_full=*)
        gw_field_full="${i#*=}"
        shift 
    ;;   
    --gw_field_L=*)
        gw_field_L_file="${i#*=}"
        shift 
    ;;   
    --gw_field_P=*)
        gw_field_P_file="${i#*=}"
        shift 
    ;;   
    --gw_field_H=*)
        gw_field_H_file="${i#*=}"
        shift 
    ;;   
    --gw_mode=*)
        gw_mode="${i#*=}"
        shift 
    ;;   
    --isocenter=*)
        isocenter="${i#*=}"
        shift
    ;;
    --gw_auto=*)
        gw_auto="${i#*=}"
        shift
    ;;
    --denoising=*)
        denoising="${i#*=}"
        shift 
    ;;    
    --gibbs=*)
        gibbs="${i#*=}"
        shift 
    ;;    
    --use_T2=*)
        use_T2="${i#*=}"
        shift 
    ;;    
    --nos2v_norepol=*)
        nos2v_norepol="${i#*=}"
        shift 
    ;; 
    --output_res=*)
        output_res="${i#*=}"
        shift 
    ;;   
    --output_size=*)
        output_size="${i#*=}"
        shift 
    ;;  
    --output_orientation=*)
        output_orientation="${i#*=}"
        shift 
    ;;  
    --smri_dir=*)
        smri_dir="${i#*=}"
        shift
    ;;
    -p=*|--proc_dir=*)
        proc_dir="${i#*=}"
        shift 
    ;;   
    -o=*)
        output_dir="${i#*=}"
        shift 
    ;;   
    *)
        fail "unrecognized command line option ${i}"
    ;;
    *)
            # unknown option
    ;;
esac
done


if [[ ${study} != "abcd" &&  ${study} != "hbcd" ]]; then
    fail "study needs to be either ABCD or HBCD"
fi

# the infant synthstrip model HBCD's skull stripping uses is not in the image
#   (nor, so far, in any public source)
synthstrip_infant=${FREESURFER_HOME}/models/synthstrip.infant.1.pt
if [ ${study} == "hbcd" ] && [ ! -f ${synthstrip_infant} ]; then
    fail "HBCD needs synthstrip.infant.1.pt, not available (${synthstrip_infant})"
fi

if [ -z ${output_dir} ]
then
    fail "output folder must be entered"
fi


curr_subj_folder=${prf}/sub-${subj_id}/ses-${ses_id}
echo "Current subject folder: " ${curr_subj_folder}

if [ -z "${proc_dir}"  ]
then
    proc_dir=${output_dir}/sub-${subj_id}/ses-${ses_id}/dwi/proc_dwi
fi
mkdir -p ${proc_dir}
echo "Temporary processing folder: " ${proc_dir}

# temporary files of the tools (SynB0, ANTs) go inside proc_dir, so they are kept or cleaned with it
export TMPDIR=${proc_dir}/tmp
mkdir -p ${TMPDIR} || fail "creating ${TMPDIR}"

# FREESURFER LICENCE: read by FreeSurfer from FS_LICENSE; no licence ships with
# the pipeline, and it is not saved with the outputs
[ -n "${fs_license}" ] || fail "--fs_license=FILE is required (FreeSurfer licence)"
[ -f "${fs_license}" ] || fail "FreeSurfer licence ${fs_license} not found"
export FS_LICENSE=${fs_license}
trap 'exit 143' TERM INT

# identity transform for ApplyTransformationToScalar2
id_xfm=${proc_dir}/id.txt
printf '%s\n' '#Insight Transform File V1.0' '#Transform 0' 'Transform: AffineTransform_double_3_3' \
    'Parameters: 1 0 0 0 1 0 0 0 1 0 0 0' 'FixedParameters: 0 0 0' > ${id_xfm} || fail "writing ${id_xfm}"


# Decide on which anatomicals to use
if [ ${study} == "abcd" ]
then
    if [ -z ${useT1} ]
    then
        useT1=1
        if [ -z ${useT2} ]
        then
            useT2=0
        fi
    fi   
else
    if [ -z ${useT2} ]
    then
        useT2=1
        if [ -z ${useT1} ]
        then
            useT1=0
        fi
    fi   
fi


# FIND THE ANATOMICALS IN THE sMRI PROCESSING OUTPUT
# Matches sub-X_ses-Y_run-NN_desc-preproc_T1w and sub-X_ses-Y_rec-norm_run-NN_desc-preproc_T1w,
# any run, independent of the DWI run. More than one match is an error.
find_anat ()
{
    local suffix=$1
    local anat_dir=${smri_dir}/sub-${subj_id}/ses-${ses_id}/anat
    local matches=( $( find ${anat_dir} -maxdepth 1 -name "sub-${subj_id}_ses-${ses_id}*_desc-preproc_${suffix}.nii*" 2>/dev/null | sort ) )
    if [ ${#matches[@]} -gt 1 ]
    then
        fail "more than one ${suffix} in ${anat_dir}: ${matches[*]}"
    fi
    found_anat=${matches[0]}
}

if [ -z "${T1W_structural}" ] && [ -n "${smri_dir}" ]
then
    find_anat T1w
    T1W_structural=${found_anat}
fi
echo "Using as T1w data: " ${T1W_structural}

# the T2w is only needed when it will be used, or as the fallback when there is no T1w
if [ -z "${T2W_structural}" ] && [ -n "${smri_dir}" ] && { [ "${useT2}" == "1" ] || [ -z "${T1W_structural}" ]; }
then
    find_anat T2w
    T2W_structural=${found_anat}
fi
echo "Using as T2w data: " ${T2W_structural}


if [[ -z ${T1W_structural} && -z ${T2W_structural}  ]]; then
    fail "no T1W or T2W anatomical images were found or entered as a parameter"
fi    


if [ ! -z ${T1W_structural} ]
then
    if [ ! -e "${T1W_structural}" ]
    then
        fail "T1W_structural ${T1W_structural} does not exist"
    fi
fi

if [ ! -z ${T2W_structural} ]
then
    if [ ! -e "${T2W_structural}" ]
    then
        fail "T2W_structural ${T2W_structural} does not exist"
    fi
fi


if [[ -z ${T2W_structural} &&  ! -z ${T1W_structural} ]];then
    useT2=0
    useT1=1
fi   
if [[ -z ${T1W_structural} &&  ! -z ${T2W_structural} ]];then
    useT2=1
    useT1=0
fi   

echo "Use the T1W image for procesing: " ${useT1}
echo "Use the T2W image for procesing: " ${useT2}


#DECIDE ON THE INPUTS BASED ON PARAMETERS

# ABCD with no -u/-d: choose the runs and field map from scans.tsv (usable runs, shim match)
if [ ${study} == "abcd" ] && [ -z "${up_nifti}" ] && [ -z "${down_nifti}" ]
then
    sel=`python3 ${py_scripts_folder}/cmig_select_bids_dwi_runs.py ${curr_subj_folder}` \
        || fail "selecting the DWI runs: ${sel}"
    eval "${sel}"
    echo "Selected: up ${up_nifti} ${up_nifti2}, down ${down_nifti} (by ${fmap_method}), shim match '${shim_match}'"
fi

if [ -z "${up_nifti}" ]
then
   if [ ! -e ${curr_subj_folder} ]
   then
       fail "${curr_subj_folder}/dwi folder does not exist"
   fi
   
   phase=""
   if [ ${study} == "hbcd" ]
   then 
       phase="AP"
   fi
   echo ${curr_subj_folder}  
   
   up_nifti=$( find ${curr_subj_folder}/dwi -maxdepth 1 -name "sub-${subj_id}_ses-${ses_id}*${phase}*_run*${run_id}*_dwi.nii.gz"  ! -name "*run*run*" )    
  
   if [ -z ${up_nifti} ]
   then
       up_nifti=$( find ${curr_subj_folder}/dwi -maxdepth 1 -name "sub-${subj_id}_ses-${ses_id}*${phase}*_run*${run_id}*_dwi.nii"  ! -name "*run*run*" )     
   fi
   up_nifti=$( realpath ${up_nifti} )
fi
echo "Using as up data: " ${up_nifti}

if [ -z "${down_nifti}" ]
then
    if [ ${study} == "abcd" ]
    then
        phase="dwi_dir-AP"
        down_nifti=$( find ${curr_subj_folder}/fmap -maxdepth 1 -name "sub-${subj_id}_ses-${ses_id}*${phase}*_run*${run_id}*_epi.nii.gz"  ! -name "*run*run*" ) 
        if [ -z ${down_nifti} ]
        then
            down_nifti=$( find ${curr_subj_folder}/fmap -maxdepth 1 -name "sub-${subj_id}_ses-${ses_id}*${phase}*_run*${run_id}*_epi.nii" ! -name "*run*run*" ) 
        fi
    else
        phase="PA"
        down_nifti=$( find ${curr_subj_folder}/dwi -maxdepth 1 -name "sub-${subj_id}_ses-${ses_id}*${phase}*_run*${run_id}*_dwi.nii.gz" ! -name "*run*run*" ) 
        if [ -z ${down_nifti} ]
        then
            down_nifti=$( find ${curr_subj_folder}/dwi -maxdepth 1 -name "sub-${subj_id}_ses-${ses_id}*${phase}*_run*${run_id}*_dwi.nii" ! -name "*run*run*" ) 
        fi
    fi
    down_nifti=$( realpath ${down_nifti} )
fi
echo "Using as down data: " ${down_nifti}

#INPUT ERROR CHECKS

if [ ! -e "${up_nifti}" ]
then
    fail "up NIFTI ${up_nifti} does not exist"
fi
if [ ! -e "${down_nifti}" ]
then
    fail "down NIFTI ${down_nifti} does not exist"
fi


if [[ "$up_nifti" == *".nii.gz"* ]]; then
  up_ext=nii.gz
else
  up_ext=nii  
fi
if [[ "$down_nifti" == *".nii.gz"* ]]; then
  down_ext=nii.gz
else
  down_ext=nii  
fi





# BIDS URI of a rawdata input (bids:rawdata:, the DatasetLinks name of the input root),
# from its file name: datatype fmap for an _epi, else dwi
raw_uri ()
{
    local base=`basename $1` dt=dwi
    [[ ${base} == *_epi.nii* ]] && dt=fmap
    echo bids:rawdata:sub-${subj_id}/ses-${ses_id}/${dt}/${base}
}
up_raw_uri=`raw_uri ${up_nifti}`
down_raw_uri=`raw_uri ${down_nifti}`
up2_raw_uri=""

# STAGE THE INPUTS INTO proc_dir
# The DWI input directory's dwi/fmap folders are read-only, so the original in-place
# bvals/bvecs handling that lived here could not write next to the NIFTIs.
# Everything the pipeline touches is copied into proc_dir first.
up_src=`dirname ${up_nifti}`
up_base=`basename ${up_nifti} .${up_ext}`
install -m 664 ${up_nifti} ${up_src}/${up_base}.json ${proc_dir}/ || fail "staging ${up_nifti}"
install -m 664 ${up_bvals:-${up_src}/${up_base}.bval} ${proc_dir}/${up_base}.bval || fail "staging up bval"
install -m 664 ${up_bvecs:-${up_src}/${up_base}.bvec} ${proc_dir}/${up_base}.bvec || fail "staging up bvec"
up_nifti=${proc_dir}/`basename ${up_nifti}`
up_bvals=${proc_dir}/${up_base}.bval
up_bvecs=${proc_dir}/${up_base}.bvec

# DERIVE THE GRADIENT NONLINEARITY SETTINGS FROM THE JSON SIDECAR
# Skipped when the caller supplied a gradwarp source explicitly, so anything
# passed on the command line still wins. cmig_get_bids_gradwarpinfo.py emits KEY=VALUE
# (gw_mode, isocenter, gw_field_L/P/H) and prints nothing when no gradwarp
# applies -- Philips data, or a scan the scanner already corrected in 3D.
if [ "${gw_auto}" -eq 1 ] && [ -z "${gw_coeffs_file}" ] && [ -z "${gw_field_full}" ] && [ -z "${gw_field_L_file}" ]
then
    gw_derived=$( python3 ${py_scripts_folder}/cmig_get_bids_gradwarpinfo.py ${proc_dir}/${up_base}.json )
    if [ $? -ne 0 ]
    then
        fail "could not determine gradwarp settings from ${proc_dir}/${up_base}.json"
    elif [ -n "${gw_derived}" ]
    then
        eval "${gw_derived}"
        gw_field_L_file=${gw_field_L}
        gw_field_P_file=${gw_field_P}
        gw_field_H_file=${gw_field_H}
        echo "Gradwarp derived from sidecar: mode=${gw_mode} isocenter=${isocenter}"
        echo "Gradwarp fields: `basename ${gw_field_L_file}` `basename ${gw_field_P_file}` `basename ${gw_field_H_file}`"
    else
        echo "No gradwarp correction applies to this scan"
    fi
fi

down_src=`dirname ${down_nifti}`
down_base=`basename ${down_nifti} .${down_ext}`
# run entity of the reverse-PE scan, for the B0dx name
down_run=`echo ${down_base} | grep -o '_run-[0-9]*_' | tr -d _`
[ -n "${down_run}" ] || fail "no run entity in ${down_base}"
install -m 664 ${down_nifti} ${proc_dir}/down.${down_ext} || fail "staging ${down_nifti}"
install -m 664 ${down_src}/${down_base}.json ${proc_dir}/down.json || fail "staging ${down_src}/${down_base}.json"
down_nifti=${proc_dir}/down.${down_ext}
# bval/bvec of the down series, as proc_dwi/down.bval/.bvec (TORTOISE reads them beside down.nii):
# -db/-dv if given, else beside the NIfTI, else zeros for a BIDS field map (b=0 only)
[ -z "${down_bvals}" ] && [ -e ${down_src}/${down_base}.bval ] && down_bvals=${down_src}/${down_base}.bval
[ -z "${down_bvecs}" ] && [ -e ${down_src}/${down_base}.bvec ] && down_bvecs=${down_src}/${down_base}.bvec
if [ -z "${down_bvals}" ] && [ -z "${down_bvecs}" ] && [[ ${down_base} == *_epi ]]
then
    down_nvol=`fslinfo ${down_nifti} | awk '$1=="dim4"{print $2}'`
    [ -n "${down_nvol}" ] || fail "could not read the number of volumes of ${down_nifti}"
    zeros=""
    for ((i=1; i<=${down_nvol}; i+=1)); do zeros="0 ${zeros}"; done
    echo ${zeros} > ${proc_dir}/down.bval || fail "staging down bval"
    { echo ${zeros}; echo ${zeros}; echo ${zeros}; } > ${proc_dir}/down.bvec || fail "staging down bvec"
else
    [ -n "${down_bvals}" ] && [ -n "${down_bvecs}" ] \
        || fail "no bval/bvec for the down series ${down_base}; give them with -db/-dv"
    install -m 664 ${down_bvals} ${proc_dir}/down.bval || fail "staging down bval"
    install -m 664 ${down_bvecs} ${proc_dir}/down.bvec || fail "staging down bvec"
fi
down_bvals=${proc_dir}/down.bval
down_bvecs=${proc_dir}/down.bvec

for f in ${up_bvals} ${up_bvecs} ${down_bvals} ${down_bvecs}
do
    if [ ! -e ${f} ]
    then
        fail "missing after staging: ${f}"
    fi
done

echo "Using as up data:   " ${up_nifti} ${up_bvals} ${up_bvecs}
echo "Using as down data: " ${down_nifti} ${down_bvals} ${down_bvecs}





# CHECK IF ABCD AND PHILIPS or Prisma
up_json=`dirname ${up_nifti}`/`basename ${up_nifti} .${up_ext}`.json
sidecar_vals=`python3 ${py_scripts_folder}/cmig_read_bids_sidecar.py ${up_json} Manufacturer ManufacturersModelName` \
    || fail "reading ${up_json}"
eval "${sidecar_vals}"
isABCDPhilips=0
if [ ${study} == "abcd" ] && [[ "${Manufacturer,,}" == *"philips"* ]]; then
    isABCDPhilips=1
fi
isPrisma=0
if [[ "${ManufacturersModelName,,}" == *"prisma"* ]]; then
    isPrisma=1
fi



# run entity of the input DWI (run-01 for Philips), for the output name
up_run=`basename ${up_nifti} | grep -o '_run-[0-9]*_' | tr -d _`
[ -n "${up_run}" ] || fail "no run entity in `basename ${up_nifti}`"

# staged run-01 (or only) DWI, the source image in its sidecar
up_nifti1=${up_nifti}

# IF ABCD and Philips handle the second series
if [ $isABCDPhilips -eq 1 ] 
then 
    echo "Handling ABCD Philips data"
    phase=""

    if [  -z "${up_nifti2}"  ]
    then
        up_nifti2=$( realpath $( find ${curr_subj_folder}/dwi -maxdepth 1 -name "sub-${subj_id}_ses-${ses_id}_run-02*_dwi.nii*"  ) )                               
    fi
    
    if [[ "$up_nifti2" == *".nii.gz"* ]]; then
        up_ext2=nii.gz
    else
        up_ext2=nii
    fi

    # STAGE RUN-02 INTO proc_dir, as for run-01 above. FSLBVecsToTORTOISEBmatrix
    # writes its .bmtxt beside the bval/bvec, so they must not stay in the DWI input directory.
    up2_src=`dirname ${up_nifti2}`
    up2_base=`basename ${up_nifti2} .${up_ext2}`
    up2_raw_uri=`raw_uri ${up_nifti2}`
    install -m 664 ${up_nifti2} ${proc_dir}/ || fail "staging ${up_nifti2}"
    install -m 664 ${up2_src}/${up2_base}.json ${proc_dir}/ || fail "staging ${up2_src}/${up2_base}.json"
    up2_run=`echo ${up2_base} | grep -o '_run-[0-9]*_' | tr -d _`
    [ -n "${up2_run}" ] || fail "no run entity in ${up2_base}"
    install -m 664 ${up_bvals2:-${up2_src}/${up2_base}.bval} ${proc_dir}/${up2_base}.bval || fail "staging run-02 bval"
    install -m 664 ${up_bvecs2:-${up2_src}/${up2_base}.bvec} ${proc_dir}/${up2_base}.bvec || fail "staging run-02 bvec"
    up_nifti2=${proc_dir}/`basename ${up_nifti2}`

    up1_nvol=`fslinfo ${up_nifti} | grep ^dim4 | awk '{print $2}'`
    up2_nvol=`fslinfo ${up_nifti2} | grep ^dim4 | awk '{print $2}'`
    [ -n "${up1_nvol}" ] || fail "could not read the number of volumes of ${up_nifti}"
    [ -n "${up2_nvol}" ] || fail "could not read the number of volumes of ${up_nifti2}"
        
    my_RUN "FSLBVecsToTORTOISEBmatrix `dirname ${up_nifti}`/`basename ${up_nifti} .${up_ext}`.bval   `dirname ${up_nifti}`/`basename ${up_nifti} .${up_ext}`.bvec"
    my_RUN "FSLBVecsToTORTOISEBmatrix `dirname ${up_nifti2}`/`basename ${up_nifti2} .${up_ext2}`.bval   `dirname ${up_nifti2}`/`basename ${up_nifti2} .${up_ext2}`.bvec"
    
        
    cp `dirname ${up_nifti}`/`basename ${up_nifti} .${up_ext}`.json  `dirname ${up_nifti}`/sub-${subj_id}_ses-${ses_id}_run-01_run-02_dwi.json || fail "copying the run-01 json for the combined Philips series"
    
    combined_data_name=`dirname ${up_nifti}`/sub-${subj_id}_ses-${ses_id}_run-01_run-02_dwi.nii
    my_RUN "CombineDWIsWithBmatrix ${combined_data_name} ${up_nifti} ${up_nifti2}"
    up_nifti=${combined_data_name}
    up_ext=nii           
fi



# OUTPUT THE ENTERED COMMAND TO LOG
cmd_file=${proc_dir}/cmd.log
if [  -e "${cmd_file}" ]
then
   rm ${cmd_file}
fi
echo ${cmd} > ${cmd_file} 








# PROCESS THE STRUCTURAL IMAGES
# first reorient to LPS axial
echo "Processing the anatomicals..."
if [ ! -z ${T1W_structural}   ]
then
    T1W_smri=${T1W_structural}
    install -m 664 ${T1W_structural} ${proc_dir}/ || fail "staging ${T1W_structural}"
    # the sMRI T1w's name without _T1w.nii[.gz] (e.g. sub-X_ses-Y_rec-norm_run-01), for its derived files
    t1_stem=`basename ${T1W_structural} | sed 's/_T1w\.nii\(\.gz\)\{0,1\}$//'`
    [ "${t1_stem}" != "`basename ${T1W_structural}`" ] || fail "`basename ${T1W_structural}` does not end in _T1w.nii[.gz]"
    T1W_structural=${proc_dir}/`basename ${T1W_structural}`
    oname=${proc_dir}/${t1_stem}_T1w_reoriented.nii.gz
    my_RUN "ReorientImage3D -i ${T1W_structural} -o ${oname}"
    T1W_structural=${oname}
fi
if [ "${useT2}" == "1" ]
then
    T2W_smri=${T2W_structural}
    install -m 664 ${T2W_structural} ${proc_dir}/ || fail "staging ${T2W_structural}"
    # the sMRI T2w's name without _T2w.nii[.gz] (e.g. sub-X_ses-Y_rec-norm_run-01), for its derived files
    t2_stem=`basename ${T2W_structural} | sed 's/_T2w\.nii\(\.gz\)\{0,1\}$//'`
    [ "${t2_stem}" != "`basename ${T2W_structural}`" ] || fail "`basename ${T2W_structural}` does not end in _T2w.nii[.gz]"
    T2W_structural=${proc_dir}/`basename ${T2W_structural}`
    oname=${proc_dir}/${t2_stem}_T2w_reoriented.nii.gz
    my_RUN "ReorientImage3D -i ${T2W_structural} -o ${oname}"
    T2W_structural=${oname}    
fi

#then brain mask
if [ ! -z ${T1W_structural}  ]
then
    if [ ${study} == "abcd" ]
    then
        cmd="mri_synthstrip -i ${T1W_structural} -o ${proc_dir}/${t1_stem}_T1w_reoriented_brain.nii.gz -m ${proc_dir}/${t1_stem}_T1w_reoriented_brain_mask.nii.gz"
        my_RUN "${cmd}"
    else
        cmd="mri_synthstrip -i ${T1W_structural} -o ${proc_dir}/${t1_stem}_T1w_reoriented_brain.nii.gz -m ${proc_dir}/${t1_stem}_T1w_reoriented_brain_mask.nii.gz --model ${synthstrip_infant}"
        my_RUN "${cmd}"
    fi
fi
if [ "${useT2}" == "1" ]
then
    if [ ${study} == "abcd" ]
    then
        cmd="mri_synthstrip -i ${T2W_structural} -o ${proc_dir}/${t2_stem}_T2w_reoriented_brain.nii.gz -m ${proc_dir}/${t2_stem}_T2w_reoriented_brain_mask.nii.gz"
        my_RUN "${cmd}"
    else
        cmd="mri_synthstrip -i ${T2W_structural} -o ${proc_dir}/${t2_stem}_T2w_reoriented_brain.nii.gz -m ${proc_dir}/${t2_stem}_T2w_reoriented_brain_mask.nii.gz --model ${synthstrip_infant}"
        my_RUN "${cmd}"
    fi
fi


# for T2s, do smooth nonbrain tissue transition
if [ "${useT2}" == "1" ]
then
    if [ ${study} == "abcd" ]
    then
        midslice=100
    else
        midslice=50
    fi
    my_RUN "SmoothTransitionHCPStructural ${T2W_structural} ${midslice} ${proc_dir}/${t2_stem}_T2w_reoriented_smoothtrans.nii.gz ${proc_dir}/${t2_stem}_T2w_reoriented_brain.nii.gz "
        
fi

if [ ${useT1} -eq 1 ]
#if [ ${study} == "abcd" ]
then
    echo "ABCD/HBCD MESSAGE: performing T1 to T2 conversion initial steps"
    mkdir -p ${proc_dir}/docker
    mkdir -p ${proc_dir}/docker/INPUTS
    mkdir -p ${proc_dir}/docker/OUTPUTS
    
    my_RUN "SelectBestB0 ${up_nifti} ${up_bvals} ${proc_dir}/docker/INPUTS/b0_lowres.nii"
    cp ${T1W_structural}  ${proc_dir}/docker/INPUTS/T1.nii.gz || fail "copying the T1w into ${proc_dir}/docker/INPUTS"
    cp ${proc_dir}/${t1_stem}_T1w_reoriented_brain.nii.gz ${proc_dir}/docker/INPUTS/T1_brain.nii.gz || fail "copying the T1w brain into ${proc_dir}/docker/INPUTS"
    
    #handle gradient nonlinearities
    
    current_gradwarp_field=""
    if [ ! -z "${gw_coeffs_file}"  ]
    then
        gradnonlin_parameter="${gradnonlin_parameter} --gradnonlin_coeff_file ${gw_coeffs_file}"
        
        echo "Creating gradient nonlinearity displacement field..."
        my_RUN "CreateNonlinearityDisplacementMap ${gw_coeffs_file}  ${proc_dir}/docker/INPUTS/b0_lowres.nii  ${proc_dir}/field_inv.nii  ${isocenter}"
        current_gradwarp_field_inv=${proc_dir}/field_inv.nii    
                
    else
        if [ ! -z "${gw_field_full}"  ]
        then
            current_gradwarp_field_inv=${gw_field_full}                            
        else
	        if [ ! -z "${gw_field_L_file}"  ]
	        then	    
	            echo "Creating gradient nonlinearity displacement field..."
		        echo "ABCD MESSAGE:  CONVERTING THE FREESURFER FIELDS TO TORTOISE FIELDS"
        		my_RUN "mri_convert ${gw_field_L_file}  ${proc_dir}/field_L.nii"
	        	my_RUN "mri_convert ${gw_field_P_file}  ${proc_dir}/field_P.nii"
		        my_RUN "mri_convert ${gw_field_H_file}  ${proc_dir}/field_H.nii"     
			            
        		my_RUN "CreateDisplacementField 3 0 ${proc_dir}/field_L.nii ${proc_dir}/field_P.nii ${proc_dir}/field_H.nii ${proc_dir}/field_inv.nii"
        		current_gradwarp_field_inv=${proc_dir}/field_inv.nii						  		        		                  
	        fi
        fi
    fi    
    if [ !  -z "${current_gradwarp_field_inv}"  ]   
    then
        temp_ge_field=${proc_dir}/`basename ${current_gradwarp_field_inv} .nii`_ge.nii
        cp ${current_gradwarp_field_inv} ${temp_ge_field} || fail "copying the gradwarp field ${current_gradwarp_field_inv}"
            
        if [ ${isocenter} -eq 1 ]
        then      
            echo "Setting gradient nonlinearity field to isocenter"  
            my_RUN "SetFieldToIsocenter ${temp_ge_field} ${proc_dir}/docker/INPUTS/b0_lowres.nii ${temp_ge_field}"
            if [[ "${gw_mode}" == *"2d"* ]]
            then 
                echo "Setting gradient nonlinearity field to 2D"  
                my_RUN "ZeroFieldDims ${temp_ge_field} 2 ${temp_ge_field}"
            fi    
            
            if [[ "${gw_mode}" == *"1d"* ]]
            then 
                echo "Setting gradient nonlinearity field to 1D"  
                my_RUN "ZeroFieldDims ${temp_ge_field} 0,1 ${temp_ge_field}"
            fi       
        else
            if [[ "${gw_mode}" == *"2d"* ]]
            then 
                echo "Setting gradient nonlinearity field to 2D"  
                my_RUN "ZeroFieldDims ${temp_ge_field} 2 ${temp_ge_field}"
            fi    
            if [[ "${gw_mode}" == *"1d"* ]]
            then 
                echo "Setting gradient nonlinearity field to 1D"  
                my_RUN "ZeroFieldDims ${temp_ge_field} 0,1 ${temp_ge_field}"
            fi          
        fi   
        current_gradwarp_field_inv=${temp_ge_field}
    fi
    
    echo "ABCD/HBCD MESSAGE:  RIGIDLY REGISTERING THE B0 IMAGE TO T1W IMAGE"
    # epi_reg starts from b=0 -> raw T1w (headers, same scanner frame) then raw T1w -> sMRI T1w
    # (rigid FLIRT), not from its own search, which can settle on a wrong fit
    t1_run=`basename ${T1W_smri} | grep -o '_run-[0-9]*_'`
    raw_t1w=`ls ${curr_subj_folder}/anat/sub-${subj_id}_ses-${ses_id}*${t1_run}T1w.nii* 2>/dev/null`
    [ `echo ${raw_t1w} | wc -w` -eq 1 ] || fail "need one raw T1w${t1_run}in ${curr_subj_folder}/anat to initialise epi_reg, found: ${raw_t1w:-none}"
    reg_dir=${proc_dir}/docker/INPUTS
    my_RUN "flirt -in ${reg_dir}/b0_lowres.nii -ref ${raw_t1w} -usesqform -applyxfm -omat ${reg_dir}/b0_to_rawT1.mat -out ${reg_dir}/b0_on_rawT1.nii.gz"
    my_RUN "flirt -in ${raw_t1w} -ref ${reg_dir}/T1.nii.gz -usesqform -applyxfm -omat ${reg_dir}/rawT1_hdr.mat -out ${reg_dir}/rawT1_hdr.nii.gz"
    my_RUN "flirt -in ${raw_t1w} -ref ${reg_dir}/T1.nii.gz -dof 6 -init ${reg_dir}/rawT1_hdr.mat -searchrx -30 30 -searchry -30 30 -searchrz -30 30 -omat ${reg_dir}/rawT1_to_T1.mat -out ${reg_dir}/rawT1_on_T1.nii.gz"
    my_RUN "convert_xfm -omat ${reg_dir}/b0_start.mat -concat ${reg_dir}/rawT1_to_T1.mat ${reg_dir}/b0_to_rawT1.mat"
    my_RUN "epi_reg --epi=${reg_dir}/b0_lowres.nii --t1=${reg_dir}/T1.nii.gz --t1brain=${reg_dir}/T1_brain.nii.gz --out=${reg_dir}/b0.nii.gz --init=${reg_dir}/b0_start.mat"
    # epi_reg exits 0 when it fails
    [ -s ${reg_dir}/b0.mat ] || fail "epi_reg produced no b0.mat (its errors are in the output above)"

    # check the fit by its BBR cost (no optimisation); above bbr_bad it is likely misregistered
    # (e.g. the head moved between the T1w's and the DWI's studies), so also register through the
    # atlases (b=0 -> atlas T2, T1w -> atlas T1, combined and reduced to rigid, then BBR) and keep
    # the fit with the lower cost
    bbr_bad=0.72
    bbr_cost () { flirt -in ${reg_dir}/b0_lowres.nii -ref ${reg_dir}/T1.nii.gz -wmseg ${reg_dir}/b0_fast_wmseg.nii.gz -cost bbr -init $1 -schedule ${FSLDIR}/etc/flirtsch/measurecost1.sch -omat /dev/null | head -1 | awk '{print $1}'; }
    reg_cost=`bbr_cost ${reg_dir}/b0.mat`
    [ -n "${reg_cost}" ] || fail "measuring the BBR cost of epi_reg's fit"
    reg_route=raw_t1w
    echo "ABCD/HBCD MESSAGE: BBR cost of the b=0 -> T1w fit: ${reg_cost}"
    if awk "BEGIN {exit !(${reg_cost} > ${bbr_bad})}"
    then
        echo "ABCD/HBCD MESSAGE: BBR cost above ${bbr_bad}, registering through the atlases"
        atlas=${CMIG_PIPELINE_ATLASDIR}/T1T2
        my_RUN "python3 ${py_scripts_folder}/cmig_dwi_atlas_b0_start.py com ${reg_dir}/b0_lowres.nii ${atlas}/atlasT2.nii.gz ${atlas}/maskT2.nii.gz ${reg_dir}/b0_to_atlas_com.mat"
        my_RUN "flirt -in ${reg_dir}/b0_lowres.nii -ref ${atlas}/atlasT2.nii.gz -refweight ${atlas}/maskT2.nii.gz -dof 9 -init ${reg_dir}/b0_to_atlas_com.mat -omat ${reg_dir}/b0_to_atlas.mat"
        my_RUN "flirt -in ${reg_dir}/T1.nii.gz -ref ${atlas}/atlasT1.nii.gz -usesqform -applyxfm -omat ${reg_dir}/T1_to_atlas_hdr.mat"
        my_RUN "flirt -in ${reg_dir}/T1.nii.gz -ref ${atlas}/atlasT1.nii.gz -inweight ${proc_dir}/${t1_stem}_T1w_reoriented_brain_mask.nii.gz -dof 12 -init ${reg_dir}/T1_to_atlas_hdr.mat -omat ${reg_dir}/T1_to_atlas.mat"
        my_RUN "convert_xfm -omat ${reg_dir}/atlas_to_T1.mat -inverse ${reg_dir}/T1_to_atlas.mat"
        my_RUN "convert_xfm -omat ${reg_dir}/b0_atlas_affine.mat -concat ${reg_dir}/atlas_to_T1.mat ${reg_dir}/b0_to_atlas.mat"
        my_RUN "python3 ${py_scripts_folder}/cmig_dwi_atlas_b0_start.py rigid ${reg_dir}/b0_atlas_affine.mat ${reg_dir}/b0_lowres.nii ${reg_dir}/b0_atlas_start.mat"
        my_RUN "flirt -in ${reg_dir}/b0_lowres.nii -ref ${reg_dir}/T1.nii.gz -dof 6 -cost bbr -wmseg ${reg_dir}/b0_fast_wmseg.nii.gz -init ${reg_dir}/b0_atlas_start.mat -omat ${reg_dir}/b0_atlas.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch"
        atlas_cost=`bbr_cost ${reg_dir}/b0_atlas.mat`
        [ -n "${atlas_cost}" ] || fail "measuring the BBR cost of the atlas fit"
        echo "ABCD/HBCD MESSAGE: BBR cost of the atlas fit: ${atlas_cost}"
        if awk "BEGIN {exit !(${atlas_cost} < ${reg_cost})}"
        then
            # SynB0-DISCO reads the b=0 in T1w space (b0.nii.gz), so rewrite it from the atlas fit
            mv ${reg_dir}/b0.mat ${reg_dir}/b0_raw_t1w.mat || fail "keeping epi_reg's fit as b0_raw_t1w.mat"
            cp ${reg_dir}/b0_atlas.mat ${reg_dir}/b0.mat || fail "copying the atlas fit to b0.mat"
            my_RUN "applywarp -i ${reg_dir}/b0_lowres.nii -r ${reg_dir}/T1.nii.gz -o ${reg_dir}/b0.nii.gz --premat=${reg_dir}/b0.mat --interp=spline"
            reg_cost=${atlas_cost}
            reg_route=atlas
        fi
        awk "BEGIN {exit !(${reg_cost} > ${bbr_bad})}" && echo "ABCD/HBCD WARNING: b=0 -> T1w fit still has BBR cost ${reg_cost} (above ${bbr_bad}): likely misregistered"
    fi
    echo "ABCD/HBCD MESSAGE: b=0 -> T1w fit from the ${reg_route} route, BBR cost ${reg_cost}"
    
    #T1 to T2 conversion with SynB0
    echo "ABCD MESSAGE: SYNB0-DISCO"
    # SynB0-DISCO carries on after a failed step and still exits 0, so check for its result;
    # clear its outputs first so a b0_u.nii.gz left from an earlier run cannot pass the check
    rm -rf ${proc_dir}/docker/OUTPUTS || fail "clearing ${proc_dir}/docker/OUTPUTS"
    mkdir -p ${proc_dir}/docker/OUTPUTS || fail "creating ${proc_dir}/docker/OUTPUTS"
    my_RUN "bash ${script_folder}/cmig_synb0_disco.sh ${proc_dir}/docker"
    [ -s ${proc_dir}/docker/OUTPUTS/b0_u.nii.gz ] || fail "SynB0-DISCO produced no b0_u.nii.gz (its errors are in the output above)"

    if [  -z "${current_gradwarp_field_inv}"  ]   
    then
        gunzip -f ${proc_dir}/docker/OUTPUTS/b0_u.nii.gz || fail "unzipping the SynB0 output b0_u.nii.gz"
        cp ${proc_dir}/docker/OUTPUTS/b0_u.nii ${proc_dir}/T2W_synth_str.nii || fail "copying the SynB0 output b0_u.nii"
    else    
        my_RUN "ApplyTransformationToScalar2 ${proc_dir}/docker/OUTPUTS/b0_u.nii.gz  ${current_gradwarp_field_inv}  ${proc_dir}/T2W_synth_str.nii  ${proc_dir}/docker/OUTPUTS/b0_u.nii.gz   BSP"
    fi
    
    cp ${proc_dir}/T2W_synth_str.nii ${proc_dir}/T2W_synth_input.nii || fail "copying T2W_synth_str.nii"
    cmd="mri_synthstrip -i ${proc_dir}/T2W_synth_input.nii -o ${proc_dir}/T2W_synth_input_brain.nii -m ${proc_dir}/T2W_synth_input_brain_mask.nii "
    my_RUN "${cmd}"
    my_RUN "ImageMath 3 ${proc_dir}/T2W_synth_input_skull.nii  - ${proc_dir}/T2W_synth_input.nii  ${proc_dir}/T2W_synth_input_brain.nii"
    my_RUN "ImageMath 3 ${proc_dir}/T2W_synth_input_skullm.nii m ${proc_dir}/T2W_synth_input_skull.nii 0.7"
    my_RUN "ImageMath 3 ${proc_dir}/T2W_synth_input_brainm.nii + ${proc_dir}/T2W_synth_input_brain.nii ${proc_dir}/T2W_synth_input_skullm.nii"
    
    structural_for_TORTOISE=${proc_dir}/T2W_synth_input_brainm.nii
else
    structural_for_TORTOISE=${proc_dir}/${t2_stem}_T2w_reoriented_smoothtrans.nii.gz
fi



#TORTOISE PART
echo "ABCD/HBCD MESSAGE: TORTOISE correction"

if [ ! -z "${gw_coeffs_file}"  ]
then
    gradnonlin_param="--grad_nonlin  ${gw_coeffs_file}[${isocenter},${gw_mode^^}]"
else
    if [ ! -z "${gw_field_full}"  ]
    then
        gradnonlin_param="--grad_nonlin  ${gw_field_full}[${isocenter},${gw_mode^^}]"                           
    else
	    if [ ! -z "${gw_field_L_file}"  ]
	    then	      		
    		gradnonlin_param="--grad_nonlin  ${proc_dir}/field_inv.nii[${isocenter},${gw_mode^^}]"				  		        		                  
	    fi
    fi
fi


if [  ${denoising} -eq 1 ]
then
    denoising_param="--denoising for_final"
else
    denoising_param="--denoising for_reg"
fi


if [  ${isPrisma} -eq 0 ]
then
    denoising_param="${denoising_param} --denoising_kernel_size 7"
fi

gibbs_param="--gibbs ${gibbs}"
rot_eddy_center_param="--rot_eddy_center center_slice "
s2v_param="--s2v 1"
repol_param="--repol 1 --outlier_replacement_mode aggressive"
if [ ${study} == "abcd" ]
then
    drift_param="--drift linear"
else
    drift_param=""
fi    



if [ -z "${output_dir}"  ]  
then
   final_folder=`dirname ${up_nifti}`/final_data
else
   final_folder=${output_dir}/sub-${subj_id}/ses-${ses_id}/dwi
   mkdir -p ${output_dir}/sub-${subj_id}/ses-${ses_id}/dwi
fi
mkdir -p ${final_folder}


if [ ${study} == "abcd" ]
then
    outputname=sub-${subj_id}_ses-${ses_id}_dwi.nii
    output_param="-o ${proc_dir}/${outputname}"
else
    outputname=sub-${subj_id}_ses-${ses_id}_space-ACPC_desc-preproc_dwi.nii
    output_param="-o ${proc_dir}/${outputname}"
fi






if [ ! -z "${output_res}"  ]   
then
    if [[ "${output_res}" == *"x"*  ]]
    then
        res=(${output_res//x/ })
        output_res_param="--output_res ${res[0]} ${res[1]} ${res[2]}"     
    else
        output_res_param="--output_res ${output_res} ${output_res} ${output_res}"    
    fi      
fi
if [  -z "${output_size}"  ]   
then
    if [ ${study} == "hbcd" ]
    then
        output_nvoxels_param="--output_voxels 74 89 72"      
    else
        output_nvoxels_param="--output_voxels 140 140 140"      
    fi

    output_nvoxels_param=""
else
    if [[ "${output_size}" == *"x"*  ]]
    then
        size=(${output_size//x/ })
        output_nvoxels_param="--output_voxels ${size[0]} ${size[1]} ${size[2]}"      
    fi      
fi

if [ -z ${output_orientation}  ]
then
    output_orientation_param=""
else
    output_orientation_param="--output_orientation ${output_orientation}"
fi
# the B0dx is signed along the output's +j, so j must run anterior-posterior (TORTOISE's
# default is LPS)
case "${output_orientation:-LPS}" in
    ?[AP]?) ;;
    *) fail "--output_orientation=${output_orientation}: its second letter must be A or P" ;;
esac
output_gradnonlin_Bmtxt_type_param="--output_gradnonlin_Bmtxt_type grad_dev"
   
if [ ${gpu} -eq 1 ]
then
    TORTOISE_exec="TORTOISEProcess_cuda"
else
    TORTOISE_exec="TORTOISEProcess"
fi

if [ ${study} == "abcd" ]
then
    IsDirectionPositiveDet ${up_nifti}
    det_rc=$?
    [ ${det_rc} -le 1 ] || fail "IsDirectionPositiveDet failed (exit ${det_rc}) on ${up_nifti}"
    if [ ${det_rc} -eq 1 ]
    then 
        flip_param="--flipY 1"
    else    
        flip_param=""
    fi
fi

if [ ${study} == "abcd" ]
then
    data_combination_param="--output_data_combination JacSep"
else
    data_combination_param="--output_data_combination JacConcat"
fi    

# down bvals/bvecs are found by basename: staging names everything down.*
my_RUN "${TORTOISE_exec}  \
-u ${up_nifti}  -d ${down_nifti} \
-t ${proc_dir} -s ${structural_for_TORTOISE} \
${gradnonlin_param}  ${denoising_param} \
${gibbs_param} ${rot_eddy_center_param} \
${s2v_param} ${repol_param} ${drift_param} \
${output_param} ${output_res_param} ${output_nvoxels_param} \
${output_orientation_param} ${output_gradnonlin_Bmtxt_type_param} \
${data_combination_param} ${flip_param} \
--DRBUDDI_rigid_for_heavily_distorted 1"



# structural_0 is the SynB0 synthetic T2 when the T1w branch ran, else the real T2w
if [ ${useT1} -eq 1 ]; then
    t2w_name=sub-${subj_id}_ses-${ses_id}_res-dwi_desc-synth_T2w.nii
else
    t2w_name=${t2_stem}_res-dwi_T2w.nii
fi
cp ${proc_dir}/structural_0.nii ${proc_dir}/${t2w_name} || fail "copying the TORTOISE structural_0.nii"
gzip -nf ${proc_dir}/${t2w_name} || fail "gzip ${proc_dir}/${t2w_name}"



# If Philips, split the merged series back into its input runs, published in dwi/;
# the merged series stays in proc_dwi/
if [ ! -z "${up_nifti2}"  ]
then
    for r in 1 2; do
        if [ ${r} -eq 1 ]; then run=${up_run}; v0=0; v1=$((up1_nvol - 1))
        else run=${up2_run}; v0=${up1_nvol}; v1=$((up1_nvol + up2_nvol - 1)); fi
        out=${final_folder}/sub-${subj_id}_ses-${ses_id}_${run}_desc-preproc_dwi
        my_RUN "ExtractDWISubset -i ${proc_dir}/${outputname} -o ${out}.nii -v ${v0}-${v1}"
        my_RUN "TORTOISEBmatrixToFSLBVecs ${out}.bmtxt"
        mv ${out}.bvals ${out}.bval && mv ${out}.bvecs ${out}.bvec && rm ${out}.bmtxt \
            || fail "renaming the ${run} gradient files"
        gzip -nf ${out}.nii || fail "gzip ${out}.nii"
    done
fi

gzip -nf ${proc_dir}/${outputname} || fail "gzip ${proc_dir}/${outputname}"

# PUBLISH: a single series is copied from proc_dwi/ (the Philips runs are already in dwi/)
proc_stem=${proc_dir}/`basename ${outputname} .nii`
if [ -z "${up_nifti2}" ]; then
    out_stem=sub-${subj_id}_ses-${ses_id}_${up_run}_desc-preproc_dwi
    cp ${proc_stem}.nii.gz ${final_folder}/${out_stem}.nii.gz || fail "copying the output DWI to ${final_folder}"
    cp ${proc_stem}.bvals ${final_folder}/${out_stem}.bval || fail "copying the output bvals to ${final_folder}"
    cp ${proc_stem}.bvecs ${final_folder}/${out_stem}.bvec || fail "copying the output bvecs to ${final_folder}"
fi

# Per-volume TORTOISE files of the processed series; each published series takes its rows.
# moteddy_transformations.txt has one bracketed list per volume; layout as documented for
# TORTOISE v1.0.5 _rpd.transformations (tortoisedti.nichd.nih.gov/transformation-files.html):
# elements 1-3 translations x,y,z (mm), 4-6 rotations x,y,z (radians), 7-24 eddy terms
tortoise_stem=${proc_dir}/`basename ${up_nifti} .${up_ext}`_proc
n_vol=`fslinfo ${proc_stem}.nii.gz | awk '$1=="dim4"{print $2}'`
for t in moteddy_transformations outlier_map slice_resids; do
    n_rows=`grep -c . ${tortoise_stem}_${t}.txt` || fail "reading ${tortoise_stem}_${t}.txt"
    [ "${n_rows}" == "${n_vol}" ] || fail "${tortoise_stem}_${t}.txt has ${n_rows} lines, the output DWI ${n_vol} volumes"
done

# SpatialReference of the outputs: the sMRI image whose space they are in; bids:: when it is in
# this output tree (the sMRI stream writes there), bids:rawdata: when in the input root (-T1/-T2)
[ ${useT1} -eq 1 ] && ref_smri=${T1W_smri} || ref_smri=${T2W_smri}
ref_abs=$( realpath `dirname ${ref_smri}` )/`basename ${ref_smri}`
out_abs=$( realpath ${output_dir} )
raw_abs=$( realpath ${prf} )
if [[ ${ref_abs} == "${out_abs}"/* ]]; then
    spatial_ref=bids::${ref_abs#"${out_abs}"/}
elif [[ ${ref_abs} == "${raw_abs}"/* ]]; then
    spatial_ref=bids:rawdata:${ref_abs#"${raw_abs}"/}
else
    echo "WARNING: ${ref_smri} is in neither ${output_dir} nor ${prf}; SpatialReference is its file URI"
    spatial_ref=file://${ref_abs}
fi

# MOTION TSV (_desc-motion_timeseries.tsv) and sidecar of one published series, from rows first..first+nvol-1 (1-based).
# Its outlier map and slice residuals are TORTOISE's own for a single series, else written
# from those rows as proc_dwi/sub-X_ses-Y_run-NN_dwi_proc_{outlier_map,slice_resids}.txt.
# The motion reader wants degrees (mmil_mean_motion.m)
write_series_files ()
{
    local run=$1 in_json=$2 source=$3 first=$4 nvol=$5
    shift 5
    local sources="$*"
    local out=${final_folder}/sub-${subj_id}_ses-${ses_id}_${run}
    local rows=${tortoise_stem} t
    local last=$((first + nvol - 1))
    local n=`fslinfo ${out}_desc-preproc_dwi.nii.gz | awk '$1=="dim4"{print $2}'`
    [ "${n}" == "${nvol}" ] || fail "${out}_desc-preproc_dwi.nii.gz has ${n} volumes, expected ${nvol}"
    sed -n "${first},${last}p" ${tortoise_stem}_moteddy_transformations.txt | tr -d '[],' | awk -v OFS='\t' '
        BEGIN { print "t_indx", "rot_z", "rot_x", "rot_y", "trans_z", "trans_x", "trans_y"; deg = 45 / atan2(1, 1) }
        { print NR - 1, $6 * deg, $4 * deg, $5 * deg, $3, $1, $2 }' > ${out}_desc-motion_timeseries.tsv \
        || fail "writing ${out}_desc-motion_timeseries.tsv"
    # TR into the NIfTI header, b-values below 5 as 0 and the rest whole, bvecs as BIDS reads them
    local rep_time=`jq -r '.RepetitionTime // empty' ${in_json}`
    [ -n "${rep_time}" ] || fail "no RepetitionTime in ${in_json}"
    my_RUN "python3 ${py_scripts_folder}/cmig_finish_bids_dwi_series.py --nifti ${out}_desc-preproc_dwi.nii.gz --tr ${rep_time} --bval ${out}_desc-preproc_dwi.bval --bvec ${out}_desc-preproc_dwi.bvec"
    if [ ${nvol} -ne ${n_vol} ]; then
        rows=${proc_dir}/sub-${subj_id}_ses-${ses_id}_${run}_dwi_proc
        for t in outlier_map slice_resids; do
            sed -n "${first},${last}p" ${tortoise_stem}_${t}.txt > ${rows}_${t}.txt || fail "writing ${rows}_${t}.txt"
        done
    fi
    my_RUN "python3 ${py_scripts_folder}/cmig_write_bids_dwi_sidecar.py --in_json ${in_json} --source ${source} --dwi ${out}_desc-preproc_dwi.nii.gz --gradwarp ${gw_applied} --outlier_map ${rows}_outlier_map.txt --slice_resids ${rows}_slice_resids.txt --bval ${out}_desc-preproc_dwi.bval --spatial_reference ${spatial_ref} ${shim_match:+--shim_match ${shim_match}} ${fmap_method:+--fmap_selection ${fmap_method}} --report ${proc_dir}/processing_report.json --kind dwi --sources ${sources} -o ${out}_desc-preproc_dwi.json"
}
[ -n "${gradnonlin_param}" ] && gw_applied=1 || gw_applied=0
# Sources of each series: its rawdata scan, then the other inputs it was corrected with
if [ -n "${up_nifti2}" ]; then
    write_series_files ${up_run} ${up_json} ${up_nifti1} 1 ${up1_nvol} ${up_raw_uri} ${up2_raw_uri} ${down_raw_uri}
    write_series_files ${up2_run} ${proc_dir}/${up2_base}.json ${up_nifti2} $((up1_nvol + 1)) ${up2_nvol} ${up2_raw_uri} ${up_raw_uri} ${down_raw_uri}
else
    write_series_files ${up_run} ${up_json} ${up_nifti1} 1 ${n_vol} ${up_raw_uri} ${down_raw_uri}
fi

if [ ! -z ${T1W_structural} ]
then
    my_RUN "ApplyTransformationToScalar2 ${T1W_structural} ${id_xfm} ${proc_dir}/${t1_stem}_res-dwi_T1w.nii.gz ${proc_dir}/${t2w_name}.gz  BSP"
fi

#finally generating an extra image for susceptibility fields in Hz
sidecar_vals=`python3 ${py_scripts_folder}/cmig_read_bids_sidecar.py ${up_json} EffectiveEchoSpacing AcquisitionMatrixPE SliceThickness` \
    || fail "reading ${up_json}"
eval "${sidecar_vals}"
echo_spacing=${EffectiveEchoSpacing}
mat_pe_size=${AcquisitionMatrixPE}
res=${SliceThickness}

#readout=`echo " $mat_pe_size * $echo_spacing" | bc`

#hz_to_mm_scale=`echo " $readout * $res" | bc`


my_RUN "DisplacementFieldOperations ${proc_dir}/deformation_temp.nii - ${proc_dir}/deformation_FINV.nii.gz ${proc_dir}/deformation_MINV.nii.gz"
my_RUN "DisplacementFieldOperations ${proc_dir}/deformation_temp.nii m ${proc_dir}/deformation_temp.nii 0.5"

my_RUN "ExtractComponentFromDisplacementField -i ${proc_dir}/deformation_temp.nii -o ${proc_dir}/deformation_temp_phase.nii -v 1"
# Hz = mm / (EffectiveEchoSpacing x PE field of view); PE is j, so the FOV is dim2 x pixdim2 of the DWI
pe_nvox=`fslinfo ${up_nifti} | awk '$1=="dim2"{print $2}'`
pe_vox_mm=`fslinfo ${up_nifti} | awk '$1=="pixdim2"{print $2}'`
[ -n "${pe_nvox}" ] && [ -n "${pe_vox_mm}" ] || fail "reading dim2/pixdim2 of ${up_nifti}"
hz_to_mm_scale=`echo "${EffectiveEchoSpacing} * ${pe_nvox} * ${pe_vox_mm}" | bc -l`
my_RUN "ImageMath 3 ${proc_dir}/deformation_field_Hz.nii.gz / ${proc_dir}/deformation_temp_phase.nii ${hz_to_mm_scale}"

# B0DX: mm map moved into the output space for the QC stream's b0_rms.
# The map is on TORTOISE's native-DWI grid; b0_to_str_rigidtrans.hdf5 takes it to the
# structural space, and it is resampled onto the output grid (structural_0)
fmap_folder=`dirname ${final_folder}`/fmap
mkdir -p ${fmap_folder} || fail "creating ${fmap_folder}"
b0dx=${fmap_folder}/sub-${subj_id}_ses-${ses_id}_acq-dwi_desc-AP${down_run#run-}PA${up_run#run-}_B0dx.nii.gz
my_RUN "ApplyTransformationToScalar2 ${proc_dir}/deformation_temp_phase.nii ${proc_dir}/b0_to_str_rigidtrans.hdf5 ${b0dx} ${proc_dir}/structural_0.nii Lin"
# the map is the displacement towards posterior (TORTOISE's LPS y); sign it along the
# output's +j, which points anterior when the output orientation is ?A? (e.g. LAS)
if [ "${output_orientation:1:1}" == "A" ]; then
    my_RUN "ImageMath 3 ${b0dx} m ${b0dx} -1"
fi

# the published DWI series (bids:: URIs): the field map's IntendedFor and the B0dx's Sources
dwi_uris=(bids::sub-${subj_id}/ses-${ses_id}/dwi/sub-${subj_id}_ses-${ses_id}_${up_run}_desc-preproc_dwi.nii.gz)
[ -n "${up_nifti2}" ] && dwi_uris+=(bids::sub-${subj_id}/ses-${ses_id}/dwi/sub-${subj_id}_ses-${ses_id}_${up2_run}_desc-preproc_dwi.nii.gz)

# Corrected reverse-PE field map, on the DWI grid, for the QC stream's b0_rms (it handles
# acq-dwi only if processed fmap/*_acq-dwi_dir-*_desc-preproc_epi scans exist and masks from them);
# named as the rawdata scan with desc-preproc added, as cmig_proc_bids_func names its field maps
if [[ ${down_base} == *_epi ]]
then
    down_out=${down_base%_epi}_desc-preproc_epi
    gzip -nc ${proc_dir}/down_proc_TORTOISE_final.nii > ${fmap_folder}/${down_out}.nii.gz \
        || fail "writing ${fmap_folder}/${down_out}.nii.gz"
    rep_time=`jq -r '.RepetitionTime // empty' ${proc_dir}/down.json`
    [ -n "${rep_time}" ] || fail "no RepetitionTime in ${proc_dir}/down.json"
    my_RUN "python3 ${py_scripts_folder}/cmig_finish_bids_dwi_series.py --nifti ${fmap_folder}/${down_out}.nii.gz --tr ${rep_time}"
    # its sidecar as the DWI's: input sidecar fields (TotalReadoutTime, PhaseEncodingDirection, ...)
    # and the processed image's grid and processing flags; Sources its rawdata scan and the
    # DWI series it was corrected with, IntendedFor the published DWI series
    my_RUN "python3 ${py_scripts_folder}/cmig_write_bids_dwi_sidecar.py --in_json ${proc_dir}/down.json --source ${down_nifti} --dwi ${fmap_folder}/${down_out}.nii.gz --gradwarp ${gw_applied} --spatial_reference ${spatial_ref} --kind fmap --sources ${down_raw_uri} ${up_raw_uri} ${up2_raw_uri} --intended_for ${dwi_uris[*]} -o ${fmap_folder}/${down_out}.json"
fi

# B0DX SIDECAR: Sources are the published images the field was estimated from (the corrected
# reverse-PE field map, if published, and the DWI series), in the order of the name's AP, PA
b0dx_src=()
[[ ${down_base} == *_epi ]] && b0dx_src+=(bids::sub-${subj_id}/ses-${ses_id}/fmap/${down_out}.nii.gz)
b0dx_src+=("${dwi_uris[@]}")
b0dx_desc="B0 displacement field estimated by DRBUDDI (TORTOISEProcess) from the DWI series and the reverse phase-encoding field map, as half the difference of their two displacement fields, in mm along the phase-encode axis (the output's j axis)"
printf '%s\n' "${b0dx_src[@]}" | jq -R . | jq -s --arg d "${b0dx_desc}" --arg r "${spatial_ref}" \
    '{Description: $d, Sources: ., SpatialReference: $r}' > ${b0dx%.nii.gz}.json \
    || fail "writing ${b0dx%.nii.gz}.json"


# COMPRESS: gzip (without name and time in the header, -n) every uncompressed NIfTI left in proc_dwi/ and the output dwi/ folder.
# Last step: TORTOISE and the steps above read the .nii files, nothing after this does.
echo "Compressing NIfTI files..."
# sort -u: proc_dir is inside final_folder unless --proc_dir puts it elsewhere
find ${proc_dir} ${final_folder} -type f -name "*.nii" -print0 | sort -zu \
    | xargs -0 -r -n 1 -P ${OMP_NUM_THREADS} gzip -nf \
    || fail "compressing the NIfTI files in ${proc_dir} and ${final_folder}"

# processing_report.json: paths relative to proc_dwi/, where it sits, and .nii -> .nii.gz
# after the compression above, which compresses every .nii in proc_dwi/
report=${proc_dir}/processing_report.json
sed -i -e "s#\"${proc_dir}/#\"#g" -e 's#\.nii"#.nii.gz"#g' ${report} || fail "rewriting ${report}"
! grep -q '": "/' ${report} || fail "absolute paths left in ${report}"
