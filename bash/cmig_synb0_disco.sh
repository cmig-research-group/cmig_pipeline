#!/bin/bash

# cmig_synb0_disco.sh: Synb0-DISCO for the dwi step: registers the b=0 and
#   T1w to the MNI atlas (prepare_input.sh), averages the five folds of
#   Synb0-DISCO's inference into a synthetic undistorted b=0, and moves it
#   back to the b=0's space; run by cmig_proc_dwi.bash
#
# usage: cmig_synb0_disco.sh DATA_DIR [-i|--notopup] [-s|--stripped]
#   Synb0-DISCO from SYNB0_DISCO_DIR (set in the CMIG BIDS app image)
#
# Created:  11/02/2026 by Don Hagler
# Last mod: 11/02/2026 by Don Hagler
#
# Based on Synb0-DISCO's src/pipeline.sh (MASILab, Vanderbilt University;
#   MIT License), as pipelineOK.sh in the dwi code (abcd_preproc 81210b9)
data_dir=$1

SYNB0_dir=${SYNB0_DISCO_DIR}

TOPUP=0
MNI_T1_1_MM_FILE=${SYNB0_dir}/atlases/mni_icbm152_t1_tal_nlin_asym_09c.nii.gz

for arg in "$@"
do
    case $arg in
        -i|--notopup)
            TOPUP=0
	        ;;
    	-s|--stripped)
	        MNI_T1_1_MM_FILE=${SYNB0_dir}/atlases/mni_icbm152_t1_tal_nlin_asym_09c_mask.nii.gz
            ;;
    esac
done

# Set path for executable
#export PATH=$PATH:/extra

# Set up freesurfer
#export FREESURFER_HOME=/extra/freesurfer
#source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Set up FSL
#. /extra/fsl/etc/fslconf/fsl.sh
#export PATH=$PATH:/extra/fsl/bin
#export FSLDIR=/extra/fsl

# Set up ANTS
#export ANTSPATH=/extra/ANTS/bin/ants/bin/
#export PATH=$PATH:$ANTSPATH:/extra/ANTS/ANTs/Scripts

# Prepare input
${SYNB0_dir}/data_processing/prepare_input.sh ${data_dir}/INPUTS/b0.nii.gz ${data_dir}/INPUTS/T1.nii.gz $MNI_T1_1_MM_FILE ${SYNB0_dir}/atlases/mni_icbm152_t1_tal_nlin_asym_09c_2_5.nii.gz ${data_dir}/OUTPUTS

# Run inference
NUM_FOLDS=5
for i in $(seq 1 $NUM_FOLDS);
  do echo Performing inference on FOLD: "$i"
  python3 ${SYNB0_dir}/src/inference.py ${data_dir}/OUTPUTS/T1_norm_lin_atlas_2_5.nii.gz ${data_dir}/OUTPUTS/b0_d_lin_atlas_2_5.nii.gz ${data_dir}/OUTPUTS/b0_u_lin_atlas_2_5_FOLD_"$i".nii.gz ${SYNB0_dir}/src/train_lin/num_fold_"$i"_total_folds_"$NUM_FOLDS"_seed_1_num_epochs_100_lr_0.0001_betas_\(0.9\,\ 0.999\)_weight_decay_1e-05_num_epoch_*.pth
done

# Take mean
echo Taking ensemble average
fslmerge -t ${data_dir}/OUTPUTS/b0_u_lin_atlas_2_5_merged.nii.gz ${data_dir}/OUTPUTS/b0_u_lin_atlas_2_5_FOLD_*.nii.gz
fslmaths ${data_dir}/OUTPUTS/b0_u_lin_atlas_2_5_merged.nii.gz -Tmean ${data_dir}/OUTPUTS/b0_u_lin_atlas_2_5.nii.gz

# Apply inverse xform to undistorted b0
echo Applying inverse xform to undistorted b0
antsApplyTransforms -d 3 -i ${data_dir}/OUTPUTS/b0_u_lin_atlas_2_5.nii.gz -r ${data_dir}/INPUTS/b0.nii.gz -n BSpline -t [${data_dir}/OUTPUTS/epi_reg_d_ANTS.txt,1] -t [${data_dir}/OUTPUTS/ANTS0GenericAffine.mat,1] -o ${data_dir}/OUTPUTS/b0_u.nii.gz

# Smooth image
echo Applying slight smoothing to distorted b0
fslmaths ${data_dir}/INPUTS/b0.nii.gz -s 1.15 ${data_dir}/OUTPUTS/b0_d_smooth.nii.gz



# Done
echo FINISHED!!!
