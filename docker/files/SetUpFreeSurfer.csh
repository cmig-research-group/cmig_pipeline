#!/bin/tcsh -f
# Container version of the host's $PUBSH/bin/SetUpFreeSurfer.csh: the image
#   holds one FreeSurfer, so any requested version gets it. MATLAB sources
#   this and reads FREESURFER_HOME from its output (mmil_get_fshome), so it
#   must print nothing.

setenv FREESURFER_VER `basename $FREESURFER_HOME | sed 's/^Rocky8-R//'`
setenv FREESURFER_HOME ${PUBSW}/packages/freesurfer/Rocky8-R${FREESURFER_VER}
if (! $?SUBJECTS_DIR) setenv SUBJECTS_DIR $FREESURFER_HOME/subjects
setenv NO_FSFAST 1
setenv NO_FSL 1
setenv FS_FREESURFERENV_NO_OUTPUT 1
setenv FS_OVERRIDE 1
# as FreeSurfer's own SetUpFreeSurfer.csh, for a single install tree;
#   FreeSurferEnv.csh does not set it, and mri_mcadura_seg reads its model
#   from it
unsetenv FREESURFER_FSPYTHON
setenv FREESURFER_HOME_FSPYTHON $FREESURFER_HOME
source $FREESURFER_HOME/FreeSurferEnv.csh
if ("$path" !~ *$FREESURFER_HOME/bin*) set path = ($FREESURFER_HOME/bin $path)
