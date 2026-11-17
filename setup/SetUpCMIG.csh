#!/bin/tcsh -f

################################################################################
#   SetUpCMIG.csh: environment for cmig_pipeline
#
#   Usage: source SetUpCMIG.csh [version]
#     (default: the highest version installed)
#
#   Sets CMIG_PIPELINE_DIR, CMIG_PIPELINE_VERSION, CMIG_PIPELINE_ATLASDIR,
#   GRADIENT_DISTORTION_FIELDS_DIR and SURFACE_MORPHOMETRY_TOOLBOX_DIR, and
#   puts cmig_pipeline's bin/, csh/ and python/ on PATH and python/ on
#   PYTHONPATH, from the installs ${CMIG_PARENT}/<repo>/<repo>.<version>.
#   The version is cmig_pipeline_dev's; the versions of the other three
#   installs are read from its setup/versions.txt.
#
#   Runs by itself or after SetUpMMPS.csh; after it, every variable and the
#   cmig_pipeline entries of PATH and PYTHONPATH are reset to these installs
#   (the entries keep their place, below MMPS's).
#
#   For MATLAB, startup_CMIG.m (in setup/) adds the toolbox, the gradient
#   maps and cmig_pipeline's matlab/ to the path.
#
#   Created:  2026-10-27 dhagler
################################################################################

setenv CMIG_PARENT ${PUBSW}/packages/MMPS

set cmig_ok = 1
if ($#argv >= 1) then
  set cmig_ver = $argv[1]
else
  # the highest cmig_pipeline_dev.<X.Y.Z> installed, compared as versions
  set cmig_ver = `ls ${CMIG_PARENT}/cmig_pipeline_dev |& sed -n 's/^cmig_pipeline_dev\.\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\)$/\1/p' | sort -t. -k1,1n -k2,2n -k3,3n | tail -1`
  if ("$cmig_ver" == "") then
    echo "SetUpCMIG: ERROR: no version installed in ${CMIG_PARENT}/cmig_pipeline_dev; environment unchanged"
    set cmig_ok = 0
  endif
endif

################################################################################
### the installs

set cmig_dir = ${CMIG_PARENT}/cmig_pipeline_dev/cmig_pipeline_dev.${cmig_ver}
set cmig_verfile = ${cmig_dir}/setup/versions.txt
if ($cmig_ok && ! -e $cmig_verfile) then
  echo "SetUpCMIG: ERROR: $cmig_verfile not found; environment unchanged"
  set cmig_ok = 0
endif
if ($cmig_ok) then
  foreach repo ( surface_morphometry_toolbox gradient_distortion_fields cmig_pipeline_atlases )
    set ver = `awk -v r=$repo '$1 == r {print $2}' $cmig_verfile`
    set dir = ${CMIG_PARENT}/${repo}/${repo}.${ver}
    if ("$ver" == "" || ! -d $dir) then
      echo "SetUpCMIG: ERROR: $repo version '$ver' not installed; environment unchanged"
      set cmig_ok = 0
    endif
    set cmig_$repo = $dir
  end
endif

################################################################################
### PATH and PYTHONPATH: replace earlier cmig_pipeline entries in place,
###   add the ones missing at the front

if ($cmig_ok) then
  if ($?CMIG_PIPELINE_DIR) then
    set cmig_old = $CMIG_PIPELINE_DIR
    set path=`echo $path | tr ' ' '\n' | awk -v o="$cmig_old" -v n="$cmig_dir" 'index($0,o) == 1 {$0 = n substr($0,length(o)+1)} {print}' | paste -sd' '`
    if ($?PYTHONPATH) then
      setenv PYTHONPATH `echo $PYTHONPATH | tr ':' '\n' | awk -v o="$cmig_old" -v n="$cmig_dir" 'index($0,o) == 1 {$0 = n substr($0,length(o)+1)} {print}' | paste -sd':'`
    endif
    unset cmig_old
  endif
  foreach dir ( ${cmig_dir}/python ${cmig_dir}/csh ${cmig_dir}/bin )
    if (-e $dir && "$path" !~ *$dir*) set path = ($dir $path)
  end
  if (! $?PYTHONPATH) then
    setenv PYTHONPATH ${cmig_dir}/python
  else if ("$PYTHONPATH" !~ *${cmig_dir}/python*) then
    setenv PYTHONPATH ${cmig_dir}/python:${PYTHONPATH}
  endif

  setenv CMIG_PIPELINE_DIR $cmig_dir
  setenv CMIGVER $cmig_ver
  setenv CMIG_PIPELINE_VERSION $cmig_ver
  setenv CMIG_PIPELINE_ATLASDIR $cmig_cmig_pipeline_atlases
  setenv GRADIENT_DISTORTION_FIELDS_DIR $cmig_gradient_distortion_fields
  setenv SURFACE_MORPHOMETRY_TOOLBOX_DIR $cmig_surface_morphometry_toolbox
endif

unset cmig_ver cmig_ok cmig_dir cmig_verfile repo ver dir
unset cmig_surface_morphometry_toolbox cmig_gradient_distortion_fields cmig_cmig_pipeline_atlases
