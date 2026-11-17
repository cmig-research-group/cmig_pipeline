# CMIG BIDS app: the cmig bids pipeline and its MATLAB modules, compiled as
# the single app cmig_run_bids_module, with the python wrappers in front,
# Python 3.12 with pinned packages, and FSL, AFNI, FreeSurfer, the MATLAB
# Runtime, the dwi tools (TORTOISE, ANTs, c3d_affine_tool, Synb0-DISCO)
# and the bids-validator at pinned versions.
#
# Stages: one per package (fsl, afni, mcr, freesurfer, ants, c3d, synb0,
# tortoise, bidsvalidator, fetch); base (the OS, the system packages, Python 3.12 and its
# packages, nothing under /usr/pubsw), which docker build --target base
# gives as the python image the hosts run the pipeline python in (the
# build script's --python_sing_flag); and the app image, FROM base.
#
# Everything sits under /usr/pubsw, as on the CMIG hosts, so the python
# defaults (appdir, mrtdir) and the MATLAB shell strings ($PUBSH, $PUBSW)
# work unchanged; the code finds its own files through CMIG_PIPELINE_DIR,
# CMIG_PIPELINE_ATLASDIR and GRADIENT_DISTORTION_FIELDS_DIR.
#
# Build with python/cmig_build_run_bids_module.py (--docker_flag 1), which
# assembles the context below from the repo tree; a plain docker build of
# docker/ does not work, the context holding more than docker/.
#
# Build context (assembled by the build script):
#   Dockerfile
#   files/SetUpFreeSurfer.csh   container version of the host setup script
#   files/prune_keep.sh,        each package stage prunes its tree to its
#   files/keep_<pkg>.txt        keep-list (fsl, afni, mcr, freesurfer), the
#                               files the pipeline uses; PRUNE=0 skips it;
#                               tortoise is always pruned to keep_tortoise.txt
#   files/requirements.txt      the python packages, pinned with sha256
#   files/synb0_disco.patch     the dwi step's edits to Synb0-DISCO
#   documentation/              LICENSE.txt and licenses/, to /documentation;
#                               the third-party license texts are added
#                               from the packages during the build
#   cmig_pipeline/              from the cmig_pipeline tree: python/,
#                               parms/, bin/ (the binaries the modules run:
#                               reg, applyUnwarpB0, epic, with RUNPATH
#                               $ORIGIN/../lib), lib/ (their libraries:
#                               gsl 1.6, Intel cce 10.1.013 and MKL
#                               10.0.2.018, libg2c, and licenses/)
#   local/                      the pieces not downloaded by tag (empty
#                               where a tag is pinned below):
#     app/cmig_run_bids_module/   the compiled app, with
#                                 cmig_run_bids_module_help.json
#     atlases/                    the atlas folders the modules read
#     gradient_distortion_fields/ the gradient nonlinearity distortion maps
#     tortoise/                   the TORTOISE programs and settings/ (the
#                                 build script's --tortoise_dir)
#
# Pinned downloads (option B): the release step writes the tag and sha256 of
# each piece into the ARG defaults below, so a build from a checked-out
# release gets the same files with no arguments; an empty tag takes that
# piece from local/ instead (a build from cmig_pipeline_dev between
# releases):
#   CMIG_PIPELINE_TAG, APP_SHA256   cmig_pipeline release asset
#                                   cmig_run_bids_module_<tag>.tar.gz
#   ATLAS_TAG, ATLAS_SHA256         cmig_pipeline release asset
#                                   cmig_pipeline_atlases_<tag>.tar.gz
#   GDF_TAG, GDF_SHA256             the gradient_distortion_fields tag's
#                                   source archive (its .mgh maps)
#   TORTOISE_TAG, TORTOISE_ASSET,   TORTOISEV4 release asset TORTOISE_ASSET
#   TORTOISE_SHA256                 of release TORTOISE_TAG
#
# Run (a BIDS app):
#   docker run --rm -v /data/rawdata:/data/rawdata:ro \
#     -v /data/derivatives:/data/derivatives \
#     -v /path/license.txt:/license.txt:ro \
#     cmig_bids /data/rawdata /data/derivatives participant \
#       [--participant_label L ...] [--session_label S ...] [--module NAME] \
#       [options] --fs_license /license.txt
# participant runs the rawqc and full pipelines for each session (or the
# module or pipeline NAME alone), group the tabulated data and the
# concatenation across sessions; run.py --help gives the options and
# names, run.py --module NAME --help the options of NAME. The FreeSurfer
# license key is supplied at run time, never built in: bound in and named by
# --fs_license, which run.py checks and sets as FS_LICENSE (an FS_LICENSE
# set with -e also works).

ARG FSL_VERSION=6.0.4-centos7
ARG FSL_TARBALL=fsl-6.0.4-centos7_64.tar.gz
ARG FSL_SHA256=99d4a93c7ce01c355d919b5f80e3d1ffe7e50ccec6ff0bb10cfce4d1255e7916
ARG AFNI_VERSION=AFNI_26.2.06
# the commit the tag points to; a tag can be moved, a commit cannot
ARG AFNI_COMMIT=d202535f6bd424be64e14a2dd1db52ba831e7d76
# the GitHub archive of that commit; GitHub does not promise these bytes
#   stay the same, so a mismatch here first means re-checking the commit
ARG AFNI_SHA256=940ba6d7a47e34b0b6df30da25d0d23dc6b4c5ec513c0ce82e1d389c2910559c
ARG FS_VERSION=8.2.0
ARG FS_VER=820
ARG FS_RPM=freesurfer-Rocky8-8.2.0-1.x86_64.rpm
# recorded 10/13/2026 (10,173,878,856 bytes); none is published
ARG FS_SHA256=4e8bf45662ded87a37450b8a54541d74cfaa9e2f1312534ed7a63afc9eba5b95
# SynthStrip's infant model (not in the RPM), for the dwi step's HBCD sessions
ARG SYNTHSTRIP_INFANT_URL=https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/requirements/synthstrip.infant.1.pt
ARG SYNTHSTRIP_INFANT_SHA256=3fbedf345d37b34126d25f289fdee57a487a1f9107e30ef6a951544cc9b215ab
ARG MCR_VERSION="R2023a Update 4"
ARG MCR_ZIP=MATLAB_Runtime_R2023a_Update_4_glnxa64.zip
# recorded 10/13/2026 (4,596,033,993 bytes); none is published
ARG MCR_SHA256=124d3ccf0fb11f0275bcb76d7957891b70e9edc211cb91b3c000d9a224962403
# the dwi tools: ANTs, the 7 programs the dwi step runs, from the release's
#   Alma Linux 8 build; its license from the source tag
ARG ANTS_VERSION=2.6.5
ARG ANTS_ZIP=ants-2.6.5-almalinux8-X64-gcc.zip
ARG ANTS_SHA256=cea0617ddbdfc4d926f1ce5dc171df5b2627b125d50aa32d0c6ad402b2ec911c
ARG ANTS_LICENSE_SHA256=c6596eb7be8581c18be736c846fb9173b69eccf6ef94c5135893ec56bd92ba08
# Convert3D's c3d_affine_tool, from the 1.1.0 Linux binary package (the
#   newer packages need glibc 2.29 or later); its license (GPL-3.0,
#   unchanged since the commit below) from the source
ARG C3D_VERSION=1.1.0
ARG C3D_TARBALL=c3d-1.1.0-Linux-x86_64.tar.gz
ARG C3D_SHA256=c1afdd57af7892a4f3a1674969e61fc1976d4b3073c38d10fd346681fd0c1e41
ARG C3D_LICENSE_COMMIT=93acef7192f1f5688313dcc8b7e760b229291bff
ARG C3D_LICENSE_SHA256=8ceb4b9ee5adedde47b31e975c1d90c73ad27b6b165a1dcd80c7c545eb65b903
# Synb0-DISCO at a commit (v3.1-5-gd813f94; the model weights are in the
#   repo); the GitHub archive of that commit, as for AFNI
ARG SYNB0_VERSION=v3.1-5-gd813f94
ARG SYNB0_COMMIT=d813f9435cffb5b663f15b52d55ed0a27181f0b2
ARG SYNB0_SHA256=94c4eaa766cf0ca1de1ec188dc1e9a3b8852f0a1aa4c0b0709d027fdce6ffa0c
# the bids-validator (run.py --module bids_validator): the deno-compiled
#   binary of the schema validator, from the PyPI wheel bids-validator-deno
#   (the GitHub releases have no assets); 2.4.1 is the last release with the
#   binary in the wheel (3.x wheels need the deno package) and has the BIDS
#   1.11.1 schema; the binary (157,251,383 bytes) embeds the Deno 2.6.10
#   runtime, whose license (MIT) is added from its source tag
ARG BIDS_VALIDATOR_VERSION=2.4.1
ARG BIDS_VALIDATOR_WHEEL=bids_validator_deno-2.4.1-py2.py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl
ARG BIDS_VALIDATOR_URL=https://files.pythonhosted.org/packages/28/15/4c2b1758ab7a628d9f1f401e6883ead26087a15940705a575e37450a0637/bids_validator_deno-2.4.1-py2.py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl
ARG BIDS_VALIDATOR_SHA256=5087fcd13341f529b5404792a74e22f570e50f17e5ea4cabe9c91b1caf74cbf4
ARG DENO_VERSION=2.6.10
ARG DENO_LICENSE_SHA256=f62497fffecc0852960c8d3e6934b9db86d16396e9b604072e923892cae3a588
# Python 3.12 from Rocky 8 AppStream, at a fixed build
ARG PY312_RPM_VERSION=3.12.15-2.el8_10
# the base of every stage, pinned by digest (the multi-platform index of
#   rockylinux/rockylinux:8.10, Rocky's own repository on Docker Hub, on
#   10/29/2026; the official rockylinux image stops at 8.9)
ARG BASE_IMAGE=rockylinux/rockylinux:8.10@sha256:e8a49c5403b687db05d4d67333fa45808fbe74f36e683cec7abb1f7d0f2338c6
# git describe of the cmig_pipeline tree, set by the build script
# the release step writes the release tag here, as it does the pins below
ARG CMIG_PIPELINE_VERSION=unknown
# the image's name:tag, set by the build script; ENV CMIG_IMAGE_TAG, the
#   status files' app_version
ARG CMIG_IMAGE_TAG=
ARG GITHUB_URL=https://github.com/cmig-research-group
ARG CMIG_PIPELINE_TAG=1.0.0
ARG APP_SHA256=bf37be8250ec89296543e7c39ab37af8176361df322db9c059edbdb8b578facd
ARG ATLAS_TAG=1.0.0
ARG ATLAS_SHA256=e3511fc2a0e427d7a6158319273fa18d5cdfcf656292c89606865990b116557a
ARG GDF_TAG=1.0.0
ARG GDF_SHA256=acf5c81cc0cb936fec14f76571ffb6e546694a396b20d2abc1160b96b16eabcc
# the TORTOISE release (asset and sha256); an empty tag takes the programs
#   from local/tortoise
ARG TORTOISE_URL=https://github.com/QMICodeBase/TORTOISEV4
ARG TORTOISE_TAG=TORTOISEV4.3.0
ARG TORTOISE_ASSET=TORTOISEV4.3.0.tar.gz
ARG TORTOISE_SHA256=745cfcd4f78c42c384569dde66af40e098707cbe9acb7cb65f3df1043aa73358
# 1: prune FSL, AFNI, the MATLAB Runtime and FreeSurfer to their keep-lists
#   (files/keep_<pkg>.txt); 0: keep them whole (the build script's
#   --prune_flag)
ARG PRUNE=1

########################################################################
### FSL: the vendor tarball, verified, unpacked and renamed
########################################################################
FROM ${BASE_IMAGE} AS fsl
ARG FSL_VERSION
ARG FSL_TARBALL
ARG FSL_SHA256
ARG PRUNE
RUN dnf install -y curl tar gzip findutils && dnf clean all
WORKDIR /fetch
# the archive unpacks to fsl/, not fsl-<version>/
RUN curl -fsSL -o ${FSL_TARBALL} https://fsl.fmrib.ox.ac.uk/fsldownloads/${FSL_TARBALL} \
 && echo "${FSL_SHA256}  ${FSL_TARBALL}" | sha256sum -c - \
 && mkdir -p /usr/pubsw/packages/fsl \
 && tar -xzf ${FSL_TARBALL} -C /usr/pubsw/packages/fsl \
 && mv /usr/pubsw/packages/fsl/fsl /usr/pubsw/packages/fsl/fsl-${FSL_VERSION} \
 && rm -f ${FSL_TARBALL} \
 && mkdir -p /licenses/FSL \
 && cp /usr/pubsw/packages/fsl/fsl-${FSL_VERSION}/LICENCE /licenses/FSL/
# pruned in a step of its own, so a keep-list change does not repeat the
#   download
COPY files/prune_keep.sh files/keep_fsl.txt /prune/
RUN if [ "${PRUNE}" = 1 ]; then \
      bash /prune/prune_keep.sh /usr/pubsw/packages/fsl/fsl-${FSL_VERSION} /prune/keep_fsl.txt ; \
    fi
# fsl_abspath (called by epi_reg with --wmseg or existing _wmseg/_wmedge
#   files) is a fslpython script; replaced so fslpython can be pruned
COPY files/fsl_abspath /usr/pubsw/packages/fsl/fsl-${FSL_VERSION}/bin/fsl_abspath

########################################################################
### AFNI: built from the pinned source tag
########################################################################
FROM ${BASE_IMAGE} AS afni
ARG AFNI_VERSION
ARG AFNI_COMMIT
ARG AFNI_SHA256
ARG PRUNE
# AFNI's rocky_8 prerequisites without its R part, plus gcc-c++ and cmake,
# which dcm2niix_afni needs
RUN dnf install -y \
      findutils tcsh gcc gcc-c++ motif-devel gsl-devel libXpm-devel \
      mesa-libGLw-devel mesa-libGLU-devel libXi-devel glib2-devel \
      zlib-devel expat-devel libXt-devel libXmu-devel \
      make cmake bc git curl tar gzip \
 && dnf clean all
WORKDIR /build
RUN curl -fsSL -o /tmp/afni.tar.gz https://github.com/afni/afni/archive/${AFNI_COMMIT}.tar.gz \
 && if [ -n "${AFNI_SHA256}" ]; then echo "${AFNI_SHA256}  /tmp/afni.tar.gz" | sha256sum -c - ; \
    else echo "AFNI sha256 not pinned:" && sha256sum /tmp/afni.tar.gz ; fi \
 && tar xzf /tmp/afni.tar.gz --strip-components=1 \
 && rm -f /tmp/afni.tar.gz
RUN cd src \
 && cp Makefile.linux_rocky_8 Makefile \
 && /bin/bash -c 'set -o pipefail; make itall 2>&1 | tail -60' \
 && test -x linux_rocky_8/3dvolreg \
 && test -e linux_rocky_8/libmri.so
RUN mkdir -p /licenses/AFNI \
 && cp LICENSE.txt doc/README/README.copyright /licenses/AFNI/
COPY files/prune_keep.sh files/keep_afni.txt /prune/
RUN if [ "${PRUNE}" = 1 ]; then \
      bash /prune/prune_keep.sh /build/src/linux_rocky_8 /prune/keep_afni.txt ; \
    fi

########################################################################
### MATLAB Runtime R2023a Update 4, at the host's path
########################################################################
FROM ${BASE_IMAGE} AS mcr
ARG MCR_ZIP
ARG MCR_SHA256
ARG PRUNE
RUN dnf install -y curl unzip findutils && dnf clean all
WORKDIR /fetch
RUN curl -fsSL -o ${MCR_ZIP} \
      https://ssd.mathworks.com/supportfiles/downloads/R2023a/Release/4/deployment_files/installer/complete/glnxa64/${MCR_ZIP} \
 && if [ -n "${MCR_SHA256}" ]; then echo "${MCR_SHA256}  ${MCR_ZIP}" | sha256sum -c - ; \
    else echo "MCR sha256 not pinned:" && sha256sum ${MCR_ZIP} ; fi \
 && unzip -q ${MCR_ZIP} -d installer \
 && ./installer/install -mode silent -agreeToLicense yes \
      -destinationFolder /usr/pubsw/packages/matlab \
 && test -d /usr/pubsw/packages/matlab/R2023a/runtime/glnxa64 \
 && mkdir -p /licenses/MATLAB_Runtime \
 && cp installer/matlabruntime_license_agreement.pdf /licenses/MATLAB_Runtime/ \
 && rm -rf ${MCR_ZIP} installer
COPY files/prune_keep.sh files/keep_mcr.txt /prune/
RUN if [ "${PRUNE}" = 1 ]; then \
      bash /prune/prune_keep.sh /usr/pubsw/packages/matlab/R2023a /prune/keep_mcr.txt ; \
    fi

########################################################################
### FreeSurfer: the Rocky 8 RPM, installed and linked at the host's path
########################################################################
FROM ${BASE_IMAGE} AS freesurfer
ARG FS_VERSION
ARG FS_VER
ARG FS_RPM
ARG FS_SHA256
ARG SYNTHSTRIP_INFANT_URL
ARG SYNTHSTRIP_INFANT_SHA256
ARG PRUNE
RUN dnf install -y curl findutils && dnf clean all
RUN curl -fsSL -o /tmp/${FS_RPM} \
      https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/${FS_VERSION}/${FS_RPM} \
 && if [ -n "${FS_SHA256}" ]; then echo "${FS_SHA256}  /tmp/${FS_RPM}" | sha256sum -c - ; \
    else echo "FreeSurfer sha256 not pinned:" && sha256sum /tmp/${FS_RPM} ; fi \
 && dnf install -y /tmp/${FS_RPM} \
 && dnf clean all \
 && rm -f /tmp/${FS_RPM} \
 && fsdir=$(dirname $(rpm -qa | grep -i '^freesurfer' | head -1 | xargs rpm -ql \
      | grep '/FreeSurferEnv.csh$' | head -1)) \
 && test -d "$fsdir" \
 && mkdir -p /usr/pubsw/packages/freesurfer /licenses/FreeSurfer \
 && ln -s "$fsdir" /usr/pubsw/packages/freesurfer/Rocky8-R${FS_VER} \
 && cp "$fsdir"/docs/license*.txt /licenses/FreeSurfer/ \
 && curl -fsSL -o "$fsdir"/models/synthstrip.infant.1.pt ${SYNTHSTRIP_INFANT_URL} \
 && echo "${SYNTHSTRIP_INFANT_SHA256}  $fsdir/models/synthstrip.infant.1.pt" | sha256sum -c -
COPY files/prune_keep.sh files/keep_freesurfer.txt /prune/
RUN if [ "${PRUNE}" = 1 ]; then \
      bash /prune/prune_keep.sh $(readlink -f /usr/pubsw/packages/freesurfer/Rocky8-R${FS_VER}) \
        /prune/keep_freesurfer.txt ; \
    fi

########################################################################
### ANTs: the 7 programs the dwi step runs, from the release zip
########################################################################
FROM ${BASE_IMAGE} AS ants
ARG ANTS_VERSION
ARG ANTS_ZIP
ARG ANTS_SHA256
ARG ANTS_LICENSE_SHA256
RUN dnf install -y curl unzip findutils && dnf clean all
WORKDIR /fetch
RUN progs="antsApplyTransforms antsRegistration antsRegistrationSyNQuick.sh \
           CreateDisplacementField ImageMath N4BiasFieldCorrection PrintHeader" \
 && curl -fsSL -o ${ANTS_ZIP} \
      https://github.com/ANTsX/ANTs/releases/download/v${ANTS_VERSION}/${ANTS_ZIP} \
 && echo "${ANTS_SHA256}  ${ANTS_ZIP}" | sha256sum -c - \
 && mkdir -p /usr/pubsw/packages/ants/${ANTS_VERSION}/bin \
 && unzip -j ${ANTS_ZIP} $(for p in $progs ; do echo ants-${ANTS_VERSION}/bin/$p ; done) \
      -d /usr/pubsw/packages/ants/${ANTS_VERSION}/bin \
 && rm -f ${ANTS_ZIP} \
 && for p in $progs ; do test -x /usr/pubsw/packages/ants/${ANTS_VERSION}/bin/$p ; done \
 && mkdir -p /licenses/ANTs \
 && curl -fsSL -o /licenses/ANTs/COPYING.txt \
      https://raw.githubusercontent.com/ANTsX/ANTs/v${ANTS_VERSION}/COPYING.txt \
 && echo "${ANTS_LICENSE_SHA256}  /licenses/ANTs/COPYING.txt" | sha256sum -c -

########################################################################
### Convert3D: c3d_affine_tool only, from the binary package
########################################################################
FROM ${BASE_IMAGE} AS c3d
ARG C3D_VERSION
ARG C3D_TARBALL
ARG C3D_SHA256
ARG C3D_LICENSE_COMMIT
ARG C3D_LICENSE_SHA256
RUN dnf install -y curl tar gzip && dnf clean all
WORKDIR /fetch
# SourceForge redirects to a mirror; the sha256 fixes the bytes
RUN curl -fsSL -o ${C3D_TARBALL} \
      https://sourceforge.net/projects/c3d/files/c3d/Experimental/${C3D_TARBALL}/download \
 && echo "${C3D_SHA256}  ${C3D_TARBALL}" | sha256sum -c - \
 && mkdir -p /usr/pubsw/packages/c3d/${C3D_VERSION}/bin \
 && tar -xzf ${C3D_TARBALL} --strip-components=2 \
      -C /usr/pubsw/packages/c3d/${C3D_VERSION}/bin \
      ${C3D_TARBALL%.tar.gz}/bin/c3d_affine_tool \
 && rm -f ${C3D_TARBALL} \
 && test -x /usr/pubsw/packages/c3d/${C3D_VERSION}/bin/c3d_affine_tool \
 && mkdir -p /licenses/Convert3D \
 && curl -fsSL -o /licenses/Convert3D/COPYING.txt \
      https://raw.githubusercontent.com/pyushkevich/c3d/${C3D_LICENSE_COMMIT}/COPYING.txt \
 && echo "${C3D_LICENSE_SHA256}  /licenses/Convert3D/COPYING.txt" | sha256sum -c -

########################################################################
### Synb0-DISCO: the scripts, model weights and atlases the dwi step
###   uses, from the pinned commit, with the dwi step's edits
########################################################################
FROM ${BASE_IMAGE} AS synb0
ARG SYNB0_VERSION
ARG SYNB0_COMMIT
ARG SYNB0_SHA256
RUN dnf install -y curl tar gzip patch && dnf clean all
WORKDIR /build
# the edits (files/synb0_disco.patch): normalize_T1.sh runs mri_nu_correct.mni
#   with --ants-n4; prepare_input.sh calls normalize_T1.sh by its own
#   directory and takes the T1 brain mask beside the T1 (T1_brain.nii.gz) in
#   place of running bet; inference.py runs on the CPU
COPY files/synb0_disco.patch /tmp/synb0_disco.patch
RUN curl -fsSL -o /tmp/synb0.tar.gz \
      https://github.com/MASILab/Synb0-DISCO/archive/${SYNB0_COMMIT}.tar.gz \
 && echo "${SYNB0_SHA256}  /tmp/synb0.tar.gz" | sha256sum -c - \
 && tar -xzf /tmp/synb0.tar.gz --strip-components=1 \
 && rm -f /tmp/synb0.tar.gz \
 && patch -p1 --forward --batch < /tmp/synb0_disco.patch \
 && D=/usr/pubsw/packages/synb0-disco/${SYNB0_VERSION} \
 && mkdir -p $D/src/train_lin $D/data_processing $D/atlases \
 && cp -p src/inference.py src/model.py src/util.py $D/src/ \
 && cp -p src/train_lin/*.pth $D/src/train_lin/ \
 && test $(ls $D/src/train_lin/*.pth | wc -l) = 5 \
 && cp -p data_processing/prepare_input.sh data_processing/normalize_T1.sh $D/data_processing/ \
 && chmod 755 $D/data_processing/*.sh \
 && cp -p atlases/*.nii.gz atlases/cmd.txt $D/atlases/ \
 && cp -p README.md synb0-license-MIT.txt $D/ \
 && mkdir -p /licenses/Synb0-DISCO \
 && cp synb0-license-MIT.txt /licenses/Synb0-DISCO/

########################################################################
### TORTOISE: the programs the dwi step runs and settings/, from the
###   release asset, or from local/tortoise when TORTOISE_TAG is empty
########################################################################
FROM ${BASE_IMAGE} AS tortoise
ARG TORTOISE_URL
ARG TORTOISE_TAG
ARG TORTOISE_ASSET
ARG TORTOISE_SHA256
RUN dnf install -y curl tar gzip unzip findutils && dnf clean all
COPY local/tortoise /local/tortoise
COPY files/prune_keep.sh files/keep_tortoise.txt /prune/
# the package's top directory is the one with bin/TORTOISEProcess; pruned to
#   the keep-list whatever PRUNE is (the package holds ~100 programs)
RUN D=/usr/pubsw/packages/tortoise/${TORTOISE_TAG:-local} \
 && if [ -n "${TORTOISE_TAG}" ] ; then \
      test -n "${TORTOISE_ASSET}" -a -n "${TORTOISE_SHA256}" \
        || { echo "** TORTOISE: tag ${TORTOISE_TAG} has no asset or sha256"; exit 1; } ; \
      curl -fsSL -o /tmp/${TORTOISE_ASSET} \
        ${TORTOISE_URL}/releases/download/${TORTOISE_TAG}/${TORTOISE_ASSET} \
      && echo "${TORTOISE_SHA256}  /tmp/${TORTOISE_ASSET}" | sha256sum -c - \
      && mkdir /tmp/pkg \
      && case ${TORTOISE_ASSET} in \
           *.zip) unzip -q /tmp/${TORTOISE_ASSET} -d /tmp/pkg ;; \
           *) tar -xf /tmp/${TORTOISE_ASSET} -C /tmp/pkg ;; \
         esac \
      && rm -f /tmp/${TORTOISE_ASSET} \
      && top=$(dirname $(dirname $(find /tmp/pkg -path '*/bin/TORTOISEProcess' | head -1))) ; \
    else top=/local/tortoise ; fi \
 && test -x $top/bin/TORTOISEProcess \
 && mkdir -p /licenses/TORTOISE \
 && for f in $top/LICENSE* $top/COPYING* ; do \
      if [ -f "$f" ] ; then cp "$f" /licenses/TORTOISE/ ; fi ; done \
 && bash /prune/prune_keep.sh $top /prune/keep_tortoise.txt \
 && mkdir -p $(dirname $D) \
 && mv $top $D \
 && rm -rf /tmp/pkg /local

########################################################################
### bids-validator: the deno-compiled binary from the PyPI wheel, installed
###   as bids-validator
########################################################################
FROM ${BASE_IMAGE} AS bidsvalidator
ARG BIDS_VALIDATOR_VERSION
ARG BIDS_VALIDATOR_WHEEL
ARG BIDS_VALIDATOR_URL
ARG BIDS_VALIDATOR_SHA256
ARG DENO_VERSION
ARG DENO_LICENSE_SHA256
RUN dnf install -y curl unzip && dnf clean all
WORKDIR /fetch
# the wheel's binary is bids-validator-deno; renamed to the name run.py
#   calls; its license (MIT) from the wheel's dist-info
RUN W=bids_validator_deno-${BIDS_VALIDATOR_VERSION} \
 && D=/usr/pubsw/packages/bids-validator/${BIDS_VALIDATOR_VERSION} \
 && curl -fsSL -o ${BIDS_VALIDATOR_WHEEL} ${BIDS_VALIDATOR_URL} \
 && echo "${BIDS_VALIDATOR_SHA256}  ${BIDS_VALIDATOR_WHEEL}" | sha256sum -c - \
 && mkdir -p $D/bin /licenses/bids-validator \
 && unzip -j ${BIDS_VALIDATOR_WHEEL} $W.data/scripts/bids-validator-deno -d $D/bin \
 && mv $D/bin/bids-validator-deno $D/bin/bids-validator \
 && chmod 755 $D/bin/bids-validator \
 && unzip -j ${BIDS_VALIDATOR_WHEEL} $W.dist-info/licenses/LICENSE -d /licenses/bids-validator \
 && rm -f ${BIDS_VALIDATOR_WHEEL} \
 && curl -fsSL -o /licenses/bids-validator/LICENSE_deno.md \
      https://raw.githubusercontent.com/denoland/deno/v${DENO_VERSION}/LICENSE.md \
 && echo "${DENO_LICENSE_SHA256}  /licenses/bids-validator/LICENSE_deno.md" | sha256sum -c -

########################################################################
### the compiled app, the atlases and the gradient maps: downloaded by
###   tag and checked against their sha256, or taken from local/
########################################################################
FROM ${BASE_IMAGE} AS fetch
ARG GITHUB_URL
ARG CMIG_PIPELINE_TAG
ARG APP_SHA256
ARG ATLAS_TAG
ARG ATLAS_SHA256
ARG GDF_TAG
ARG GDF_SHA256
RUN dnf install -y curl tar gzip findutils && dnf clean all
COPY local /local
# get NAME TAG SHA256 URL: download URL to /tmp/NAME.tar.gz and check it,
#   failing the build on a missing sha256 or a mismatch
RUN mkdir -p /out/app /out/atlases /out/gradient_distortion_fields \
 && get() { test -n "$3" || { echo "** $1: tag $2 has no sha256"; exit 1; } ; \
      curl -fsSL -o /tmp/$1.tar.gz "$4" \
      && echo "$3  /tmp/$1.tar.gz" | sha256sum -c - ; } \
 && if [ -n "${CMIG_PIPELINE_TAG}" ] ; then \
      get app ${CMIG_PIPELINE_TAG} "${APP_SHA256}" \
        ${GITHUB_URL}/cmig_pipeline/releases/download/${CMIG_PIPELINE_TAG}/cmig_run_bids_module_${CMIG_PIPELINE_TAG}.tar.gz \
      && tar -xzf /tmp/app.tar.gz -C /out/app ; \
    else cp -a /local/app/. /out/app/ ; fi \
 && test -x /out/app/cmig_run_bids_module/run_cmig_run_bids_module.sh \
 && if [ -n "${ATLAS_TAG}" ] ; then \
      get atlases ${ATLAS_TAG} "${ATLAS_SHA256}" \
        ${GITHUB_URL}/cmig_pipeline/releases/download/${ATLAS_TAG}/cmig_pipeline_atlases_${ATLAS_TAG}.tar.gz \
      && tar -xzf /tmp/atlases.tar.gz -C /out/atlases ; \
    else cp -a /local/atlases/. /out/atlases/ ; fi \
 && for a in T1_Atlas T1T2 DTI_Atlas fuzzy_clusters ; do \
      test -d /out/atlases/$a || { echo "** atlas $a missing"; exit 1; } ; done \
 && if [ -n "${GDF_TAG}" ] ; then \
      get gdf ${GDF_TAG} "${GDF_SHA256}" \
        ${GITHUB_URL}/gradient_distortion_fields/archive/refs/tags/${GDF_TAG}.tar.gz \
      && mkdir /tmp/gdf && tar -xzf /tmp/gdf.tar.gz -C /tmp/gdf \
      && find /tmp/gdf -mindepth 2 -maxdepth 2 -name '*.mgh' \
           -exec cp -p {} /out/gradient_distortion_fields/ \; ; \
    else cp -a /local/gradient_distortion_fields/. /out/gradient_distortion_fields/ ; fi \
 && test -n "$(ls /out/gradient_distortion_fields/*.mgh)" \
 && rm -rf /tmp/*.tar.gz /tmp/gdf

########################################################################
### base: the OS, the system packages and Python 3.12 with its packages,
###   nothing under /usr/pubsw; the app image is built FROM it, and it is
###   the python image the hosts run the pipeline python in, with the host's
###   /usr/pubsw bound in
########################################################################
FROM ${BASE_IMAGE} AS base
ARG PY312_RPM_VERSION

# system packages: tcsh (mmil_unix runs every command with tcsh -c, and
#   recon-all runs under /bin/csh), the libraries of the compiled app and
#   cmig_pipeline binaries, the FSL and AFNI runtime libraries, python 3.12,
#   jq and time (the dwi step), and the FreeSurfer RPM's requirements (rpm
#   -qR, mapped to packages; the RPM itself is installed in its own stage)
RUN dnf update -y \
 && dnf install -y \
      tcsh bc tzdata which procps-ng findutils gzip tar unzip file \
      fontconfig glibc-langpack-en perl jq time \
      libXt libnsl libgomp libXp libXpm libXmu mesa-libGLU libquadmath \
      openblas libgfortran libstdc++ libgcc zlib \
      motif gsl mesa-libGLw libXi glib2 expat \
      python3.12-${PY312_RPM_VERSION} python3.12-libs-${PY312_RPM_VERSION} \
      python3.12-pip \
      bzip2-libs freetype krb5-libs libffi libglvnd-glx libjpeg-turbo \
      libpng12 libpng15 libICE libSM libtiff libX11 libxcb libxcrypt \
      libXext libXft libXScrnSaver perl-interpreter qt5-qtbase \
      qt5-qtbase-gui qt5-qtx11extras tcl xz-libs \
 && dnf clean all

# python3 and python are python3.12: the alternatives, and links in
#   /usr/local/bin, which the host wrapper puts first on PATH
RUN alternatives --set python3 /usr/bin/python3.12 \
 && alternatives --set python /usr/bin/python3.12 \
 && ln -s /usr/bin/python3.12 /usr/local/bin/python3 \
 && ln -s /usr/bin/python3.12 /usr/local/bin/python \
 && test "$(python3 -c 'import sys; print(sys.version_info[:2])')" = "(3, 12)" \
 && test "$(python -c 'import sys; print(sys.version_info[:2])')" = "(3, 12)"

# the python packages, every file pinned by sha256 (files/requirements.txt),
#   into the system python (/usr/local/lib*/python3.12/site-packages)
COPY files/requirements.txt /tmp/requirements.txt
RUN pip3.12 install --no-cache-dir --require-hashes -r /tmp/requirements.txt \
 && rm -f /tmp/requirements.txt \
 && python3 -c 'import numpy, scipy, pandas, nibabel, matplotlib, torch, torchvision'

# the license texts of the python packages, from each package's dist-info
#   as installed (LICENSE*, COPYING*, NOTICE*, AUTHORS* and licenses/)
RUN L=/documentation/licenses/python \
 && mkdir -p $L \
 && for d in /usr/local/lib/python3.12/site-packages/*.dist-info \
             /usr/local/lib64/python3.12/site-packages/*.dist-info ; do \
      test -d "$d" || continue ; \
      n=$(basename $d .dist-info) ; \
      mkdir -p $L/$n ; \
      ( cd $d && find . -type f \( -iname 'LICEN[CS]E*' -o -iname 'COPYING*' \
          -o -iname 'NOTICE*' -o -iname 'AUTHORS*' -o -path './licenses/*' \) -print0 \
        | xargs -0 -r cp --parents -t $L/$n ) ; \
      test -n "$(ls -A $L/$n 2>/dev/null)" \
        || { echo "** $n: no license file in its dist-info"; exit 1; } ; \
    done

ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    TZ=America/Los_Angeles

########################################################################
### the app image
########################################################################
FROM base
ARG BASE_IMAGE
ARG FSL_VERSION
ARG AFNI_VERSION
ARG AFNI_COMMIT
ARG FS_VERSION
ARG FS_VER
ARG FS_RPM
ARG FS_SHA256
ARG MCR_VERSION
ARG CMIG_PIPELINE_VERSION
ARG CMIG_PIPELINE_TAG
ARG APP_SHA256
ARG ATLAS_TAG
ARG ATLAS_SHA256
ARG GDF_TAG
ARG GDF_SHA256
ARG ANTS_VERSION
ARG C3D_VERSION
ARG SYNB0_VERSION
ARG SYNB0_COMMIT
ARG TORTOISE_TAG
ARG TORTOISE_ASSET
ARG TORTOISE_SHA256
ARG BIDS_VALIDATOR_VERSION
ARG DENO_VERSION

COPY --from=fsl /usr/pubsw/packages/fsl /usr/pubsw/packages/fsl
COPY --from=afni /build/src/linux_rocky_8 /usr/pubsw/packages/afni/${AFNI_VERSION}/linux_rocky_8
COPY --from=mcr /usr/pubsw/packages/matlab /usr/pubsw/packages/matlab
COPY --from=afni /licenses/AFNI /documentation/licenses/AFNI
COPY --from=mcr /licenses/MATLAB_Runtime /documentation/licenses/MATLAB_Runtime
COPY --from=freesurfer /usr/local/freesurfer /usr/local/freesurfer
COPY --from=freesurfer /usr/pubsw/packages/freesurfer /usr/pubsw/packages/freesurfer
COPY --from=freesurfer /licenses/FreeSurfer /documentation/licenses/FreeSurfer
COPY --from=fsl /licenses/FSL /documentation/licenses/FSL
COPY --from=ants /usr/pubsw/packages/ants /usr/pubsw/packages/ants
COPY --from=ants /licenses/ANTs /documentation/licenses/ANTs
COPY --from=c3d /usr/pubsw/packages/c3d /usr/pubsw/packages/c3d
COPY --from=c3d /licenses/Convert3D /documentation/licenses/Convert3D
COPY --from=synb0 /usr/pubsw/packages/synb0-disco /usr/pubsw/packages/synb0-disco
COPY --from=synb0 /licenses/Synb0-DISCO /documentation/licenses/Synb0-DISCO
COPY --from=tortoise /usr/pubsw/packages/tortoise /usr/pubsw/packages/tortoise
COPY --from=tortoise /licenses/TORTOISE /documentation/licenses/TORTOISE
COPY --from=bidsvalidator /usr/pubsw/packages/bids-validator /usr/pubsw/packages/bids-validator
COPY --from=bidsvalidator /licenses/bids-validator /documentation/licenses/bids-validator

# cmig_pipeline, its atlases and gradient maps, and the setup script MATLAB
#   sources for FreeSurfer
COPY files/SetUpFreeSurfer.csh /usr/pubsw/bin/SetUpFreeSurfer.csh
COPY cmig_pipeline /usr/pubsw/packages/cmig_pipeline
COPY --from=fetch /out/atlases /usr/pubsw/packages/cmig_pipeline_atlases
COPY --from=fetch /out/gradient_distortion_fields /usr/pubsw/packages/gradient_distortion_fields
COPY documentation /documentation

# the license texts of the third-party software, copied from the packages
#   as installed, so they match the versions in the image (the python
#   packages' in the base stage)
RUN L=/documentation/licenses \
 && cp -r /usr/pubsw/packages/cmig_pipeline/lib/licenses/GSL \
      /usr/pubsw/packages/cmig_pipeline/lib/licenses/Intel \
      /usr/pubsw/packages/cmig_pipeline/lib/licenses/libg2c $L/ \
 && rm $L/GSL/gsl-1.6.tar.gz \
 && test -s $L/FSL/LICENCE && test -s $L/GSL/COPYING \
 && test -s $L/ANTs/COPYING.txt && test -s $L/Convert3D/COPYING.txt \
 && test -s $L/Synb0-DISCO/synb0-license-MIT.txt && test -s $L/TORTOISE/LICENSE \
 && test -s $L/bids-validator/LICENSE && test -s $L/bids-validator/LICENSE_deno.md \
 && test -s $L/MNI_ICBM152/COPYING && test -n "$(ls $L/python)" \
 && test -s $L/PolyForm-Noncommercial-1.0.0.md \
 && ls $L/FreeSurfer/license.freesurfer_SLA.txt
COPY --from=fetch /out/app/cmig_run_bids_module /usr/pubsw/packages/MMPS/apps/cmig_run_bids_module

# the BIDS app front end: run.py is cmig_run_bids_app.py, whose directory
#   python then puts on sys.path
RUN chmod 755 /usr/pubsw/packages/cmig_pipeline/python/cmig_run_bids_app.py \
 && ln -s /usr/pubsw/packages/cmig_pipeline/python/cmig_run_bids_app.py /run.py

# the environment the host sets in its cshrc and the SetUp scripts;
#   the MATLAB Runtime is not on LD_LIBRARY_PATH: run_cmig_run_bids_module.sh
#   adds it for the app alone
ENV PUBSW=/usr/pubsw \
    PUBSH=/usr/pubsw \
    SCRATCHROOT=/tmp \
    MCR_CACHE_ROOT=/tmp \
    CMIG_DOCKER=1 \
    CMIG_PIPELINE_DIR=/usr/pubsw/packages/cmig_pipeline \
    CMIG_PIPELINE_VERSION=${CMIG_PIPELINE_VERSION} \
    CMIG_PIPELINE_ATLASDIR=/usr/pubsw/packages/cmig_pipeline_atlases \
    GRADIENT_DISTORTION_FIELDS_DIR=/usr/pubsw/packages/gradient_distortion_fields \
    FSLDIR=/usr/pubsw/packages/fsl/fsl-${FSL_VERSION} \
    FSLOUTPUTTYPE=NIFTI \
    FSLMULTIFILEQUIT=TRUE \
    AFNIDIR=/usr/pubsw/packages/afni/${AFNI_VERSION}/linux_rocky_8 \
    AFNI_PLUGINPATH=/usr/pubsw/packages/afni/${AFNI_VERSION}/linux_rocky_8 \
    AFNI_NOSPLASH=YES \
    AFNI_VERSION_CHECK=NO \
    AFNI_IMSAVE_WARNINGS=NO \
    AFNI_DONT_LOGFILE=YES \
    FREESURFER_VER=${FS_VER} \
    FREESURFER_HOME=/usr/pubsw/packages/freesurfer/Rocky8-R${FS_VER} \
    SUBJECTS_DIR=/usr/pubsw/packages/freesurfer/Rocky8-R${FS_VER}/subjects \
    FS_FREESURFERENV_NO_OUTPUT=1 \
    ANTSPATH=/usr/pubsw/packages/ants/${ANTS_VERSION}/bin/ \
    C3DDIR=/usr/pubsw/packages/c3d/${C3D_VERSION} \
    SYNB0_DISCO_DIR=/usr/pubsw/packages/synb0-disco/${SYNB0_VERSION} \
    TORTOISEDIR=/usr/pubsw/packages/tortoise/${TORTOISE_TAG:-local} \
    BIDS_VALIDATOR_DIR=/usr/pubsw/packages/bids-validator/${BIDS_VALIDATOR_VERSION} \
    PYTHONPATH=/usr/pubsw/packages/cmig_pipeline/python
# the values FreeSurfer's setup gives in the image, so FreeSurfer scripts
#   run without it (mri_mcadura_seg reads FREESURFER_HOME_FSPYTHON) never
#   take a value singularity passes in from the host
ENV FREESURFER=${FREESURFER_HOME} \
    FREESURFER_HOME_FSPYTHON=${FREESURFER_HOME} \
    FSFAST_HOME=${FREESURFER_HOME}/fsfast \
    FUNCTIONALS_DIR=${FREESURFER_HOME}/sessions \
    LOCAL_DIR=${FREESURFER_HOME}/local \
    MINC_BIN_DIR=${FREESURFER_HOME}/mni/bin \
    MINC_LIB_DIR=${FREESURFER_HOME}/mni/lib \
    MNI_DATAPATH=${FREESURFER_HOME}/mni/data \
    MNI_DIR=${FREESURFER_HOME}/mni \
    MNI_PERL5LIB=${FREESURFER_HOME}/mni/share/perl5 \
    PERL5LIB=${FREESURFER_HOME}/mni/share/perl5
ENV PATH=/usr/pubsw/packages/cmig_pipeline/bin:${FSLDIR}/bin:${AFNIDIR}:${FREESURFER_HOME}/bin:${ANTSPATH}:${C3DDIR}/bin:${TORTOISEDIR}/bin:${TORTOISEDIR}/bin/DRTAMAS:${BIDS_VALIDATOR_DIR}/bin:${PATH} \
    LD_LIBRARY_PATH=${AFNIDIR}

# fail the build rather than ship an image whose tools do not start
RUN for t in topup applytopup 3dvolreg 3dDeconvolve 3dREMLfit 3dFWHMx \
             3dDespike 3dTstat mri_convert mri_synthstrip recon-all reg \
             applyUnwarpB0 epic \
             mri_nu_correct.mni mri_normalize \
             applywarp convert_xfm epi_reg fast flirt fsl_abspath fslhd \
             fslinfo fslmerge imrm imtest remove_ext fslmaths \
             antsApplyTransforms antsRegistration antsRegistrationSyNQuick.sh \
             CreateDisplacementField ImageMath N4BiasFieldCorrection \
             PrintHeader c3d_affine_tool \
             ApplyTransformationToScalar2 CombineDWIsWithBmatrix \
             CreateNonlinearityDisplacementMap DisplacementFieldOperations \
             ExtractComponentFromDisplacementField ExtractDWISubset \
             ExtractImage FSLBVecsToTORTOISEBmatrix IsDirectionPositiveDet \
             ReorientImage3D SelectBestB0 SetFieldToIsocenter \
             SmoothTransitionHCPStructural TORTOISEBmatrixToFSLBVecs \
             TORTOISEProcess ZeroFieldDims jq bids-validator ; do \
      command -v $t > /dev/null || { echo "** $t not on PATH"; exit 1; } ; \
    done \
 && for t in reg applyUnwarpB0 epic ; do \
      test $(command -v $t) = /usr/pubsw/packages/cmig_pipeline/bin/$t \
        || { echo "** $t not from cmig_pipeline/bin"; exit 1; } ; \
      ! env -u LD_LIBRARY_PATH ldd $(command -v $t) | grep -E 'not found|/opt/intel|/MMPS/' \
        || { echo "** $t libraries missing or outside cmig_pipeline/lib"; exit 1; } ; \
      ! env -u LD_LIBRARY_PATH $(command -v $t) 2>&1 | grep 'error while loading' \
        || { echo "** $t does not load its libraries"; exit 1; } ; \
    done \
 && test "$(fsl_abspath /usr/pubsw/../pubsw)" = /usr/pubsw \
 && test -d ${TORTOISEDIR}/settings/mecc_settings \
 && test -s ${FREESURFER_HOME}/models/synthstrip.infant.1.pt \
 && test -x ${SYNB0_DISCO_DIR}/data_processing/prepare_input.sh \
 && test -x /usr/pubsw/packages/MMPS/apps/cmig_run_bids_module/run_cmig_run_bids_module.sh \
 && test "$(python3 -c 'import sys; print(sys.version_info[:2])')" = "(3, 12)" \
 && python3 -c 'import numpy, scipy, pandas, nibabel, matplotlib, torch, torchvision' \
 && test $(command -v bids-validator) = ${BIDS_VALIDATOR_DIR}/bin/bids-validator \
 && bids-validator --version \
 && python3 /run.py --help > /dev/null

# the package versions (versions.json: name -> version), for the status
#   files' version columns and dataset_description.json GeneratedBy
#   (cmig_sing_utils.read_versions); python's from the interpreter; an
#   empty value fails the build
ARG CMIG_IMAGE_TAG
RUN set -- \
      cmig_pipeline="${CMIG_PIPELINE_VERSION}" \
      cmig_pipeline_tag="${CMIG_PIPELINE_TAG:-local}" \
      image_tag="${CMIG_IMAGE_TAG}" \
      freesurfer="${FS_VERSION}" fsl="${FSL_VERSION}" \
      afni="${AFNI_VERSION}" afni_commit="${AFNI_COMMIT}" \
      ants="${ANTS_VERSION}" c3d="${C3D_VERSION}" mcr="${MCR_VERSION}" \
      tortoise="${TORTOISE_TAG:-local}" synb0_disco="${SYNB0_VERSION}" \
      bids_validator="${BIDS_VALIDATOR_VERSION}" deno="${DENO_VERSION}" \
      atlases="${ATLAS_TAG:-local}" \
      gradient_distortion_fields="${GDF_TAG:-local}" \
      python="$(python3 -c 'import platform; print(platform.python_version())')" \
 && for kv in "$@" ; do \
      test -n "${kv#*=}" || { echo "** versions.json: ${kv%%=*} is empty"; exit 1; } ; \
    done \
 && mkdir -p /opt/cmig \
 && python3 -c 'import json, sys; json.dump(dict(a.split("=", 1) for a in sys.argv[1:]), open("/opt/cmig/versions.json", "w"), indent=2)' "$@" \
 && python3 -c 'import json; json.load(open("/opt/cmig/versions.json"))' \
 && cat /opt/cmig/versions.json
ENV CMIG_VERSIONS_FILE=/opt/cmig/versions.json \
    CMIG_IMAGE_TAG=${CMIG_IMAGE_TAG} \
    CMIG_PIPELINE_TAG=${CMIG_PIPELINE_TAG:-local}

RUN useradd -ms /bin/tcsh appuser
USER appuser
WORKDIR /home/appuser

# where the compiled app, atlases, gradient maps and TORTOISE came from: a
#   tag and sha256, or local; the bids-validator version (all the package
#   versions in /opt/cmig/versions.json)
LABEL org.opencontainers.image.title="CMIG BIDS app" \
      org.opencontainers.image.base.name="${BASE_IMAGE}" \
      org.opencontainers.image.version="${CMIG_PIPELINE_VERSION}" \
      org.opencontainers.image.source="https://github.com/cmig-research-group/cmig_pipeline" \
      cmig.app.tag="${CMIG_PIPELINE_TAG:-local}" \
      cmig.app.sha256="${APP_SHA256}" \
      cmig.atlases.tag="${ATLAS_TAG:-local}" \
      cmig.atlases.sha256="${ATLAS_SHA256}" \
      cmig.gradient_distortion_fields.tag="${GDF_TAG:-local}" \
      cmig.gradient_distortion_fields.sha256="${GDF_SHA256}" \
      cmig.tortoise.tag="${TORTOISE_TAG:-local}" \
      cmig.tortoise.asset="${TORTOISE_ASSET}" \
      cmig.tortoise.sha256="${TORTOISE_SHA256}" \
      cmig.synb0_disco.commit="${SYNB0_COMMIT}" \
      cmig.bids_validator.version="${BIDS_VALIDATOR_VERSION}"

ENTRYPOINT ["/run.py"]
