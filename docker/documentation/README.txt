Using the CMIG BIDS app

Contents
1. What the image holds
2. FreeSurfer license
3. Running the image
4. Pipelines and modules
5. Inputs and outputs
6. Output names
7. Examples
8. Licenses
9. Acknowledgments


1. What the image holds

The CMIG BIDS app processes and analyzes BIDS-formatted MRI data (sMRI,
task and resting-state fMRI, and dMRI), one participant and session at a
time, and combines the results across sessions. The image holds:
  the cmig bids MATLAB modules, compiled as the single app
    cmig_run_bids_module, and the python wrappers and pipelines that run them
  AFNI 26.2.06, FreeSurfer 8.2.0, FSL 6.0.4 and the MATLAB Runtime R2023a
    Update 4, so nothing else needs to be installed
  the BIDS validator 2.4.1 (bids-validator, BIDS 1.11.1), run on request
    at the group level (section 4)
  the brain atlases and gradient nonlinearity maps the modules read

Docker image: cmig_bids:<version>
Singularity image: cmig_bids-<version>.sif

Software is installed under /usr/pubsw, as on the CMIG hosts; these files
are in /documentation.


2. FreeSurfer license

FreeSurfer needs a license file, free from
  https://surfer.nmr.mgh.harvard.edu/registration.html
It is not included in the image. Every participant-level run needs it, not
only FreeSurfer recon-all (mri_synthstrip and the dMRI processing use
FreeSurfer too): bind it into the container and name it with --fs_license,
or set FS_LICENSE to its path inside the container. It is checked before
anything runs, and the run stops if neither names a readable file. The
group level does not need it.


3. Running the image

The entry point, run.py, follows the BIDS app conventions:
  run.py bids_dir output_dir {participant,group} [options]

  bids_dir         the BIDS raw data directory (read only); it must have
                   a session level, sub-<label>/ses-<label>/ (a dataset
                   without ses-* directories is refused)
  output_dir       the directory holding the output roots (section 5)
  participant      for each session, the raw data QC pipeline, then the
                   full pipeline (which runs even if raw data QC failed)
  group            across sessions: the tabulated data, then the
                   concatenated analysis outputs (section 4)

  --participant_label L ...   participants, with or without sub- (default:
                              every sub-* in bids_dir)
  --session_label S ...       sessions, with or without ses- (default: every
                              ses-* of each participant)
  --module NAME               run one module or pipeline alone (section 4)
  --fs_license FILE           the FreeSurfer license file (section 2)
  --work_dir DIR              temporary files (default /tmp in the
                              container; see below)
  --dataset_name NAME         the dataset name in the per-session output
                              roots' dataset_description.json (section 5)
  --forceflag 1               rerun steps whose outputs exist (default 0
                              at the participant level, 1 at the group
                              level)
  --verbose 1                 more messages, from the steps that take it
  --dry_flag 1                print each session's steps and their
                              arguments and the roots, write nothing, and
                              stop
  --skip_bids_validator       accepted and ignored at the participant level
                              (the validator is not run per participant;
                              section 4)

  run.py --help                  the options and the module names
  run.py --module NAME --help    the options of that module or pipeline;
                                 for a MATLAB module, each option with its
                                 description and default
  run.py --version               the cmig_pipeline version

Each session runs on its own: a failed session does not stop the others,
and the exit status is 1 if any failed. --forceflag and --verbose reach
both pipelines. Any other option is passed to each step that takes it (e.g.
--task mid nback to the full pipeline, --pc_flag 1 to the raw data QC
pipeline); an option no step to be run takes is an error, before anything
runs. A value beginning with - is given as
--option=value (e.g. --fsurf_options=-autorecon1). Options not given take
the step's own defaults.

--work_dir DIR holds the temporary files: the unzipped copy of each
compressed volume the MATLAB modules read (up to ~800 MB for a 4D fMRI
series, deleted after it is read), the temporary volumes of motion
estimation and dMRI model fitting, and the group run's list of visits (the
SCRATCHROOT of the modules); and the MATLAB Runtime's cache, DIR/mcr_cache
(the extracted app, several hundred MB, reused; MCR_CACHE_ROOT). It is
created if missing. The default, /tmp in the container, is inside the
container's own filesystem with docker; to put these files on a disk of
your choice, bind a directory and name it with --work_dir.

Directories are passed as paths inside the container, so bind each input
and output directory: bids_dir read-only, and output_dir on its own (not a
directory that holds both).

Docker:
  docker run --rm \
    -v /data/rawdata:/data/rawdata:ro \
    -v /data/derivatives:/data/derivatives \
    -v /path/to/license.txt:/license.txt:ro \
    cmig_bids:<version> /data/rawdata /data/derivatives participant \
      [options] --fs_license /license.txt

Docker runs as the image's user, appuser; the output directories must be
writable by it (or add --user $(id -u):$(id -g)).

Singularity:
  singularity run --cleanenv --no-home \
    -B /data/rawdata:/data/rawdata:ro \
    -B /data/derivatives:/data/derivatives \
    -B /path/to/license.txt:/license.txt:ro \
    cmig_bids-<version>.sif /data/rawdata /data/derivatives participant \
      [options] --fs_license /license.txt

Singularity runs as the calling user. --cleanenv keeps the host's
environment variables out of the container (by default Singularity passes
them in, and host settings for FreeSurfer, FSL or python can override the
image's), and --no-home keeps out the home directory (the shells the app
starts would otherwise read the user's ~/.cshrc, and python the user's
~/.local packages). With --no-home, a directory under the home directory is
not visible unless bound with -B.


4. Pipelines and modules

Pipelines (one participant and session; a failed step skips only the steps
that need its output):
  cmig_pipeline_bids         sMRI, fMRI and dMRI: proc anat, func and dwi,
                             FreeSurfer recon-all, analyze anat, analyze
                             func for each task with data (rest, nback,
                             mid, sst; or those named with --task), analyze
                             dwi, postqc anat, fsurf, func, dwi and fmap;
                             the behavioral metrics from the events files
                             first
  cmig_pipeline_bids_rawqc   raw data QC for anat, func, dwi and fmap, and
                             with --pc_flag 1 protocol compliance

Modules (each one step of one session, at the participant level):
  cmig_proc_bids_anat, cmig_proc_bids_func, cmig_proc_bids_dwi
  cmig_fsurf_bids            FreeSurfer recon-all
  cmig_analyze_bids_anat, cmig_analyze_bids_func, cmig_analyze_bids_dwi
  cmig_analyze_bids_behav    behavioral metrics of the task events
  cmig_rawqc_bids (--modality anat, func, dwi or fmap; or
    cmig_rawqc_bids_anat, _func, _dwi, _fmap)
  cmig_pc_bids_abcd          ABCD protocol compliance
  cmig_check_event_pcqc      protocol compliance merged with raw data QC
  cmig_postqc_bids_anat, cmig_postqc_bids_func, cmig_postqc_bids_dwi,
    cmig_postqc_bids_fsurf, cmig_postqc_bids_fmap

Group level (across sessions; without --module, tabdat then concat):
  cmig_compile_bids_tabdat   the per-session tables of the raw data QC,
                             analysis and post-processing QC roots, one
                             table per release table, and the new 8.0
                             tables beside them (e.g. mr_y_qc__raw__fmri,
                             the series tables mr_y_qc__raw__series and
                             mr_y_qc__post__series); needs the release's
                             per-table json directory, bids_dir/phenotype,
                             or another named with --phenotype_dir; with
                             --module, run alone (minutes)
  cmig_concat_bids_sessions  the analysis outputs of many sessions, one
                             file per measure, for each datatype (or those
                             named with --datatype), in the release 7.0
                             layout (below); with --module, run alone
                             (hours)
  bids_validator             the BIDS validator on bids_dir and on each
                             per-session output root given by its option
                             (section 5; not tabulated or concatenated);
                             run only with --module bids_validator, never
                             in the default group run
--participant_label and --session_label, or --id_events_file FILE (a csv
with participant_id and session_id columns), limit tabdat and concat to the
visits they select.

The group level runs with --forceflag 1 by default: a rerun (with more
sessions, or other labels) rebuilds the tables and the concatenated files;
--forceflag 0 keeps those that exist. A concat datatype with no data (e.g.
a dataset without dMRI) is skipped with a WARNING, whether or not it is
named with --datatype; no folder is made for it.

Each table is written to tabulated_dir as <table>.tsv, named for the
table (e.g. mr_y_qc__raw__dmri.tsv): the tabulated and concatenated files
are not BIDS derivatives, and the BIDS file name rules (section 6) apply
only to the per-session roots. Beside the tables, compile_tabdat_report.txt
lists each table's rows and its release and new columns, the release
columns not found, the columns placed in no table, the files skipped and
any conflicting values.

The column names follow the release 8.0 data dictionary, some of them
renamed from release 7.0 (listed in the CHANGELOG). One release 7.0 name
holds a different value: the 12 fMRI FWHM columns
mr_y_qc__raw__{rsfmri,tfmri__<task>}__{x,y,z}dir_fwhm are now the
smoothness of the data as imported; the detrended value they held in
release 7.0 is __dtr__{x,y,z}dir_fwhm (with its __ok and __comp forms,
which take the dtr__ too). Tables from release 7.0 and 8.0 are not to be
merged on these names.

The concatenated files follow the release 7.0 layout under concat_dir,
one folder per datatype, with the fMRI correlations and time series split
by task:
  smri_measures    vertexwise/smri     smri_surfaces    vertexwise/surf
  tfmri_stats      vertexwise/tfmri    dmri_dti         vertexwise/dti
  dmri_isdti       vertexwise/isdti    dmri_rsi         vertexwise/rsi
  fmri_corr, fmri_tseries
                   corrmat/rsfmri (rest), corrmat/tfmri (the other tasks)
Each folder has one vol_info.mat, the session rows of every file in it;
fmri_corr and fmri_tseries share those of corrmat/, so they are to be run
on the same sessions.

Each concatenated <stem>.mat has a <stem>_status.tsv beside it (errcode,
errmsg, and the sessions with and without an input file, and those padded
or truncated on the run and time axes), and each concat run writes
concat_<datatype>_status.tsv in concat_dir, with its errcode, errmsg and
the numbers of files written, kept and without data.

The BIDS validator writes its report, for each directory it checks, to
  output_dir/logs/bids_validator_<root>_<YYMMDDHHMM>.json
and the exit status is its error code (0 if every directory passed).
Validating all of a large dataset's raw data takes hours and much memory;
run it as a job of its own.


5. Inputs and outputs

Input: a BIDS raw data directory (bids_dir), sub-<label>/ses-<label>/ with
anat, func, dwi and fmap; task fMRI analysis also reads the events files
(_events.tsv). The session level is required: every module, output path,
file name and status file is keyed by session, so a bids_dir without ses-*
directories is refused.

Outputs, under output_dir unless named by its option; each per-session
root is a BIDS derivatives directory (with its dataset_description.json):
  cmig_rawqc      --rawqc_rootdir     raw data QC and protocol compliance
                  --pc_rootdir        protocol compliance (default: the
                                      rawqc root)
  cmig_proc       --proc_rootdir      processed images
  freesurfer      --fsurf_rootdir     FreeSurfer recons, one per session
  cmig_analysis   --analysis_rootdir  analysis results
  cmig_postqc     --postqc_rootdir    post-processing QC
  tabulated       --tabulated_dir     group level: the compiled tables
  concatenated    --concat_dir        group level: the concatenated files
Each step is given the roots it reads and the one it writes under its own
option names. A root that does not exist is created for the outputs; an
input root named by its option must exist.

The group roots, tabulated and concatenated, are not laid out as BIDS
datasets (no sub-*/ses-*): they get no dataset_description.json, and the
BIDS validator does not check them.

A per-session root without a dataset_description.json is given one (an
existing one is kept), whose Name is "<dataset_name> - <root>", e.g. "ABCD
Study - cmig_proc", with --dataset_name; without it, the Name in bids_dir's
dataset_description.json, with a trailing " - Raw data" removed, takes the
place of <dataset_name>. Its GeneratedBy lists cmig (the cmig_pipeline
version; the freesurfer root lists freesurfer first) and then one entry per
package of /opt/cmig/versions.json (below), each with its Name and Version;
it records the image that first wrote the root, since an existing
dataset_description.json is kept.

The image's package versions are in /opt/cmig/versions.json (named by
$CMIG_VERSIONS_FILE), written when the image is built: cmig_pipeline,
cmig_pipeline_tag, image_tag, freesurfer, fsl, afni, afni_commit, ants, c3d,
mcr (the MATLAB Runtime), tortoise, synb0_disco, bids_validator, deno,
atlases, gradient_distortion_fields and python ("local" for a piece built
from an untagged copy).

Each step writes a status file beside its outputs,
  <participant>_<session>_..._desc-<step>_status.tsv
whose errcode is 0 when the step succeeded and whose errmsg says why it did
not (n/a if nothing went wrong). A step whose outputs exist is not rerun
unless forced (--forceflag 1, for the steps that take it). The columns after
errmsg record what wrote the file, n/a where a value does not apply:
  app_name            the module, or the program that did the work (e.g.
                      recon-all)
  app_version         the image tag (e.g. cmig_bids:<version>)
  app_sha256          the sha256 of the .sif, where the job that ran it was
                      given one (CMIG_IMAGE_SHA256)
  pipeline_name       cmig
  software_name       cmig_pipeline
  software_version    the cmig_pipeline version
  <package>_version   one per package of versions.json: freesurfer_version,
                      fsl_version, afni_version, ants_version, c3d_version,
                      mcr_version, tortoise_version, synb0_version,
                      bids_validator_version, deno_version, atlases_version,
                      gradient_distortion_fields_version, python_version
Run outside the image (no versions.json), the package columns are n/a.


6. Output names

File names follow the BIDS v1.11.1 derivative rules: <entities>_<suffix>.<ext>,
the entities in BIDS order (sub, ses, task, acq, rec, dir, run, hemi, space,
atlas, ..., desc), desc- always last, values alphanumeric. A file keeps the
entities of the scan it comes from; outputs of a whole session carry only
sub- and ses-, and those of a task also task-. Below, <s> stands for
sub-<label>_ses-<label>.

  desc-preproc     a processed image whose name would otherwise be a raw
                   data name: _desc-preproc_T1w, _T2w, _bold, _epi, _dwi
                   (the QC modules drop it from the names they derive)
  desc-<cond>      the condition of a task fMRI statistic, e.g.
                   task-nback_run-01_..._desc-2b_beta
  atlas-<label>    the atlas of ROI and surface data:
                     gp         Gordon parcels
                     gpnet      Gordon networks (groups of gp parcels)
                     gpaseg     Gordon parcels and FreeSurfer aseg ROIs
                     gpnetaseg  Gordon networks and aseg ROIs
                     dsk        Desikan-Killiany cortical parcels
                     dskaseg    Desikan-Killiany and aseg ROIs
                     dst        Destrieux cortical parcels
                     dstaseg    Destrieux and aseg ROIs
                     ico<N>     icosahedral surface mesh of order N
                   (also aseg, at for AtlasTrack fiber tracts, fzy for
                   fuzzy-cluster parcels, gpdstdskaseg)
  hemi-L, hemi-R   the hemisphere of surface data
  GIfTI            compound extensions by content: .shape.gii (thickness,
                   area, sulcal depth), .func.gii (sampled intensities,
                   fMRI statistics, coverage masks)

Documented deviations from BIDS:
  .mat             MATLAB files: per session saved -v7 (readable with
                   scipy.io.loadmat), the concatenated files -v7.3 (HDF5,
                   readable with h5py); .tsv or .gii copies exist where
                   practical (_corrmat.tsv, the .shape.gii and .func.gii
                   maps)
  custom suffixes  status, the QC tables (rawqc, postqc, pc, pcqc,
                   manualrawqc, manualpostqc; the module is the suffix and
                   the table type is in desc-, e.g. desc-session_postqc),
                   corr (a tabulated row) and corrmat (a matrix), tseries,
                   beta, sem, fstat, deconv, var, motionqc, behavioral, the
                   measures (thk, area, sulc, t1gm, dtiFA, rsiRND, ...),
                   B0dx, and the jpg names (the render or plot type, e.g.
                   desc-3view_bold, motionfd, regT1, relRMSE, tsnr)
  atlas-ico<N>     fsaverage truncated to its first 10*4^N+2 vertices, in
                   place of space-/den-. ico5 and ico6 data display on
                   fsaverage5 and fsaverage6 (the same vertices); ico3 and
                   ico4 data do not display on fsaverage3 and fsaverage4
                   (another vertex order); ico7 is fsaverage
  no space-        ROI tables are in the participant's own space
  run-             per-run files carry the BIDS run number; the columns of
                   per-run data carry the slot (r01, r02), as in the ABCD
                   release 7.0 names
  session root     QC files of a whole session that span datatypes (pc,
                   pcqc, session info, manual QC, the raw data QC status)
                   sit in sub-<label>/ses-<label>/ outside any datatype
                   directory; listed in the root's .bidsignore
  freesurfer       FreeSurfer's own layout, one directory per session,
                   sub-<label>_ses-<label>/ (mri/, surf/, label/, stats/
                   ...), and fsaverage/; listed in its .bidsignore

Each per-session root has a .bidsignore, written with its
dataset_description.json, listing the root's non-BIDS content.

Image and surface outputs have a JSON sidecar with Description, Sources
and SpatialReference (Type and SkullStripped for masks). SpatialReference
is a BIDS URI or a template name: fsaverage5 or fsaverage6 for ico5 and
ico6 data, and fsaverage for the other orders, with the truncation stated
in the Description (e.g. "fsaverage's first 10242 vertices (ico5)").
The atlas sidecars (cmig_proc atlases/) name the file of the atlas package
they come from, e.g. Sources
["bids:cmig_pipeline_atlases:T1_Atlas/T1_atlas_ABCD.mgz"]. A processed
field map's IntendedFor names the processed bold or dwi it corrects.

BIDS URIs name another file: bids::<path> within the sidecar's own root,
bids:<root>:<path> in another root, e.g.
bids:cmig_proc:sub-01/ses-01/anat/sub-01_ses-01_run-01_desc-preproc_T1w.nii.gz
The root labels (rawdata, cmig_proc, freesurfer, cmig_analysis,
cmig_rawqc, cmig_postqc, cmig_pc) resolve through the DatasetLinks of
each root's dataset_description.json, which lists the roots that exist;
cmig_proc's also links cmig_pipeline_atlases to the atlas package, the
release asset of the image's atlases version (none for an untagged one).

Tables are tab-separated with a header row, without quotes; an empty or
missing value is n/a.

Job logs: jobs run by the CMIG rounds scheduler (cmig_run_bids_rounds.py)
have their .out and .err copied to each root they write, as
  <root>/logs/sub-<label>/ses-<label>/<s>_desc-<step><YYMMDDHHMM>.{out,err}
(e.g. desc-pipelinerawqc2611121530.out; <step> is the script name without
cmig_ and _bids); every run's logs are kept.

Examples (paths under each root, from sub-<label>/ses-<label>/):
  cmig_rawqc
    func/<s>_task-rest_run-01_desc-3view_bold.jpg
    func/<s>_task-rest_run-01_motionfd.jpg
    dwi/<s>_run-01_desc-DTIb1000_relRMSE.jpg
    func/<s>_desc-bold_rawqc.tsv
    <s>_desc-series_pc.tsv, <s>_desc-reviews_manualrawqc.tsv
  cmig_proc
    anat/<s>_run-01_desc-preproc_T1w.nii.gz (and .json)
    anat/<s>_run-01_space-T1w_desc-brain_mask.nii.gz
    func/<s>_task-rest_run-01_desc-preproc_bold.nii.gz
    func/<s>_task-rest_run-01_desc-motion_timeseries.tsv
    func/<s>_task-nback_run-01_events.tsv (a copy of the raw data's, read
      by task fMRI analysis)
    fmap/<s>_acq-func_dir-AP_run-01_desc-preproc_epi.nii.gz (and its
      _desc-motion_timeseries.tsv)
    fmap/<s>_acq-func_desc-AP01PA02_B0dx.nii.gz
    dwi/<s>_run-01_desc-preproc_dwi.nii.gz (.bval, .bvec, .json)
    atlases/atlas-ABCD_T1w.nii.gz (at the root)
  freesurfer (at the root)
    <s>/surf/lh.white, <s>/<s>_desc-fsurf_status.tsv
  cmig_analysis
    anat/<s>_atlas-dsk_thk.tsv
    anat/<s>_hemi-L_atlas-ico7_desc-fwhm<F>_thk.shape.gii
    func/<s>_task-nback_run-01_atlas-dsk_desc-2b_deconv.tsv
    func/<s>_task-nback_hemi-L_atlas-ico5_desc-2b_beta.func.gii (and .json)
    func/<s>_task-rest_atlas-gpaseg_corrmat.tsv
    dwi/<s>_atlas-at_dtiFA.tsv, dwi/<s>_motionqc.tsv
    atlas-at_dseg.tsv (at the root)
  cmig_postqc
    func/<s>_desc-bold_postqc.tsv, func/<s>_desc-inter_tsnr.jpg
    func/<s>_space-T1w_desc-meanNorm_boldref.nii.gz
    fmap/<s>_acq-func_desc-APvsPA_diff.jpg
    <s>_desc-session_manualpostqc.tsv


7. Examples

Raw data QC and the whole pipeline for one session:
  docker run --rm -v /data/rawdata:/data/rawdata:ro \
    -v /data/derivatives:/data/derivatives \
    -v ~/license.txt:/license.txt:ro \
    cmig_bids:<version> /data/rawdata /data/derivatives participant \
      --participant_label 01 --session_label 01 \
      --fs_license /license.txt

Raw data QC with protocol compliance, for every session, with the
temporary files on a scratch disk:
  docker run --rm -v /data/rawdata:/data/rawdata:ro \
    -v /data/derivatives:/data/derivatives \
    -v /scratch/cmig:/work \
    -v ~/license.txt:/license.txt:ro \
    cmig_bids:<version> /data/rawdata /data/derivatives participant \
      --module cmig_pipeline_bids_rawqc --pc_flag 1 \
      --work_dir /work --fs_license /license.txt

One module, resting-state fMRI analysis with 5 and 10 minute trimmed
correlations:
  singularity run --cleanenv --no-home \
    -B /data/rawdata:/data/rawdata:ro \
    -B /data/derivatives:/data/derivatives \
    -B ~/license.txt:/license.txt:ro \
    cmig_bids-<version>.sif /data/rawdata /data/derivatives participant \
      --participant_label 01 --session_label 01 \
      --module cmig_analyze_bids_func \
      --task rest --trim_minutes 5 --rest_trim_minutes 10 \
      --fs_license /license.txt

Combining the resting-state correlations of all analyzed sessions:
  docker run --rm -v /data/rawdata:/data/rawdata:ro \
    -v /data/derivatives:/data/derivatives \
    cmig_bids:<version> /data/rawdata /data/derivatives group \
      --module cmig_concat_bids_sessions --datatype fmri_corr \
      --tasklist rest --atlaslist gpaseg

Validating the raw data and the processed images:
  docker run --rm -v /data/rawdata:/data/rawdata:ro \
    -v /data/derivatives:/data/derivatives \
    cmig_bids:<version> /data/rawdata /data/derivatives group \
      --module bids_validator --proc_rootdir /data/derivatives/cmig_proc


8. Licenses

See /documentation/LICENSE.txt and /documentation/licenses. The app's
code written at UC San Diego (2006 to April 2025) is Copyright The Regents
of the University of California, under UC San Diego's notice for
educational, research and non-profit use; the code written since at the
J. Craig Venter Institute is under the PolyForm Noncommercial License
1.0.0. The brain atlases (2010-2020) and the gradient nonlinearity maps
(2006-2025), made at UC San Diego, are Copyright The Regents of the
University of California, under UC San Diego's notice adapted for works,
for educational, research and non-profit use.


9. Acknowledgments

Developed at the University of California, San Diego (2006 to April 2025)
and at the J. Craig Venter Institute (JCVI; since April 2025).

Authors: Donald J. Hagler, Jr., Diliana Pecheva, Anders Dale, Vijay
Venkatraman, Payam Saisan, Ben Cipollini, Jason Sherfey, Feng Xue, Mingxiong
Huang, Jose Teruel, Dani Cornejo, and others.

Contributors: Octavio Ruiz, Mai Bui, Sean Hatton, Tyler Berkness, Mark
Fiecas, Josh Kuperman, Emma Pearson, Bader Chaarani, Nate White.

Code shared by others:
  The dMRI processing script, cmig_proc_dwi.bash, is based on the
    diffusion MRI processing script for ABCD and HBCD data by M. Okan
    Irfanoglu (National Institutes of Health), with changes by Diliana
    Pecheva (JCVI).
  cmig_synb0_disco.sh is based on Synb0-DISCO's pipeline script
    (src/pipeline.sh; MASILab, Vanderbilt University).
  abcd_resp_filter_motion.m and mmil_censor_motion.m (respiratory motion
    filter and censoring) are based on code shared by Eric Earl (OHSU).
  The E-Prime extraction (abcd_extract_eprime_*.m, and the python
    cmig_abcd_eprime_to_events.py, cmig_abcd_events_to_behavioral.py and
    cmig_abcd_events_to_timing.py derived from it) is based on code
    provided by ABCD sites: abcd_mid_extract.m by Mary Soules (University of
    Michigan; MID), Nback_scanner.py by Eric Feczko (OHSU; n-back),
    make-SST-refs.py by Michael Riedel (FIU; SST), and python code by Sage
    Hahn (UVM; SST glitch and violator variables).
