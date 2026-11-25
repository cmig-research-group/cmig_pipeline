# cmig_pipeline

The CMIG BIDS pipeline: MATLAB modules and python wrappers that process and
analyze BIDS-formatted MRI data (sMRI, task and resting-state fMRI, dMRI), one
participant and session at a time, and combine the results across sessions.
It runs as a container image, the CMIG BIDS app, which holds the compiled
modules and every program they use.

## Getting the image

Each release has a Docker image on the GitHub Container Registry:

    docker pull ghcr.io/cmig-research-group/cmig_pipeline:<version>

For Singularity/Apptainer, build a .sif from it:

    singularity build cmig_pipeline-<version>.sif docker://ghcr.io/cmig-research-group/cmig_pipeline:<version>

The image can also be built from a clone at a release tag with
`docker/Dockerfile`, which downloads the compiled executable (the release
asset) and the other components by tag and checksum.

How to run the image, its pipelines and modules, inputs, outputs and output
names: `docker/documentation/README.txt` (also in the image, in
/documentation). A FreeSurfer license file is needed (free from
https://surfer.nmr.mgh.harvard.edu/registration.html).

## Contents

- `matlab/` the MATLAB modules (`cmig_modules/`) and the functions they use
- `python/` the module wrapper (`cmig_run_bids_module.py`), the pipelines
  (`cmig_pipeline_bids.py`, `cmig_pipeline_bids_rawqc.py`), the other module
  scripts, the job and build scripts, and tools; python 3.6
- `parms/` parameter files the code reads (color lookup tables, diffusion
  directions, fiber legends, network parcellations, column abbreviations)
- `bin/` the binaries the modules run (`reg`, `applyUnwarpB0`, `epic`),
  patched to find their libraries in `lib/`
- `lib/` those libraries, with their licenses in `lib/licenses/`
- `csh/` `run_matlab_script`, which runs the MATLAB code (not compiled)
- `docker/` the Dockerfile, its files, and the documentation and license
  files the image holds

## Releases

Versions are X.Y.Z; release notes are in `CHANGELOG.md`. Each release's
compiled executable (`cmig_run_bids_module_<version>.tar.gz`) is a release
asset; the brain atlases are a release asset of the release that changed
them (`docker/Dockerfile` names the one it uses). The executable includes MATLAB functions from a
private toolbox, so the modules here cannot be compiled from this repository
alone; the release asset is the compiled form.

## License

`docker/documentation/LICENSE.txt`.
