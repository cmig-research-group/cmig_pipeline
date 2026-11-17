# cmig_pipeline_dev

Development repository of the CMIG BIDS pipeline: MATLAB modules and python
wrappers that process and analyze BIDS-formatted MRI data (sMRI, task and
resting-state fMRI, dMRI), one participant and session at a time, and the
Dockerfile of the CMIG BIDS app image. Private. Each release is copied to the
public cmig_pipeline repository, with the compiled executable as its release
asset.

## Contents

- `matlab/` the MATLAB modules (`cmig_modules/`) and the functions they use,
  in the subdirectories they have in MMPS
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
  files the image holds (`docker/documentation/README.txt` describes using
  the image)

## Dependencies

- surface_morphometry_toolbox (private): MATLAB functions compiled into the
  executable
- gradient_distortion_fields: the gradient nonlinearity distortion fields
  and `gradient_distortion_files.m`
- the brain atlases, a release asset
- FreeSurfer, FSL, AFNI and the MATLAB Runtime, which the image installs at
  pinned versions

## Releases

Versions are X.Y.Z. Release notes are kept in `CHANGELOG.md`, one section per
version, added to as changes are made.

## License

`docker/documentation/LICENSE.txt`.
