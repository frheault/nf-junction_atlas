# nf-junction-atlas

`nf-junction-atlas` is a Nextflow pipeline designed for advanced neuroimaging analysis, specifically focusing on the generation of junction signatures and connectivity analysis in MNI space. It integrates `Tractoflow` for diffusion MRI (dMRI) preprocessing and reconstruction, along with anatomical segmentation and registration workflows.

## Features

- **Tractoflow Integration**: Automated preprocessing of DWI data (denoising, Gibbs deringing, eddy current correction, topup, etc.) and reconstruction (DTI, fODF).
- **Anatomical Segmentation**: Optional FreeSurfer `recon-all` integration and automated generation of lobes parcellations.
- **Registration**: Robust registration of anatomical and diffusion data to MNI space using ANTs.
- **Tractography**: Local tracking with customizable parameters.
- **Connectivity Analysis**: Generation of junction signatures using warped tractograms and parcellations in MNI space.

## Prerequisites

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (>= 22.10.0)
- [Docker](https://docs.docker.com/engine/installation/) or [Singularity/Apptainer](https://apptainer.org/docs/user/main/quick_start.html#quick-installation)
- A FreeSurfer license file (if running `recon-all`).

## Usage

To run the pipeline, use the following command:

```bash
nextflow run main.nf \
    --input /path/to/data \
    --fs_license /path/to/license.txt \
    --mni_template /path/to/mni_template.nii.gz \
    --ants_template /path/to/ants_template.nii.gz \
    --ants_probability_map /path/to/ants_prob_map.nii.gz \
    --all_signatures /path/to/all_signatures.json \
    --signatures_mapping /path/to/mapping.json \
    -profile docker
```

### Input Data Structure

The pipeline expects a BIDS-like input directory structure:

```text
input/
├── sub-01/
│   └── ses-01/
│       ├── sub-01_ses-01_dwi.nii.gz
│       ├── sub-01_ses-01_dwi.bval
│       ├── sub-01_ses-01_dwi.bvec
│       ├── sub-01_ses-01_t1.nii.gz
│       └── freesurfer/ (optional)
└── sub-02/
    └── ...
```

## Pipeline Steps

1.  **Data Loading**: Automatically identifies DWI and T1 files for each subject/session.
2.  **Segmentation**: Runs FreeSurfer `recon-all` (optional) or uses existing outputs to generate parcellations.
3.  **Tractoflow**: Processes dMRI data, including:
    *   Denoising and Gibbs deringing.
    *   Eddy current and Topup correction.
    *   DTI and fODF reconstruction.
    *   Registration of T1 to DWI space.
4.  **Tracking**: Performs local tracking using the reconstructed fODF and WM masks.
5.  **MNI Registration**: Registers the subjects' T1 (in DWI space) to the provided MNI template.
6.  **Transformation**: Warps all relevant images (FA, MD, NUFO, AFD, labels) and the tractogram into MNI space.
7.  **Connectivity Signatures**: Generates junction signatures using the standardized data in MNI space.

## Configuration

Default parameters are defined in `nextflow.config`. You can override them via the command line or by providing a custom configuration file.

Key parameters include:
- `run_freesurfer`: Whether to run FreeSurfer (default: `true`).
- `output`: Directory to store results (default: `results`).
- `dti_max_bvalue`: Maximum b-value for DTI fitting.
- `fodf_sh_order`: Spherical harmonics order for fODF.
- `run_local_tracking`: Enable/disable local tracking.

## Credits

This pipeline is developed and maintained by the Scilus group. It builds upon several open-source tools and libraries, including Nextflow, ANTs, MRtrix3, and the nf-core framework.
