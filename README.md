
# nf-junction-atlas

A Nextflow pipeline for white matter junction signature analysis.

## Scientific Context
This project implements the **White Matter Junction Atlas**, as described in the manuscript *"A Connectivity-Driven Multi-Scale Atlas of the Brain's White Matter Junctions"*. The atlas identifies complex intersections in human white matter where major pathway systems (Association, Projection, Commissural, and Cerebellar) converge and overlap. By mapping these "junction" territories, the atlas provides a framework for understanding white matter organization beyond individually named tracts.

## Workflow Overview
Characterizing white matter junctions is a two-step process in this repository:

1.  **Core pipeline (`main.nf`):** This Nextflow pipeline automates the extensive preprocessing required, including:
    -   Anatomical and Diffusion MRI preprocessing (similar to [Tractoflow](https://github.com/scilus/tractoflow)).
    -   Cortical and subcortical parcellation (via FreeSurfer).
    -   Registration of all data to MNI152 space.
    -   Generation of voxel-wise connectivity signatures and the base junction label volume (`*__junction_labels.nii.gz`).
2.  **Post-processing (`junction_labeler.py`):** This (optional) script takes the junction label volume produced by the pipeline and parcellates it into hierarchical classes (Broad, Simplified, Full, and Anatomical Lobes) for detailed analysis.

## Requirements
- Nextflow (>= 22.10.0)
- Singularity or Docker
- FreeSurfer License (obtainable from [FreeSurfer's website](https://surfer.nmr.mgh.harvard.edu/registration.html))

## Usage

### 1. Run the Nextflow Pipeline
Use the following command to run the core pipeline. **Note:** You must provide the path to your own FreeSurfer license file.

```bash
nextflow run main.nf \
    --input ${PATH_TO_YOUR_DATA} \
    --fs_license ${PATH_TO_YOUR_LICENSE} \
    --mni_template data/mni_template.nii.gz \
    --ants_template data/ants_t1_template.nii.gz \
    --ants_probability_map ants_t1_brain_probability_map.nii.gz \
    --all_signatures data/white_matter_junction_atlas_signature_281.txt \
    -profile docker \
    -resume
```
(see the section *Input Data Structure* below)

### 2. Run the Junction Labeler
If desired, after the pipeline completes, the junction label volumes are published in the `results` directory. Run the post-processing script on these volumes:

```bash
# Example for subject sub-003, session ses-01
python junction_labeler.py \
    results/sub-003_ses-01/GENERATE_JUNCTION_SIGNATURES/sub-003_ses-01__junction_labels.nii.gz \
    --outdir results/sub-003_ses-01/JUNCTION_HIERARCHY
```

**Output Structure:**
- `class-4/`: Broad pathway systems (Asso, Proj, Comm, Cereb).
- `class-9/`: Inter-system junction classes.
- `class-31/`: Intra-system overlap and bottleneck classes.
- `class-lobes/`: 15 anatomical lobe parcellations.

## Parameters

| Parameter | Description |
| --- | --- |
| `--input` | Path to the directory containing subject data in a simplified BIDS-like structure. |
| `--fs_license` | Path to your FreeSurfer license file. |
| `--mni_template` | Path to the MNI152 template (provided in `data/`). |
| `--ants_template` | Path to the T1 template for ANTs registration (provided in `data/`). |
| `--ants_probability_map` | Path to the T1 brain probability map for ANTs (provided in `data/`). |
| `--all_signatures` | Path to the text file containing the junction signatures list (provided in `data/`). |
| `--run_freesurfer` | Set to `true` to run FreeSurfer `recon-all` (default), or `false` if you already have the outputs. |
| `--output` | Directory to publish results (default: `results`). |

## Input Data Structure
The `--input` directory should be organized as follows:
```
input/
  ├── sub-01/
  │   └── ses-01/
  │       ├── *dwi.nii.gz
  │       ├── *dwi.bval
  │       ├── *dwi.bvec
  │       ├── *t1.nii.gz
  │       └── freesurfer/ (Optional: if --run_freesurfer false)
  │           ├── label/
  │           ├── mri/
  │           ├── surf/
  │           ├── [...]
  │           └── touch/
  └── sub-02/
```

## Credits
This pipeline uses modules and subworkflows from [nf-neuro](https://github.com/nf-neuro/modules).
