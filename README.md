# nf-junction-atlas

A Nextflow pipeline for white matter junction signature analysis.

## Scientific Context
This project implements the **White Matter Junction Atlas**, as described in the manuscript *"A Connectivity-Driven Multi-Scale Atlas of the Brain's White Matter Junctions"*. The atlas identifies complex intersections in human white matter where major pathway systems (Association, Projection, Commissural, and Cerebellar) converge and overlap. By mapping these "junction" territories, the atlas provides a framework for understanding white matter organization beyond individually named tracts.

## Workflow Overview
Characterizing white matter junctions is a two-step process in this repository:

1.  **Core Pipeline (`main.nf`):** This Nextflow pipeline automates the extensive preprocessing required, including:
    -   Anatomical and Diffusion MRI preprocessing (via [Tractoflow](https://github.com/scilus/tractoflow)).
    -   Cortical and subcortical parcellation (via FreeSurfer).
    -   Registration of all data to MNI152 space.
    -   Generation of voxel-wise connectivity signatures and the base junction label volume (`*__junction_labels.nii.gz`).
2.  **Post-Processing & QA (`junction_labeler.py` & `junction_qa.py`):** These scripts take the junction label volume produced by the pipeline and parcellate it into hierarchical classes (Broad, Simplified, Full, and Anatomical Lobes), as well as generate an interpretable visual QA report for detailed analysis.

## Requirements
- Nextflow
- Container Platform: Docker, Singularity, or Apptainer (strongly recommended on HPC clusters).
- FreeSurfer License (obtainable from [FreeSurfer's website](https://surfer.nmr.mgh.harvard.edu/registration.html))

## Usage

### 1. Run the Nextflow Pipeline
Use the following command to run the core pipeline. **Note:** You must provide the path to your own FreeSurfer license file. 

If running on an HPC server, we strongly recommend using the **Apptainer** profile (`-profile apptainer`) rather than Docker.

> [!CAUTION]
> If you’re trying to run the pipeline on an HPC cluster and you haven’t downloaded and installed the containers to run the pipeline offline, we suggest doing so before running the pipeline as detailed [here](https://scilus.github.io/sf-tractomics/user_guides/usage/).
> 
> Once you have pre-installed the containers, do not forget to specify the following environment variables so that nextflow can appropriately find the installed containers/images.
> 
> ```bash
> export APPTAINER_CACHEDIR=/scratch/${USER}/sf-tractomics-containers
> export NXF_APPTAINER_CACHEDIR=/scratch/${USER}/sf-tractomics-containers/cache
> export SINGULARITY_CACHEDIR=${APPTAINER_CACHEDIR}
> export NXF_SINGULARITY_CACHEDIR=${NXF_APPTAINER_CACHEDIR}
> ```

```bash
nextflow run main.nf \
    --input test/raw \
    --fs_license /path/to/your/license.txt \
    --mni_template data/mni_template.nii.gz \
    --ants_template data/ants_t1_template.nii.gz \
    --ants_probability_map data/ants_t1_brain_probability_map.nii.gz \
    --all_signatures data/white_matter_junction_atlas_signature_281.txt \
    -profile apptainer \
    -resume
```

### 2. Run the Junction QA script (Automated Parcellation and QA)
After the pipeline completes, the results are published in the `results` directory. 

We highly recommend running the `junction_qa.py` script as your post-processing step. It automatically handles batch processing across all subjects, calling `junction_labeler.py` to generate the hierarchical parcellations, and then generating a rich visual QA PDF for each subject.

```bash
python junction_qa.py \
    --in_dir results/ \
    --out_dir results/QA/
```

> **Note**: If you only want the raw NIfTI parcellations without the QA PDFs, or you want to process a single subject manually, you can call `junction_labeler.py` directly:
> ```bash
> python junction_labeler.py \
>     results/sub-003_ses-01/GENERATE_JUNCTION_SIGNATURES/sub-003_ses-01__junction_labels.nii.gz \
>     --out_dir results/sub-003_ses-01/WM_JUNCTION
> ```

**Output Structure:**
For each subject `sub-XX`, the `junction_qa.py` script will create:
- `results/sub-XX/WM_JUNCTION/Class-4/`: Broad pathway systems (Asso, Proj, Comm, Cereb).
- `results/sub-XX/WM_JUNCTION/Class-9/`: Inter-system junction classes.
- `results/sub-XX/WM_JUNCTION/Class-31/`: Intra-system overlap and bottleneck classes.
- `results/sub-XX/WM_JUNCTION/Class-lobes/`: 15 anatomical lobe parcellations.
- `results/QA/sub-XX/`: Quality assurance PDF report and label volume statistics.

## Parameters

| Parameter | Description |
| --- | --- |
| `--input` | Path to the directory containing subject data in a BIDS-like structure. Expected files: `*dwi.nii.gz`, `*dwi.bval`, `*dwi.bvec`, `*t1.nii.gz`. |
| `--fs_license` | Path to your FreeSurfer license file. |
| `--mni_template` | Path to the MNI152 template. |
| `--ants_template` | Path to the T1 template for ANTs registration. |
| `--ants_probability_map` | Path to the T1 brain probability map for ANTs. |
| `--all_signatures` | Path to the text file containing the junction signatures list. |
| `--run_freesurfer` | Set to `true` to run FreeSurfer `recon-all` (default), or `false` if you already have the outputs. |
| `--output` | Directory to publish results (default: `results`). |

## Input Data Structure
The `--input` directory should be organized as follows (please respect the expected lowercase naming structure!):
```
input/
├── sub-01/
│   └── ses-01/
│       ├── *dwi.nii.gz
│       ├── *dwi.bval
│       ├── *dwi.bvec
│       ├── *t1.nii.gz
│       └── freesurfer/ (Optional: if --run_freesurfer false)
```

## Credits
This pipeline uses modules and subworkflows from [nf-neuro](https://github.com/nf-neuro/modules).