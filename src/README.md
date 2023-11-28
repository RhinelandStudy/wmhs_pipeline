# RS WMH Segmentation Pipeline

This pipeline segments WMHs using multi-modal MRI as input.
It uses an ensemble of 4 models constructed with the [Deepmedic](https://github.com/deepmedic/deepmedic) tool.
Check the *Confluence* page of the pipeline for details on the pipeline.

The pipeline can be run on the *basic* or *full* version.
The **basic** version (`--only_wmh` flag) reports the volume of white and gray matter hyperintensities in the supra- and infra-tentorial regions separately.
The **full** version (default) requires additional inputs and reports white and gray matter hyperintensity volumetry across different regions (see details in the **Input** and **Output** sections below). 

**Disclosure**:
The method is trained on hyperintensities located primarily on the white matter and basal ganglia.
However, the detection phase attempts to detect hyperintensities in the whole brain (ie, not restricted to a particular region).
Due to the lack of hyperintensities outside the white matter and basal ganglia in the training data, some hyperintensities outside these areas may be missed. 

## Requirements
The pipeline is running correctly on the following Python packages (versions):

- `deepMedic` (tag v0.7.3)
- `nipype` (version 1.2.0)

Additional Python packages: `nibabel`, `numpy`

The following software needs also to be installed:

- ANTs
- FSL FLIRT

## Execution

After cloning the repository, install the pipeline with `python setup.py install`

When the pipeline is installed, it can then be executed from the command line as `run_wmhs_pipeline -h / [ARGS]`


## Input

The pipeline requires the inputs for each subject to be in the following structure:

- `subject_id` (folder)  
  **Basic** 
  - `*FLAIR.nii.gz`: FLAIR image
  - `*T1*.nii.gz`: T1w image
  - `*T2*.nii.gz`: T2w image
  - `mri` (folder) 
    - `orig*gz`: T1 MRI in FreeSurfer
    - `aseg*gz`: Structural segmentation volume
      
  **Full**
  - `mri` (folder, cont.) 
    - `ribbon.mgz` 
  - `label` (folder)
    - `lh.aparc.annot`
    - `rh.aparc.annot`
  - `surf` (folder)
    - `lh.white`
    - `rh.white`
    - `lh.pial`
    - `rh.pial`

## Output

The output is stored in the specified folder and has the following structure for both the *basic* and *full* versions:

- `subject_id` (folder)  
  **Basic** (only-white-matter-hyperintensities flag)
  - `vols_wmh.json` global hyperintensity volumes in the whole brain and in the white matter  (see **Volumetry** section below)
  - `QCsnapshot.png` snapshot for QC containing the WMH segmentation
  - `WMH_FLAIR.nii.gz` WMH segmentation in FLAIR space
  - `FLAIR_n4_maskout_norm.nii.gz` pre-processed FLAIR image
  - `T1Warped_maskout_norm.nii.gz` pre-processed T1 image, warped to FLAIR space
  - `T2Warped_maskout_norm.nii.gz` pre-processed T2 image, warped to FLAIR space
  - `brainmask_FLAIR.nii.gz` brain mask in FLAIR space (used for image pre-processing)
  - `global_FLAIR.nii.gz` labelmap with tissue-specific supra- and infra-tentorial regions  
    
  **Full**  
  - `vols_wmh.json` global + regional WMH volumetry (see **Volumetry** section below)
  - `bulls_FLAIR.nii.gz` bullseye regions mask in FLAIR space
  - `bulls_FSURF.nii.gz` bullseye regions mask in FreeSurfer space
  - `lobes_FLAIR.nii.gz` lobes mask in FLAIR space
  - `lobes_FSURF.nii.gz` lobes mask in FreeSurfer space
  - `shells_FLAIR.nii.gz` depth shells mask in FLAIR space
  - `shells_FSURF.nii.gz` depth shells mask in FreeSurfer space
  - `vascter_tiss_FLAIR.nii.gz` vascular territories mask with thissue-specific labels in FLAIR space

## Volumetry

White and gray matter hyperintensity volumes (in cubic mm) are computed in different regions

### Global

White and gray matter hyperintensity volumes in the supra- and infra-tentorial regions, as depicted in the folowing picture:

![](images/global.png)


Variables:
- `MRI_GLOBAL_SUPRATENTORIAL_WMH_DM_R1` white matter hyperintensities in supra-tentorial region
- `MRI_GLOBAL_SUPRATENTORIAL_GMH_DM_R1` gray matter hyperintensities in supra-tentorial region
- `MRI_GLOBAL_INFRATENTORIAL_WMH_DM_R1` white matter hyperintensities in infra-tentorial region
- `MRI_GLOBAL_INFRATENTORIAL_GMH_DM_R1` gray matter hyperintensities in infra-tentorial region

Total volumes of these regions are provided in variables:
- `MRI_GLOBAL_{SUPRA|INFRA}TENTORIAL_{WM|GM}_DM_R1`

(note that the variable name for the total region volume contains `WM` or `GM` instead of `WMH` or `GMH`, as in the lesion volume)

### Shells

White and gray matter hyperintensity volumes across the 4 depth levels, as depicted in the following picture:

![](images/shells.png)

Variables:
- `MRI_SHELL_PERIVENT_WGMH_DM_R1` white and gray matter hyperintensity volume in periventricular region
- `MRI_SHELL_MEDPERIVENT_WGMH_DM_R1` white and gray matter hyperintensity volume in medial periventricular region
- `MRI_SHELL_MEDDEEP_WGMH_DM_R1` white and gray matter hyperintensity volume in medial deep region
- `MRI_SHELL_DEEP_WGMH_DM_R1` white and gray matter hyperintensity volume in deep region

(the names of the variables contain the particle `WGMH` because the shells include also basal ganglia and thalamus, which are gray matter)

Total volumes of these regions are provided in variables:
- `MRI_SHELL_{PERIVENT|MEDPERIVENT|MEDDEEP|DEEP}_DM_R1`

(note that the variable name for the total region volume does not contain `WGMH`, as in the lesion volume)


### Lobes

White and gray matter hyperintensity volumes across 4 bilateral lobes + basal ganglia and thalamus.
The regions are depicted in the following figures:

![](images/lobes_surf.png)
![](images/lobes_bgt.png)

Variables:
- `MRI_LOBE_FRONTAL_WMH_{LH|RH}_DM_R1` WMH volume in the left and right frontal lobes
- `MRI_LOBE_OCCIPITAL_WMH_{LH|RH}_DM_R1` WMH volume in the left and right parietal lobes
- `MRI_LOBE_TEMPORAL_WMH_{LH|RH}_DM_R1` WMH volume in the left and right temporal lobes
- `MRI_LOBE_PARIETAL_WMH_{LH|RH}_DM_R1` WMH volume in the left and right occipital lobes
- `MRI_LOBE_BASGANGTHAL_GMH_DM_R1` gray matter hyperintensity volume in basal ganglia and thalamus

(note that `BASGANGTHAL` is composed of gray matter, hence `GMH` is used instead of `WMH`.
Note also that `BASGANGTHAL` is not bilateral, hence no `LH` or `RH` is added)

Total volumes of these regions are provided in variables:
- `MRI_LOBE_{FRONTAL|OCCIPITAL|TEMPORAL|PARIETAL}_{LH|RH}_DM_R1`
- `MRI_LOBE_BASGANGTHAL_DM_R1`

(note that the variable name for the total region volume contain neither `WMH` nor `GMH`, as in the lesion volume)

### Bullseye

White and gray matter hyperintensity volumes across 36 regions resulting from the intersection of the 4 shells times the 9 lobes.
Results using this parcellation are amenable to be graphically represented using Bullseye plots.
The figure below shows an example of such plot containing hypothetical statistics of lesion volume across regions (the values could represent load for an individual, or a population, effect estimates, etc...). 

![](images/bullseye.jpg)

The variables are named by concatenating the name of the lobe + the shell:
- `MRI_BULL_{FRONTAL,PARIETAL,TEMPORAL,OCCIPITAL}_{PERIVENT,MEDPERIVENT,MEDDEEP,DEEP}_WMH_{LH,RH}_DM_R1` white matter hyperintensity volume in each parcel
- `MRI_BULL_BASGANGTHAL_{PERIVENT,MEDPERIVENT,MEDDEEP,DEEP}_GMH_DM_R1` gray matter hyperintensity volume across shells in the basal ganglia and thalamus

Total volumes of these regions are provided in variables:
- `MRI_BULL_{FRONTAL,PARIETAL,TEMPORAL,OCCIPITAL}_{PERIVENT,MEDPERIVENT,MEDDEEP,DEEP}_{LH,RH}_DM_R1` 
- `MRI_BULL_BASGANGTHAL_{PERIVENT,MEDPERIVENT,MEDDEEP,DEEP}_DM_R1` 

(note that the variable name for the total region volume contain neither `WMH` nor `GMH`, as in the lesion volume)


### Vascular territories

White and gray matter hyperintensity volumes across 5 bilateral vascular territories.
The regions are depicted in the following picture:

![](images/vascter2.png)
![](images/vascter1.png)

Variables:
- `MRI_VASCTER_{ACA,MCA,PCA,PONS,CEREBELLUM}_WMH_{RH,LH}_DM_R1` WMH volume in the vascular territories
- `MRI_VASCTER_{ACA,MCA,PCA,CEREBELLUM}_GMH_{RH,LH}_DM_R1` gray matter hyperintensity volume in the vascular territories (note that `PONS` has only white matter)

Total volumes of these regions are provided in variables:
- `MRI_VASCTER_{ACA,MCA,PCA,PONS,CEREBELLUM}_WM_{RH,LH}_DM_R1` 
- `MRI_VASCTER_{ACA,MCA,PCA,CEREBELLUM}_GM_{RH,LH}_DM_R1` 

(note that the variable name for the total region volume contains `WM` or `GM` instead of `WMH` or `GMH`, as in the lesion volume)


## Quality check

To facilitate QC, a snapshot is generated showing overlaid segmentation results at different sections of the image.
- `QCsnapshot.png`

The figure below shows a small section of the snapshot:

![](images/snapshot.png)
