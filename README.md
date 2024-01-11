# wmhs_pipeline

This repository contains a Nipype wrapper for  the White matter hyperintensities segmentation pipeline using T1/T2&amp;FLAIR images.

It uses an ensemble of 4 models constructed with the [Deepmedic](https://github.com/deepmedic/deepmedic) tool.

Disclosure: The method is trained on hyperintensities located primarily on the white matter and basal ganglia. However, the detection phase attempts to detect hyperintensities in the whole brain (ie, not restricted to a particular region). Due to the lack of hyperintensities outside the white matter and basal ganglia in the training data, some hyperintensities outside these areas may be missed.


If you use this wrapper please cite:

Lohner, V., Pehlivan, G., Sanroma, G., Miloschewski, A., Schirmer, M. D., Stöcker, T., ... & Breteler, M. M. (2022). Relation between sex, menopause, and white matter hyperintensities: the Rhineland study. Neurology, 99(9), e935-e943.
https://doi.org/10.1212/WNL.0000000000200782

```
@article{lohner2022relation,
  title={Relation between sex, menopause, and white matter hyperintensities: the Rhineland study},
  author={Lohner, Valerie and Pehlivan, G{\"o}khan and Sanroma, Gerard and Miloschewski, Anne and Schirmer, Markus D and St{\"o}cker, Tony and Reuter, Martin and Breteler, Monique MB},
  journal={Neurology},
  volume={99},
  number={9},
  pages={e935--e943},
  year={2022},
  publisher={AAN Enterprises}
}
```

## Build docker image

```bash

nvidia-docker build -t wmhs_pipeline -f docker/Dockerfile .


```

## Or pull from docker hub

```bash
docker pull dznerheinlandstudie/rheinlandstudie:wmhs_pipeline
```

## Required inputs:
The pipe requires the following input scans and some files from Freesurfer outputs.
1. T1.nii.gz
2. T2.nii.gz
3. FLAIR.nii.gz
4. subdirectory ```mri``` having these files:
 - mri/aseg.mgz
 - mri/orig_nu.mgz
 - mri/ribbon.mgz
5. subdirectory ```surf``` having these files:
 - surf/?h.white (?=both left and right hemispheres)
 - surf/?h.pial  (?=both left and right hemispheres)
6. subdirectory ```label``` having these files:
 - label/?h.aparc.annot (?=both left and right hemishperes)

 
## Run pipeline:

### Using docker
The pipeline can be run with docker by running the container as follows:


```bash

 nvidia-docker run --rm -v /path/to/license.txt:/opt/freesurfer/license.txt \
                 -v /path/to/input_scans:/input \
                 -v /path/to/work_folder:/work \
                 -v /path/to/output:/output \
        dznerheinlandstudie/rheinlandstudie:wmhs_pipeline \
        run_wmhs_pipeline \
        -s /input \
        --subjects test_subject_01 \
        -w /work \
        -o /output \ 
        -p 4 -t 2 -g 1 -gp 1 -d cuda

```

The command line options are described briefly if the pipeline is started with only ```-h``` option.

### Using Singulraity

The pipeline can be run with Singularity by running the singularity image as follows:

```bash


singularity build wmhs_pipeline.sif docker://dznerheinlandstudie/rheinlandstudie:wmhs_pipeline
```

When the singularit image is created, then it can be run as follows:

```bash

singularity run --nv -B /path/to/fs_license.txt:/opt/freesurfer/license.txt \
                     -B /path/to/inputdata:/input \
                     -B /path/to/work:/work \
                     -B /path/to/output:/output \
            wmhs_pipeline.sif "export TFNUM_THREADS=2;export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4;export GOTO_NUM_THREADS=2;\
            run_wmhs_pipeline \ 
                      -s /input \
                      -w /work \
                      -o /output \ 
                      -p 4 -t 2 -g 1 -gp 1 -d cuda"
```



