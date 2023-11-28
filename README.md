# wmhs_pipeline
White matter hyperintensities segmentation pipeline using T1/T2&amp;FLAIR images


## build docker image

```bash

nvidia-docker build -t wmhs_pipeline -f docker/Dockerfile .


```

## run docker container

```bash

 docker run --rm -v /path/license.txt:/opt/freesurfer/license.txt \
                 -v /path/input_scans:/input \
                 -v /path/work_folder:/work \
                 -v /path/output:/output \
        freesurfer6_pipeline run_fs_pipeline \
        -s /input \
        --subjects test_subject_01 \
        -w /work \
        -o /output \ 
        -a 3T qcache -fT1T2 -p 4 -t 2

```

