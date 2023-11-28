FROM nvcr.io/nvidia/cuda:10.0-cudnn7-runtime-centos7

MAINTAINER M.Shahid <mohammad.shahid@gmail.com>

# Install git, wget, python-dev, pip, BLAS + LAPACK and other dependencies
RUN yum install -y epel-release make cmake kernel-devel zlib-devel which bzip2 ca-certificates curl git wget unzip bc gcc gcc-c++ \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/*

# Set env 
ENV CUDA_ROOT=/usr/local/cuda \
    PATH=/opt/miniconda/envs/default/bin:/opt/fsl/bin:/opt/ants/bin:/opt/freesurfer/bin:/usr/local/cuda/bin:/usr/local/nvidia:$PATH \
    ANTSPATH=/opt/ants \
    FSLDIR=/opt/fsl \
    FREESURFER_HOME=/opt/freesurfer \
    FSLOUTPUTTYPE='NIFTI_GZ' \
    LD_LIBRARY_PATH=/usr/local/cuda/lib64:/usr/lib64:$LD_LIBRARY_PATH \
    STARTUP="/opt/startup.sh"

RUN echo "Getting ANTs binaries..." \
    && cd / \
    && mkdir -p /opt/ants/bin \
    && curl -s https://netcologne.dl.sourceforge.net/project/c3d/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz |tar xzf - c3d-1.0.0-Linux-x86_64/bin/c3d_affine_tool \
    && cp c3d-1.0.0-Linux-x86_64/bin/c3d_affine_tool /opt/ants/bin \
    && rm -rf c3d-1.0.0-Linux-x86_64/bin/c3d_affine_tool \
    && wget https://github.com/ANTsX/ANTs/archive/refs/tags/v2.2.0.tar.gz \
    && tar zxf v2.2.0.tar.gz && rm -f v2.2.0.tar.gz \
    && cd ANTs-2.2.0 && mkdir build install && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX=/ANTs-2.2.0/install ../ \
    && sed -i 's/ANTS_USE_GIT_PROTOCOL:BOOL=ON/ANTS_USE_GIT_PROTOCOL:BOOL=OFF/g;' CMakeCache.txt \
    && sed -i 's/BUILD_TESTING:BOOL=ON/BUILD_TESTING:BOOL=OFF/g;' CMakeCache.txt \
    && make -j8 && cd ANTS-build && make install \
    && mkdir -p /opt/ants/bin \
    && cp /ANTs-2.2.0/install/bin/antsApplyTransforms /opt/ants/bin \
    && cp /ANTs-2.2.0/install/bin/antsRegistration* /opt/ants/bin/ \
    && cp /ANTs-2.2.0/install/bin/ConvertTransformFile /opt/ants/bin/ \
    && cp /ANTs-2.2.0/install/bin/DenoiseImage /opt/ants/bin/ \
    && cp /ANTs-2.2.0/install/bin/N4BiasFieldCorrection /opt/ants/bin/ \
    && rm -rf /ANTs-2.2.0

COPY copy_bin.sh fs_binaries.txt fsl_binaries.txt /opt/

RUN echo "Getting Freesurfer binaries..." \
    && cd /opt \
    && curl -sSL --retry 5 https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz | gunzip | tar -x -v --files-from fs_binaries.txt -f -


RUN echo "Getting FSL binaries..." \
    && wget -q https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-5.0.10-sources.tar.gz \
    && tar zxf fsl-5.0.10-sources.tar.gz \
    && cd fsl \
    && ./build || true  \
    && cd / \
    && chmod a+x /opt/copy_bin.sh && ./opt/copy_bin.sh \
    && rm -rf fsl-5.0.10-sources.tar.gz && rm -rf /fsl

COPY ./src /opt/

RUN echo "Installing miniconda and pipeline..." \
    && curl -sSL -o miniconda.sh http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash miniconda.sh -b -p /opt/miniconda \
    && rm -f miniconda.sh \
    && /opt/miniconda/bin/conda update -y -n base conda \
    && /opt/miniconda/bin/conda config --add channels conda-forge \
    && /opt/miniconda/bin/conda create -y -n default python=3.7 \
    && source /opt/miniconda/bin/activate default \
    && conda clean -y --all \
    && pip install -U --no-cache-dir pip \
    && cd /opt \
    && pip install -r requirements.txt \
    && pip install --no-cache-dir https://github.com/deepmedic/deepmedic/archive/refs/tags/v0.7.4.tar.gz \
    && pip install -e .


RUN echo "Updating settings..." \
    && source /opt/miniconda/bin/activate default \
    && cp -f /opt/wmhs_pipeline/multiproc.py /opt/miniconda/envs/default/lib/python3.7/site-packages/nipype/pipeline/plugins/ \
    && cp -f /opt/wmhs_pipeline/testSession.py /opt/miniconda/envs/default/lib/python3.7/site-packages/deepmedic/frontEnd/testSession.py \
    && if [ ! -f "$STARTUP" ]; then \
       echo '#!/usr/bin/env bash' >> $STARTUP; \
       echo 'echo "startup: setting environment"...' >> $STARTUP; \
       echo 'source /docker_environment' >> $STARTUP; \
       echo 'export LD_LIBRARY_PATH=/.singularity.d/libs:/usr/lib64:/usr/local/cuda/lib64:$LD_LIBRARY_PATH' >> $STARTUP; \
       echo 'source /opt/miniconda/bin/activate default' >> $STARTUP; \
       echo 'export TF_CPP_MIN_LOG_LEVEL=2' >> $STARTUP; \
       echo 'IFS=";" read -ra COMMANDS <<< "$@"' >> $STARTUP; \
       echo 'for command in "${COMMANDS[@]}";do' >> $STARTUP; \
       echo 'echo "`date` : executing ${command}";' >> $STARTUP; \
       echo '${command}' >> $STARTUP; \
       echo 'done' >> $STARTUP; \       
       fi \
    && chmod 755 $STARTUP 


ENTRYPOINT ["/opt/startup.sh"]

