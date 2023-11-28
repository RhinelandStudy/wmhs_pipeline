#!/usr/bin/env python

"""
#Rhineland Study MRI Post-processing pipelines
#rs_wmhs_pipeline: Pipeline for segmentation of WMH from FLAIR scans using Deepmedic/FSL/ANTs in nipype
"""
import os
import sys
from glob import glob
if os.path.exists('MANIFEST'): os.remove('MANIFEST')


def main(**extra_args):
    from setuptools import setup
    setup(name='wmhs_pipeline',
          version='1.0.0',
          description='RhinelandStudy WMHS Pipeline',
          long_description="""RhinelandStudy processing for WMH in FLAIR scans """ + \
          """It also offers support for performing additional options to run post processing analyses.""" + \
          """More pipelines addition is work in progress.""",
          author= 'shahidm',
          author_email='mohammad.shahid@dzne.de, gerard.sanroma-gueell@dzne.de',
          url='http://www.dzne.de/',
          packages = ['wmhs_pipeline'],
          entry_points={
            'console_scripts': [
                             "run_wmhs_pipeline=wmhs_pipeline.run_wmhs_pipeline:main"
                              ]
                       },
          license='DZNE License',
          classifiers = [c.strip() for c in """\
            Development Status :: 1.0
            Intended Audience :: Developers
            Intended Audience :: Science/Research
            Operating System :: OS Independent
            Programming Language :: Python
            Topic :: Software Development
            """.splitlines() if len(c.split()) > 0],    
          maintainer = 'RheinlandStudy MRI/MRI-IT group, DZNE',
          maintainer_email = 'mohammad.shahid@dzne.de, gerard.sanroma-gueell@dzne.de',
          package_data = {'wmhs_pipeline':
		['model/*ckpt*','model/*.cfg','vascter_template/*.nii.gz']},
          install_requires=["nipype","nibabel","deepmedic"],
          **extra_args
         )

if __name__ == "__main__":
    main()

