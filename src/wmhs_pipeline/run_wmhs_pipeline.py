#!/usr/bin/env python

# Copyright 2023 Population Health Sciences, German Center for Neurodegenerative Diseases (DZNE)
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.


from __future__ import print_function

from .wmhs_pipeline import create_wmhs_pipeline

from nipype import config, logging

import os, sys,glob
import argparse
from itertools import chain

def wmhs_pipeline_wf(scans_dir, work_dir, outputdir,subject_ids,
                     device,num_threads, 
                     opp=False, owmh=False, wfname='wmhs_pipeline'):
    
    wf = create_wmhs_pipeline(scans_dir, work_dir, outputdir, subject_ids,
                              device, num_threads, opp, owmh, wfname)
    
    wf.inputs.inputnode.subject_ids = subject_ids
    return wf
    
    
def main():
    """
    Command line wrapper for wmh segmentation pipeline
    """
    parser = argparse.ArgumentParser(description='Run WMHS pipelines for FLAIR'
                                     ' T1 and/or T2 imaging data.',
                                     epilog='Example-1: {prog} -s '
                                     '~/data/scans -w ~/data/work -p 2 -t 2'
                                     ' --subjects subj1 subj2 '
                                     '\nExample-2: {prog} -s ~/data/scans -w '
                                     '~/data/work -o ~/data/output -p 2 -t 2 \n\n'
                                     .format(prog=os.path.basename(
                                             sys.argv[0])),
                                             formatter_class=argparse.
                                             RawTextHelpFormatter)

    parser.add_argument('-s', '--scansdir', help='Scans directory where data'
                        ' is already downloaded for each subject.',
                        required=True)
    
    parser.add_argument('-w', '--workdir', help='Work directory where data'
                        ' is processed for each subject in the workflow.', required=True)

    parser.add_argument('-o', '--outputdir', help='Output directory where '
                        'results will be stored.', required=True)

    parser.add_argument('--subjects', help='One or more subject IDs'
                        ' (space separated). If omitted, all subjects will be processed.',
                        default=None, required=False,
                        nargs='+', action='append')
    
    
    parser.add_argument('-r', '--only_preproc',
                        help="Only perform pre-processing",
                        required=False, action='store_true')

    parser.add_argument('-m', '--only_wmh', 
                        help="Only WMH segmentation (no regional volumetry)",
                        required=False, action='store_true')

    parser.add_argument('-b', '--debug', help='debug mode',
                        action='store_true')
    
    parser.add_argument('-p', '--processes',
                        help='overall number of parallel processes',
                        default=1, type=int)
    parser.add_argument('-g', '--ngpus',
                        help='number of gpus to use (emb-) parallel',
                        default=1, type=int)
    parser.add_argument('-gp', '--ngpuproc',
                        help='number of processes per gpu',
                        default=1, type=int)
    parser.add_argument('-d', '--device', 
                        help='deepmedic -dev flag', default='cpu', type=str)
    
    parser.add_argument('-t', '--threads', help='ITK threads',
                        default=1, type=int)
    
    parser.add_argument('-n', '--name', help='Pipeline workflow name', 
                        default='wmhs_pipeline')
    
    args = parser.parse_args()

    scans_dir = os.path.abspath(os.path.expandvars(args.scansdir))
    if not os.path.exists(scans_dir):
        raise IOError("Scans directory does not exist.")
        
    
    subject_ids = []
    
    if args.subjects:
        subject_ids = list(chain.from_iterable(args.subjects))
    else:
        subject_ids = glob.glob(scans_dir.rstrip('/') +'/*')
        subject_ids = [os.path.basename(s.rstrip('/')) for s in subject_ids]


    print ("Creating wmhs pipeline workflow...")
    work_dir = os.path.abspath(os.path.expandvars(args.workdir))
    outputdir = os.path.abspath(os.path.expandvars(args.outputdir))
    
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)        
    
    config.update_config({
        'logging': {'log_directory': args.workdir, 'log_to_file': True},
        'execution': {'job_finished_timeout' : 60,
                      'poll_sleep_duration' : 20,
                      'hash_method' : 'content',
                      'local_hash_check' : False,
                      'stop_on_first_crash':False,
                      'crashdump_dir': args.workdir,
                      'crashfile_format': 'txt'
                       },
                       })

    #config.enable_debug_mode()
    logging.update_logging(config)
    

    wmhs_pipeline = wmhs_pipeline_wf(scans_dir, work_dir, outputdir,
                                     subject_ids, args.device, args.threads,
                                     opp=args.only_preproc, owmh=args.only_wmh,
                                     wfname='wmhs_pipeline')

    # Visualize workflow
    if args.debug:
        wmhs_pipeline.write_graph(graph2use='colored', simple_form=True)



    wmhs_pipeline.run(
                            plugin='MultiProc',
                            plugin_args={'n_procs' : args.processes,
                                         'n_gpus': args.ngpus, 
                                         'ngpuproc': args.ngpuproc
                                         }
                           )


    print('Done WMHS pipeline!!!')

    
if __name__ == '__main__':
    sys.exit(main())
