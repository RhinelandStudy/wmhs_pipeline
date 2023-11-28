from __future__ import division

import nipype.pipeline.engine as pe
from nipype import SelectFiles
import nipype.interfaces.utility as util

from nipype.interfaces import fsl
from nipype.interfaces.fsl.utils import Reorient2Std

import nipype.interfaces.ants as ants

from nipype import IdentityInterface, DataSink

from .utils import *

from .configoptions import DM_MODEL_DIR, VASCTER_TEMPLATE_DIR
import os


def create_wmhs_pipeline(scans_dir, work_dir, outputdir, subject_ids, device,
                         num_threads, opp=False, owmh=False,
                         name='wmhs_preproc'):
    

    wmhsppwf = pe.Workflow(name=name)
    wmhsppwf.base_dir = work_dir

    inputnode = pe.Node(interface=IdentityInterface(
            fields=['subject_ids', 'outputdir']),name='inputnode')
    
    inputnode.iterables = [('subject_ids', subject_ids)]
    inputnode.inputs.subject_ids = subject_ids
    inputnode.inputs.outputdir = outputdir

    #template for input files
    templates = {"FLAIR": "{subject_id}/*FLAIR.nii.gz",
                 "T1": "{subject_id}/*T1*.nii.gz",
                 "T2": "{subject_id}/*T2*.nii.gz",
                 "T1FS": "{subject_id}/mri/orig*gz",
                 "ASEG": "{subject_id}/mri/aseg*gz"}

    # if bullseye parcels to be computed
    if not owmh:
        templates.update({"RIBBON": "{subject_id}/mri/ribbon.mgz",
                         "ANNOT_LH": "{subject_id}/label/lh.aparc.annot",
                         "ANNOT_RH": "{subject_id}/label/rh.aparc.annot",
                         "WHITE_LH": "{subject_id}/surf/lh.white",
                         "WHITE_RH": "{subject_id}/surf/rh.white",
                         "PIAL_LH": "{subject_id}/surf/lh.pial",
                         "PIAL_RH": "{subject_id}/surf/rh.pial",
                         "subject_id": "{subject_id}"})


    fileselector = pe.Node(SelectFiles(templates), name='fileselect')
    fileselector.inputs.base_directory = scans_dir


    #%% step-1a convert T1 mgz to T1.nii.gz if mgz
    convert_t1_mgz = pe.Node(interface=util.Function(
            input_names=['in_file'], output_names=['out_file'],
            function=convert_mgz), name='convert_t1_mgz')
    
    # step-1b convert aseg mgz to aseg.nii.gz if in mgz format
    convert_aseg_mgz = pe.Node(interface=util.Function(
            input_names=['in_file'], output_names=['out_file'],
            function=convert_mgz), name='convert_aseg_mgz')
    
    # step 1-c reorient 2 std
    reorient2std_fst1 = pe.Node(interface=Reorient2Std(),
                                name= 'reorient2std_fst1')
    reorient2std_aseg = pe.Node(interface=Reorient2Std(),
                                name= 'reorient2std_aseg')

    
    #%% step-3: N4BiasFieldCorrect
    #step 3a: N4 FLAIR
    n4biasfieldcorrect_fl = pe.Node(interface=ants.N4BiasFieldCorrection(),
                                    name='n4biascorrect_fl')
    n4biasfieldcorrect_fl.inputs.dimension = 3
    n4biasfieldcorrect_fl.inputs.n_iterations = [50, 50, 30, 20]
    # n4biasfieldcorrect_fl.inputs.n_iterations = [5, 5, 3, 2]
    n4biasfieldcorrect_fl.inputs.convergence_threshold = 1e-6
    n4biasfieldcorrect_fl.inputs.bspline_fitting_distance = 300
    n4biasfieldcorrect_fl.inputs.output_image='FLAIR_n4.nii.gz'
    n4biasfieldcorrect_fl.inputs.num_threads = num_threads
    
    #%% step-3b N4 T1
    n4biasfieldcorrect_t1 = pe.Node(interface=ants.N4BiasFieldCorrection(),
                                    name='n4biascorrect_t1')
    n4biasfieldcorrect_t1.inputs.dimension = 3
    n4biasfieldcorrect_t1.inputs.n_iterations = [50, 50, 30, 20]
    # n4biasfieldcorrect_t1.inputs.n_iterations = [5, 5, 3, 2]
    n4biasfieldcorrect_t1.inputs.convergence_threshold = 1e-6
    n4biasfieldcorrect_t1.inputs.bspline_fitting_distance = 300
    n4biasfieldcorrect_t1.inputs.output_image='T1_n4.nii.gz'
    n4biasfieldcorrect_t1.inputs.num_threads = num_threads

    #%% step-3c N4 T2
    n4biasfieldcorrect_t2 = pe.Node(interface=ants.N4BiasFieldCorrection(),
                                    name='n4biascorrect_t2')
    n4biasfieldcorrect_t2.inputs.dimension = 3
    n4biasfieldcorrect_t2.inputs.n_iterations = [50, 50, 30, 20]
    # n4biasfieldcorrect_t2.inputs.n_iterations = [5, 5, 3, 2]
    n4biasfieldcorrect_t2.inputs.convergence_threshold = 1e-6
    n4biasfieldcorrect_t2.inputs.bspline_fitting_distance = 300
    n4biasfieldcorrect_t2.inputs.output_image='T2_n4.nii.gz'
    n4biasfieldcorrect_t2.inputs.num_threads = num_threads


    #%% step-4: Register FST1, T1, T2 to FLAIR using FLIRT
    #step-4a flirt T1FS to FLAIR 
    t1fs_to_flair = pe.Node(interface=fsl.FLIRT(), name='t1fs_to_flair')
    t1fs_to_flair.inputs.cost = 'mutualinfo'
    t1fs_to_flair.inputs.dof = 12
    t1fs_to_flair.inputs.bins = 256
    t1fs_to_flair.inputs.searchr_x = [-25, 25]
    t1fs_to_flair.inputs.searchr_y = [-25, 25]
    t1fs_to_flair.inputs.searchr_z = [-25, 25]
    t1fs_to_flair.inputs.rigid2D=True
    t1fs_to_flair.inputs.interp='trilinear'
    t1fs_to_flair.inputs.out_matrix_file='FST1.mat'
    t1fs_to_flair.inputs.output_type = "NIFTI_GZ"
    t1fs_to_flair.inputs.out_file='FST1Warped.nii.gz'

    #%% step-4b flirt T1 to FLAIR 
    t1_to_flair = pe.Node(interface=fsl.FLIRT(), name='t1_to_flair')
    t1_to_flair.inputs.cost = 'mutualinfo'
    t1_to_flair.inputs.dof = 12
    t1_to_flair.inputs.bins = 256
    t1_to_flair.inputs.searchr_x = [-25, 25]
    t1_to_flair.inputs.searchr_y = [-25, 25]
    t1_to_flair.inputs.searchr_z = [-25, 25]
    t1_to_flair.inputs.rigid2D=True
    t1_to_flair.inputs.interp='trilinear'
    t1_to_flair.inputs.out_matrix_file='T1.mat'
    t1_to_flair.inputs.output_type = "NIFTI_GZ"
    t1_to_flair.inputs.out_file='T1Warped.nii.gz'
    
    #%% step-4c flirt T2 to FLAIR 
    t2_to_flair = pe.Node(interface=fsl.FLIRT(), name='t2_to_flair')
    t2_to_flair.inputs.cost = 'mutualinfo'
    t2_to_flair.inputs.dof = 12
    t2_to_flair.inputs.bins = 256
    t2_to_flair.inputs.searchr_x = [-25, 25]
    t2_to_flair.inputs.searchr_y = [-25, 25]
    t2_to_flair.inputs.searchr_z = [-25, 25]
    t2_to_flair.inputs.rigid2D=True
    t2_to_flair.inputs.interp='trilinear'
    t2_to_flair.inputs.out_matrix_file='T2.mat'
    t2_to_flair.inputs.output_type = "NIFTI_GZ"
    t2_to_flair.inputs.out_file='T2Warped.nii.gz'
    
    
    #%% step-4d flirt aseg to FLAIR
    aseg_to_flair = pe.Node(interface=fsl.FLIRT(), name='aseg_to_flair')
    aseg_to_flair.inputs.interp='nearestneighbour'
    aseg_to_flair.inputs.apply_xfm=True
    aseg_to_flair.inputs.output_type = "NIFTI_GZ"
    aseg_to_flair.inputs.out_file='FLAIRaseg.nii.gz'

    #%%step-5: Create brainmask in FLAIR space
    compute_mask_from_aseg = pe.Node(interface=util.Function(
            input_names=['in_file'], output_names=['out_file'],
            function=compute_mask), name='compute_mask_from_aseg')
    
    
    #%%step 14a maskout image FLAIR
    maskout_fl = pe.Node(interface=util.Function(
            input_names=['mask_file','image_file'],
            output_names=['maskoutfile'],
            function=maskout_image), name='maskout_fl')
    
    #step 14b maskout image T1warped
    maskout_t1w = pe.Node(interface=util.Function(
            input_names=['mask_file','image_file'],
            output_names=['maskoutfile'],
            function=maskout_image), name='maskout_t1w')
    
    #step 14c maskout image T2warped
    maskout_t2w = pe.Node(interface=util.Function(
            input_names=['mask_file','image_file'],
            output_names=['maskoutfile'],
            function=maskout_image), name='maskout_t2w')
    
    
    
    #%%step 15a normalize image FLAIR
    norm_fl = pe.Node(interface=util.Function(
            input_names=['mask_file','image_file'],
            output_names=['norm_outfile'],
            function=normalize_image), name='norm_fl')
    
    #step 15b normalize image T1warped
    norm_t1w = pe.Node(interface=util.Function(
            input_names=['mask_file','image_file'],
            output_names=['norm_outfile'],
            function=normalize_image), name='norm_t1w')
    
    #step 15c normalize image T2warped
    norm_t2w = pe.Node(interface=util.Function(
            input_names=['mask_file','image_file'],
            output_names=['norm_outfile'],
            function=normalize_image), name='norm_t2w')
    

    if not opp:  # if only preproc, then do not define deepmedic nodes

        #%% step-18 DeepMedic run

        create_t1_channel_config = pe.Node(interface=util.Function(
                input_names=['channel_file_path'],
                output_names=['channel_config_file'],
                function=create_deepmedic_channel_file
                ),name='create_t1_channel_config')

        create_t2_channel_config = pe.Node(interface=util.Function(
                input_names=['channel_file_path'],
                output_names=['channel_config_file'],
                function=create_deepmedic_channel_file
                ), name='create_t2_channel_config')
        

        create_flair_channel_config = pe.Node(interface=util.Function(
                input_names=['channel_file_path'],
                output_names=['channel_config_file'],
                function=create_deepmedic_channel_file
                ), name='create_flair_channel_config')

        create_roi_channel_config = pe.Node(interface=util.Function(
                input_names=['channel_file_path'],
                output_names=['channel_config_file'],
                function=create_deepmedic_channel_file
                ), name='create_roi_channel_config')
        
        #%%#
        # output predictions for each model in the ensemble
        create_pred_fl_channel_config = pe.Node(interface=util.Function(
                input_names=['channel_file_path', 'is_pred'],
                output_names=['channel_config_file'],
                function=create_deepmedic_channel_file
                ), name='create_pred_fl_channel_config')
        
        create_pred_fl_channel_config.inputs.is_pred = True

        create_pred_flt1_channel_config = pe.Node(interface=util.Function(
                input_names=['channel_file_path', 'is_pred'],
                output_names=['channel_config_file'],
                function=create_deepmedic_channel_file
                ), name='create_pred_flt1_channel_config')
        
        create_pred_flt1_channel_config.inputs.is_pred = True

        create_pred_flt2_channel_config = pe.Node(interface=util.Function(
                input_names=['channel_file_path', 'is_pred'],
                output_names=['channel_config_file'],
                function=create_deepmedic_channel_file
                ), name='create_pred_flt2_channel_config')
        
        create_pred_flt2_channel_config.inputs.is_pred = True

        create_pred_flt1t2_channel_config = pe.Node(interface=util.Function(
                input_names=['channel_file_path', 'is_pred'],
                output_names=['channel_config_file'],
                function=create_deepmedic_channel_file
                ), name='create_pred_flt1t2_channel_config')
        
        create_pred_flt1t2_channel_config.inputs.is_pred = True
        
        #%%#
        # test config for each model in the ensemble
        create_dm_test_fl_config = pe.Node(interface=util.Function(
                input_names=['model_type', 'flair_channel_file',
                             't1_channel_file','t2_channel_file',
                             'roi_channel_file','pred_channel_file'],
                             output_names=['test_config_file'],
                             function=create_deepmedic_config_file
                             ), name='create_dm_test_fl_config')
        
        create_dm_test_fl_config.inputs.model_type = 'fl'
        create_dm_test_fl_config.inputs.t1_channel_file = None
        create_dm_test_fl_config.inputs.t2_channel_file = None
        
        #%%#

        create_dm_test_flt1_config = pe.Node(interface=util.Function(
                input_names=['model_type', 'flair_channel_file',
                             't1_channel_file','t2_channel_file',
                             'roi_channel_file','pred_channel_file'],
                             output_names=['test_config_file'],
                             function=create_deepmedic_config_file
                             ), name='create_dm_test_flt1_config')
        
        create_dm_test_flt1_config.inputs.model_type = 'flt1'
        create_dm_test_flt1_config.inputs.t2_channel_file = None
        
        #%%#
        create_dm_test_flt2_config = pe.Node(interface=util.Function(
                input_names=['model_type', 'flair_channel_file',
                             't1_channel_file','t2_channel_file',
                             'roi_channel_file','pred_channel_file'],
                             output_names=['test_config_file'],
                             function=create_deepmedic_config_file
                             ), name='create_dm_test_flt2_config')
        
        create_dm_test_flt2_config.inputs.model_type = 'flt2'
        create_dm_test_flt2_config.inputs.t1_channel_file = None
        
        #%%#
        create_dm_test_flt1t2_config = pe.Node(interface=util.Function(
                input_names=['model_type', 'flair_channel_file',
                             't1_channel_file','t2_channel_file',
                             'roi_channel_file','pred_channel_file'],
                             output_names=['test_config_file'],
                             function=create_deepmedic_config_file
                             ), name='create_dm_test_flt1t2_config')
        
        create_dm_test_flt1t2_config.inputs.model_type = 'flt1t2'
        
        #%%#
        # deepmedic node for each model of the ensemble
        deepmedicrun_fl = pe.Node(interface=DeepMedic(),
                                  name='deepmedicrun_fl')
        
        deepmedicrun_fl.inputs.device = device
        deepmedicrun_fl.inputs.model_config_file = os.path.join(
                DM_MODEL_DIR, 'modelConfig_1ch.cfg')

        deepmedicrun_flt1 = pe.Node(interface=DeepMedic(),
                                    name='deepmedicrun_flt1')
        deepmedicrun_flt1.inputs.device = device
        deepmedicrun_flt1.inputs.model_config_file = os.path.join(
                DM_MODEL_DIR, 'modelConfig_2ch.cfg')

        deepmedicrun_flt2 = pe.Node(interface=DeepMedic(),
                                    name='deepmedicrun_flt2')
        deepmedicrun_flt2.inputs.device = device
        deepmedicrun_flt2.inputs.model_config_file = os.path.join(
                DM_MODEL_DIR, 'modelConfig_2ch.cfg')

        deepmedicrun_flt1t2 = pe.Node(interface=DeepMedic(),
                                      name='deepmedicrun_flt1t2')
        deepmedicrun_flt1t2.inputs.device = device
        deepmedicrun_flt1t2.inputs.model_config_file = os.path.join(
                DM_MODEL_DIR, 'modelConfig_3ch.cfg')
        

        if device=='cuda':
            #this will set the use_gpu attribute for the nodes input
            deepmedicrun_fl.inputs.use_gpu = True
            deepmedicrun_flt1.inputs.use_gpu = True
            deepmedicrun_flt2.inputs.use_gpu = True
            deepmedicrun_flt1t2.inputs.use_gpu = True
            
        #%%#
        # merge probability maps
        ensemble_segment = pe.Node(interface=util.Function(
                input_names=['fl_prob', 'flt1_prob', 'flt2_prob','flt1t2_prob'],
                output_names=['out_segmented_file'], function=merge_probmaps
                ), name='ensemble_segment')

        # global mask supra- and infra-tentorial labels for gm and wm
        global_mask = pe.Node(interface=util.Function(
                input_names=['in_file', 'include_superlist', 'fixed_id',
                             'map_pairs_list', 'out_file'],
                             output_names=['out_file'],
                             function=filter_labels), name='global_mask')
        
        global_mask.inputs.include_superlist = [[2, 41, 77],
                                                [3, 10, 11, 12, 13, 17, 18, 26,
                                                 42, 49, 50, 51, 52, 53, 54, 58],
                                                [7, 16, 46, 60, 28],
                                                [8, 47]]
        global_mask.inputs.fixed_id = None
        global_mask.inputs.map_pairs_list = [[2, 1], [3, 2], [7, 3], [8, 4]]
        global_mask.inputs.out_file = 'global_FLAIR.nii.gz'

        #%%save global volumes
        save_vols_glob = pe.Node(interface=util.Function(input_names=['wmh_file', 'global_file'], output_names=['json_file'],
                                                         function=vols_to_json), name='save_vols_glob')
        save_vols_lobes = pe.Node(interface=util.Function(input_names=['wmh_file', 'lobes_file'], output_names=['json_file'],
                                                         function=vols_to_json), name='save_vols_lobes')
        save_vols_shells = pe.Node(interface=util.Function(input_names=['wmh_file', 'shells_file'], output_names=['json_file'],
                                                         function=vols_to_json), name='save_vols_shells')
        save_vols_bulls = pe.Node(interface=util.Function(input_names=['wmh_file', 'bulls_file'], output_names=['json_file'],
                                                         function=vols_to_json), name='save_vols_bulls')
        save_vols_vter = pe.Node(interface=util.Function(input_names=['wmh_file', 'vter_file'], output_names=['json_file'],
                                                         function=vols_to_json), name='save_vols_vter')

        merge_json_basic = pe.Node(interface=util.Function(input_names=['MRI_GLOBAL'], output_names=['json_file'],
                                                         function=merge_jsons), name='merge_json_basic')
        merge_json_full = pe.Node(interface=util.Function(input_names=['MRI_GLOBAL', 'MRI_LOBE', 'MRI_SHELL', 'MRI_BULL', 'MRI_VASCTER'], output_names=['json_file'],
                                                         function=merge_jsons), name='merge_json_full')

        # snapshot for QC
        qc_snapshot = pe.Node(interface=util.Function(
                input_names=['wmh_file', 'flair_file', 'mask_file'],
                output_names=['out_fig'],
                function=plot_qc_fig), name='qc_snapshot')


        if not owmh:  # if only WMH then do not define Bullseye parcellation nodes

            # set freesurfer subjects_dir to scans_dir
            os.environ['SUBJECTS_DIR'] = scans_dir

            # lobar parcellation
            annot2label_lh = pe.Node(interface=Annot2Label(),
                                     name='annot2label_lh')
            annot2label_lh.inputs.hemi = 'lh'
            annot2label_lh.inputs.lobes = 'lobes'

            annot2label_rh = pe.Node(interface=Annot2Label(),
                                     name='annot2label_rh')
            annot2label_rh.inputs.hemi = 'rh'
            annot2label_rh.inputs.lobes = 'lobes'

            # aparc2aseg to map lobes into white matter volume
            aparc2aseg = pe.Node(interface=Aparc2Aseg(),
                                 name='aparc2aseg')
            aparc2aseg.inputs.annot = 'lobes'
            aparc2aseg.inputs.labelwm = True
            aparc2aseg.dmax = 1000
            aparc2aseg.inputs.rip = True
            aparc2aseg.inputs.hypo = True
            aparc2aseg.inputs.out_file = 'lobes+aseg.nii.gz'

            # group some lobes and discard others
            filter_lobes = pe.Node(
                interface=util.Function(
                        input_names=['in_file', 'include_superlist',
                                     'fixed_id', 'map_pairs_list'],
                                     output_names=['out_file'],
                                     function=filter_labels
                                     ), name='filter_lobes')
            
            # Here we include insula (3007, 4007) with frontal (3001, 4001)
            # We exclude the structure in the superior part spanning from anterior to posterior (3003, 4003)
            filter_lobes.inputs.include_superlist = [[3001, 3007],
                                                     [4001, 4007],
                                                     [3004],[4004],
                                                     [3005], [4005],
                                                     [3006],[4006]]#loblabels in WM
            filter_lobes.inputs.fixed_id = None
            # we give some informative label-ids
            filter_lobes.inputs.map_pairs_list = [[3001, 11], [4001, 21],
                                                  [3004, 12], [4004, 22],
                                                  [3005, 13], [4005, 23],
                                                  [3006, 14], [4006, 24]]

            # create ventricles and cortex masks
            ventricles_mask = pe.Node(
                interface=util.Function(
                        input_names=['in_file', 'include_superlist',
                                     'fixed_id', 'map_pairs_list'],
                                     output_names=['out_file'],
                                     function=filter_labels
                                     ), name='ventricles_mask')
            
            ventricles_mask.inputs.include_superlist = [[43, 4]]
            ventricles_mask.inputs.fixed_id = [1]
            ventricles_mask.inputs.map_pairs_list = None

            cortex_mask = pe.Node(interface=util.Function(
                    input_names=['in_file', 'include_superlist',
                                 'fixed_id', 'map_pairs_list'],
                                 output_names=['out_file'],
                                 function=filter_labels), name='cortex_mask')
            
            cortex_mask.inputs.include_superlist = [[1001, 2001, 1004, 2004,
                                                     1005, 2005, 1006, 2006]]  # lobar labels in cortex
            cortex_mask.inputs.fixed_id = [1]
            cortex_mask.inputs.map_pairs_list = None

            # create mask with basal ganglia + thalamus
            bgt_mask = pe.Node(interface=util.Function(
                    input_names=['in_file', 'include_superlist',
                                 'fixed_id', 'map_pairs_list'],
                                 output_names=['out_file'],
                                 function=filter_labels), name='bgt_mask')
            
            bgt_mask.inputs.include_superlist = [[10, 49, 11, 12, 50, 51,
                                                  26, 58, 13, 52]]  # basal ganglia + thalamus
            bgt_mask.inputs.fixed_id = [5]
            bgt_mask.inputs.map_pairs_list = None

            # create normalized distance map
            ndist_map = pe.Node(interface=util.Function(
                    input_names=['orig_file', 'dest_file'],
                    output_names=['out_file'],
                    function=norm_dist_map), name='ndist_map')
            

            # generate WM parcellations by filling the discarded lobes (3003, 4003) and unsegmented white matter (5001, 5002)
            gen_wmparc = pe.Node(interface=util.Function(
                    input_names=['incl_file', 'ndist_file',
                                 'label_file', 'incl_labels', 'verbose'],
                                 output_names=['out_file'],
                                 function=generate_wmparc), name='gen_wmparc')
            
            gen_wmparc.inputs.incl_labels = [3003, 4003, 5001, 5002]  # the labels that need to be 'filled'
            gen_wmparc.inputs.verbose = False

            # include bgt into wmparc to create the final lobar wmparc
            lobe_wmparc = pe.Node(interface=util.Function(
                    input_names=['in1_file', 'in2_file', 'out_file', 'intersect'],
                    output_names=['out_file'],
                    function=merge_labels), name='lobe_wmparc')
            
            lobe_wmparc.inputs.intersect = False
            lobe_wmparc.inputs.out_file = 'lobes_FSURF.nii.gz'

            # create depth shells using normalized distance maps
            depth_wmparc = pe.Node(interface=util.Function(
                    input_names=['ndist_file', 'n_shells', 'out_file', 'mask_file'],
                    output_names=['out_file'],
                    function=create_shells), name='depth_wmparc')
            
            depth_wmparc.inputs.n_shells = 4
            depth_wmparc.inputs.out_file = 'shells_FSURF.nii.gz'

            # final bullseye parcellation by intersecting depth and lobar parcellations
            bullseye_wmparc = pe.Node(interface=util.Function(
                    input_names=['in1_file', 'in2_file', 'out_file', 'intersect'],
                    output_names=['out_file'],
                    function=merge_labels), name='bullseye_wmparc')
            
            bullseye_wmparc.inputs.intersect = True
            bullseye_wmparc.inputs.out_file = 'bulls_FSURF.nii.gz'


            # nodes to bring regional parcellations to FLAIR space

            reorient2std_lobes = pe.Node(interface=Reorient2Std(),
                                         name='reorient2std_lobes')
            
            reorient2std_shells = pe.Node(interface=Reorient2Std(),
                                          name='reorient2std_shells')
            
            reorient2std_bulls = pe.Node(interface=Reorient2Std(),
                                         name='reorient2std_bulls')

            lobes_to_flair = pe.Node(interface=fsl.FLIRT(),
                                     name='lobes_to_flair')
            
            lobes_to_flair.inputs.interp = 'nearestneighbour'
            lobes_to_flair.inputs.apply_xfm = True
            lobes_to_flair.inputs.output_type = "NIFTI_GZ"
            lobes_to_flair.inputs.out_file = 'lobes_FLAIR.nii.gz'

            shells_to_flair = pe.Node(interface=fsl.FLIRT(),
                                      name='shells_to_flair')
            
            shells_to_flair.inputs.interp = 'nearestneighbour'
            shells_to_flair.inputs.apply_xfm = True
            shells_to_flair.inputs.output_type = "NIFTI_GZ"
            shells_to_flair.inputs.out_file = 'shells_FLAIR.nii.gz'

            bulls_to_flair = pe.Node(interface=fsl.FLIRT(),
                                     name='bulls_to_flair')
            
            bulls_to_flair.inputs.interp = 'nearestneighbour'
            bulls_to_flair.inputs.apply_xfm = True
            bulls_to_flair.inputs.output_type = "NIFTI_GZ"
            bulls_to_flair.inputs.out_file = 'bulls_FLAIR.nii.gz'
            
            #%%#
            # vascular territories

            # registration
            reg_vter = pe.Node(interface=ants.Registration(), name='reg_vter')
            reg_vter.inputs.fixed_image = os.path.join(
                    VASCTER_TEMPLATE_DIR,
                    'caa_flair_in_mni_template_smooth_brain_intres.nii.gz')
            
            reg_vter.inputs.output_transform_prefix = 'output_'
            reg_vter.inputs.initial_moving_transform_com = 1
            reg_vter.inputs.dimension = 3
            reg_vter.inputs.transforms = ['Rigid', 'Affine', 'SyN']
            n_steps = len(reg_vter.inputs.transforms)
            reg_vter.inputs.transform_parameters = [(0.1,), (0.1,),
                                                    (0.1, 3.0, 0.0)]
            reg_vter.inputs.number_of_iterations = [[500, 250, 125, 0],
                                                    [500, 250, 125, 0],
                                                    [70, 50, 30, 0]]
            reg_vter.inputs.winsorize_upper_quantile = 0.995
            reg_vter.inputs.winsorize_lower_quantile = 0.005
            reg_vter.inputs.float = True
            reg_vter.inputs.collapse_output_transforms = True
            # metrics
            reg_vter.inputs.metric = ['MI', 'MI', 'CC']
            reg_vter.inputs.metric_weight = [1, 1, 1]
            reg_vter.inputs.radius_or_number_of_bins = [32, 32, 4]
            reg_vter.inputs.sampling_strategy = ['Regular', 'Regular', None]
            reg_vter.inputs.sampling_percentage = [0.25, 0.25, None]
            reg_vter.inputs.convergence_threshold = [1.e-6] * n_steps
            reg_vter.inputs.convergence_window_size = [10] * n_steps
            reg_vter.inputs.smoothing_sigmas = [[3, 2, 1, 0]] * n_steps
            reg_vter.inputs.sigma_units = ['vox'] * n_steps
            reg_vter.inputs.shrink_factors = [[8, 4, 2, 1]] * n_steps
            reg_vter.inputs.write_composite_transform = True
            reg_vter.inputs.num_threads = num_threads

            # apply transforms

            vter_to_flair = pe.Node(interface=ants.ApplyTransforms(),
                                    name='vter_to_flair')
            vter_to_flair.inputs.dimension = 3
            vter_to_flair.inputs.input_image = os.path.join(
                    VASCTER_TEMPLATE_DIR,
                    'mni_vascular_territories.nii.gz')
            
            vter_to_flair.inputs.output_image = 'vascter_FLAIR.nii.gz'
            vter_to_flair.inputs.interpolation = 'NearestNeighbor'

            # create gray and white matter masks (for intersection with vascter)
            tissue_mask = pe.Node(interface=util.Function(
                    input_names=['in_file', 'include_superlist',
                                 'fixed_id', 'map_pairs_list', 'out_file'],
                                 output_names=['out_file'],
                                 function=filter_labels), name='tissue_mask')
            
            tissue_mask.inputs.include_superlist = [[2, 41, 77, 7, 
                                                     16, 46, 60, 28],
                                                    [3, 10, 11, 12, 13,
                                                     17,18, 26, 42, 49,
                                                     50,51, 52, 53, 54,
                                                     58, 8, 47]]
            tissue_mask.inputs.fixed_id = None
            tissue_mask.inputs.map_pairs_list = [[2, 1], [3, 2]]  # 1 -> WM, 2 -> GM
            tissue_mask.inputs.out_file = 'tissue_FLAIR.nii.gz'

            # intersection between vascter and tissues
            vter_tiss_parc = pe.Node(interface=util.Function(
                    input_names=['in1_file', 'in2_file', 'out_file',
                                 'intersect', 'num_digits'],
                                 output_names=['out_file'],
                                 function=merge_labels), name='vter_tiss_parc')
            
            vter_tiss_parc.inputs.intersect = True
            vter_tiss_parc.inputs.num_digits = 2
            vter_tiss_parc.inputs.out_file = 'vascter_tiss_FLAIR.nii.gz'



    #%% 18 collect outputs
    datasinkout = pe.Node(interface=DataSink(), name='datasinkout')
    datasinkout.inputs.parameterization=False

    ######## WORKFLOW CONNECTIONS #########

    #step 1a
    wmhsppwf.connect(inputnode        , 'subject_ids',      fileselector,'subject_id')
    wmhsppwf.connect(fileselector     , 'T1FS',             convert_t1_mgz, 'in_file')
    #step 1b
    wmhsppwf.connect(fileselector     , 'ASEG',             convert_aseg_mgz,'in_file')
    
    #step 1c
    wmhsppwf.connect(convert_aseg_mgz ,'out_file',          reorient2std_aseg,'in_file')
    wmhsppwf.connect(convert_t1_mgz   ,'out_file',          reorient2std_fst1,'in_file')
    
    #step 3a
    wmhsppwf.connect(fileselector     , 'FLAIR',           n4biasfieldcorrect_fl  , 'input_image')
    #step 3b
    wmhsppwf.connect(fileselector     , 'T1',           n4biasfieldcorrect_t1  , 'input_image')
    #step 3c
    wmhsppwf.connect(fileselector     , 'T2',           n4biasfieldcorrect_t2  , 'input_image')
    
    
    #step 4a
    wmhsppwf.connect(reorient2std_fst1       , 'out_file',         t1fs_to_flair  , 'in_file')
    wmhsppwf.connect(n4biasfieldcorrect_fl   , 'output_image',     t1fs_to_flair  , 'reference')
    #step 4b
    wmhsppwf.connect(n4biasfieldcorrect_t1   , 'output_image',     t1_to_flair  , 'in_file')
    wmhsppwf.connect(n4biasfieldcorrect_fl   , 'output_image',     t1_to_flair  , 'reference')
    #step 4c
    wmhsppwf.connect(n4biasfieldcorrect_t2   , 'output_image',     t2_to_flair  , 'in_file')
    wmhsppwf.connect(n4biasfieldcorrect_fl   , 'output_image',     t2_to_flair  , 'reference')
    #step 4d
    wmhsppwf.connect(reorient2std_aseg       ,'out_file',          aseg_to_flair  , 'in_file')
    wmhsppwf.connect(n4biasfieldcorrect_fl   , 'output_image',     aseg_to_flair  , 'reference')
    wmhsppwf.connect(t1fs_to_flair           , 'out_matrix_file',  aseg_to_flair  , 'in_matrix_file')
    
    #step 5
    wmhsppwf.connect(aseg_to_flair           , 'out_file',         compute_mask_from_aseg,'in_file')
    
    
    #step 14
    
    wmhsppwf.connect(compute_mask_from_aseg     , 'out_file',          maskout_fl, 'mask_file')
    wmhsppwf.connect(n4biasfieldcorrect_fl   , 'output_image',      maskout_fl, 'image_file')
        
    wmhsppwf.connect(compute_mask_from_aseg      , 'out_file',         maskout_t1w, 'mask_file')
    wmhsppwf.connect(t1_to_flair             , 'out_file',          maskout_t1w, 'image_file')
        
    wmhsppwf.connect(compute_mask_from_aseg      , 'out_file',         maskout_t2w, 'mask_file')
    wmhsppwf.connect(t2_to_flair             , 'out_file',          maskout_t2w, 'image_file')
    

    #step 15
    
    wmhsppwf.connect(compute_mask_from_aseg     , 'out_file',         norm_fl, 'mask_file')
    wmhsppwf.connect(maskout_fl                 , 'maskoutfile',      norm_fl, 'image_file')
        
    wmhsppwf.connect(compute_mask_from_aseg     , 'out_file',         norm_t1w, 'mask_file')
    wmhsppwf.connect(maskout_t1w                , 'maskoutfile',      norm_t1w, 'image_file')
        
    wmhsppwf.connect(compute_mask_from_aseg     , 'out_file',         norm_t2w, 'mask_file')
    wmhsppwf.connect(maskout_t2w                , 'maskoutfile',      norm_t2w, 'image_file')
    

    # outputs
    wmhsppwf.connect(inputnode               , 'subject_ids',       datasinkout, 'container')
    wmhsppwf.connect(inputnode               , 'outputdir',         datasinkout, 'base_directory')


    wmhsppwf.connect(norm_fl                 , 'norm_outfile',    datasinkout, '@flair')
    wmhsppwf.connect(norm_t1w                , 'norm_outfile',    datasinkout, '@t1w')
    wmhsppwf.connect(norm_t2w                , 'norm_outfile',    datasinkout, '@t2w')
    wmhsppwf.connect(compute_mask_from_aseg  , 'out_file',        datasinkout, '@brainmask')


    if not opp:  # if only preproc, then do not define deepmedic nodes

        #outputs for deepmedic
        wmhsppwf.connect(norm_fl                 , 'norm_outfile',       create_flair_channel_config, 'channel_file_path')
        wmhsppwf.connect(norm_t1w                , 'norm_outfile',       create_t1_channel_config, 'channel_file_path')
        wmhsppwf.connect(norm_t2w                , 'norm_outfile',       create_t2_channel_config, 'channel_file_path')
        wmhsppwf.connect(compute_mask_from_aseg  , 'out_file',           create_roi_channel_config, 'channel_file_path')

        #just a dummy input to create_pred_channel_config nodes
        wmhsppwf.connect(create_flair_channel_config, 'channel_config_file', create_pred_fl_channel_config, 'channel_file_path')
        wmhsppwf.connect(create_t1_channel_config  , 'channel_config_file', create_pred_flt1_channel_config, 'channel_file_path')
        wmhsppwf.connect(create_t2_channel_config  , 'channel_config_file', create_pred_flt2_channel_config, 'channel_file_path')
        wmhsppwf.connect(create_flair_channel_config  , 'channel_config_file', create_pred_flt1t2_channel_config, 'channel_file_path')


        wmhsppwf.connect(create_flair_channel_config  , 'channel_config_file', create_dm_test_fl_config, 'flair_channel_file')
        wmhsppwf.connect(create_flair_channel_config  , 'channel_config_file', create_dm_test_flt1_config, 'flair_channel_file')
        wmhsppwf.connect(create_flair_channel_config  , 'channel_config_file', create_dm_test_flt2_config, 'flair_channel_file')
        wmhsppwf.connect(create_flair_channel_config  , 'channel_config_file', create_dm_test_flt1t2_config, 'flair_channel_file')
        wmhsppwf.connect(create_t1_channel_config     , 'channel_config_file', create_dm_test_flt1_config, 't1_channel_file')
        wmhsppwf.connect(create_t1_channel_config     , 'channel_config_file', create_dm_test_flt1t2_config, 't1_channel_file')
        wmhsppwf.connect(create_t2_channel_config     , 'channel_config_file', create_dm_test_flt2_config, 't2_channel_file')
        wmhsppwf.connect(create_t2_channel_config     , 'channel_config_file', create_dm_test_flt1t2_config, 't2_channel_file')
        wmhsppwf.connect(create_roi_channel_config    , 'channel_config_file', create_dm_test_fl_config, 'roi_channel_file')
        wmhsppwf.connect(create_roi_channel_config    , 'channel_config_file', create_dm_test_flt1_config, 'roi_channel_file')
        wmhsppwf.connect(create_roi_channel_config    , 'channel_config_file', create_dm_test_flt2_config, 'roi_channel_file')
        wmhsppwf.connect(create_roi_channel_config    , 'channel_config_file', create_dm_test_flt1t2_config, 'roi_channel_file')
        wmhsppwf.connect(create_pred_fl_channel_config   , 'channel_config_file', create_dm_test_fl_config, 'pred_channel_file')
        wmhsppwf.connect(create_pred_flt1_channel_config   , 'channel_config_file', create_dm_test_flt1_config, 'pred_channel_file')
        wmhsppwf.connect(create_pred_flt2_channel_config   , 'channel_config_file', create_dm_test_flt2_config, 'pred_channel_file')
        wmhsppwf.connect(create_pred_flt1t2_channel_config   , 'channel_config_file', create_dm_test_flt1t2_config, 'pred_channel_file')

        wmhsppwf.connect(create_dm_test_fl_config        , 'test_config_file',    deepmedicrun_fl, 'test_config_file')
        wmhsppwf.connect(create_dm_test_flt1_config        , 'test_config_file',    deepmedicrun_flt1, 'test_config_file')
        wmhsppwf.connect(create_dm_test_flt2_config        , 'test_config_file',    deepmedicrun_flt2, 'test_config_file')
        wmhsppwf.connect(create_dm_test_flt1t2_config        , 'test_config_file',    deepmedicrun_flt1t2, 'test_config_file')


        wmhsppwf.connect(deepmedicrun_fl                 , 'out_probmap_file',  ensemble_segment,'fl_prob')
        wmhsppwf.connect(deepmedicrun_flt1                 , 'out_probmap_file',  ensemble_segment,'flt1_prob')
        wmhsppwf.connect(deepmedicrun_flt2                 , 'out_probmap_file',  ensemble_segment,'flt2_prob')
        wmhsppwf.connect(deepmedicrun_flt1t2                 , 'out_probmap_file',  ensemble_segment,'flt1t2_prob')

        # create global mask
        wmhsppwf.connect(aseg_to_flair, 'out_file', global_mask, 'in_file')

        # save global volumes
        wmhsppwf.connect(ensemble_segment, 'out_segmented_file', save_vols_glob, 'wmh_file')
        wmhsppwf.connect(global_mask, 'out_file', save_vols_glob, 'global_file')

        # QC snapshot
        wmhsppwf.connect(norm_fl, 'norm_outfile', qc_snapshot, 'flair_file')
        wmhsppwf.connect(compute_mask_from_aseg, 'out_file', qc_snapshot, 'mask_file')
        wmhsppwf.connect(ensemble_segment, 'out_segmented_file', qc_snapshot, 'wmh_file')

        # outputs to keep
        wmhsppwf.connect(global_mask, 'out_file', datasinkout, '@global_flair')
        wmhsppwf.connect(ensemble_segment, 'out_segmented_file',  datasinkout,'@segmentation')
        wmhsppwf.connect(qc_snapshot, 'out_fig', datasinkout, '@snapshot')


        if not owmh:  # if only WMH then do not create bullseye

            wmhsppwf.connect(fileselector, 'subject_id', annot2label_lh, 'subject')
            wmhsppwf.connect(fileselector, 'ANNOT_LH', annot2label_lh, 'in_annot')
            wmhsppwf.connect(fileselector, 'subject_id', annot2label_rh, 'subject')
            wmhsppwf.connect(fileselector, 'ANNOT_RH', annot2label_rh, 'in_annot')

            wmhsppwf.connect(annot2label_rh, 'out_annot_file', aparc2aseg, 'in_lobes_rh')
            wmhsppwf.connect(annot2label_lh, 'out_annot_file', aparc2aseg, 'in_lobes_lh')
            wmhsppwf.connect(fileselector, 'subject_id', aparc2aseg, 'subject')
            wmhsppwf.connect(aparc2aseg, 'out_file', filter_lobes, 'in_file')

            wmhsppwf.connect(aparc2aseg, 'out_file', ventricles_mask, 'in_file')
            wmhsppwf.connect(aparc2aseg, 'out_file', cortex_mask, 'in_file')
            wmhsppwf.connect(aparc2aseg, 'out_file', bgt_mask, 'in_file')

            wmhsppwf.connect(ventricles_mask, 'out_file', ndist_map, 'orig_file')
            wmhsppwf.connect(cortex_mask, 'out_file', ndist_map, 'dest_file')

            wmhsppwf.connect(aparc2aseg, 'out_file', gen_wmparc, 'incl_file')
            wmhsppwf.connect(ndist_map, 'out_file', gen_wmparc, 'ndist_file')
            wmhsppwf.connect(filter_lobes, 'out_file', gen_wmparc, 'label_file')

            wmhsppwf.connect(gen_wmparc, 'out_file', lobe_wmparc, 'in1_file')
            wmhsppwf.connect(bgt_mask, 'out_file', lobe_wmparc, 'in2_file')

            wmhsppwf.connect(ndist_map, 'out_file', depth_wmparc, 'ndist_file')
            wmhsppwf.connect(lobe_wmparc, 'out_file', depth_wmparc, 'mask_file')

            wmhsppwf.connect(lobe_wmparc, 'out_file', bullseye_wmparc, 'in1_file')
            wmhsppwf.connect(depth_wmparc, 'out_file', bullseye_wmparc, 'in2_file')

            wmhsppwf.connect(lobes_to_flair, 'out_file', save_vols_lobes, 'lobes_file')
            wmhsppwf.connect(ensemble_segment, 'out_segmented_file', save_vols_lobes, 'wmh_file')
            wmhsppwf.connect(shells_to_flair, 'out_file', save_vols_shells, 'shells_file')
            wmhsppwf.connect(ensemble_segment, 'out_segmented_file', save_vols_shells, 'wmh_file')
            wmhsppwf.connect(bulls_to_flair, 'out_file', save_vols_bulls, 'bulls_file')
            wmhsppwf.connect(ensemble_segment, 'out_segmented_file', save_vols_bulls, 'wmh_file')

            wmhsppwf.connect(lobe_wmparc, 'out_file', reorient2std_lobes, 'in_file')
            wmhsppwf.connect(reorient2std_lobes, 'out_file', lobes_to_flair, 'in_file')
            wmhsppwf.connect(n4biasfieldcorrect_fl, 'output_image', lobes_to_flair, 'reference')
            wmhsppwf.connect(t1fs_to_flair, 'out_matrix_file', lobes_to_flair, 'in_matrix_file')

            wmhsppwf.connect(depth_wmparc, 'out_file', reorient2std_shells, 'in_file')
            wmhsppwf.connect(reorient2std_shells, 'out_file', shells_to_flair, 'in_file')
            wmhsppwf.connect(n4biasfieldcorrect_fl, 'output_image', shells_to_flair, 'reference')
            wmhsppwf.connect(t1fs_to_flair, 'out_matrix_file', shells_to_flair, 'in_matrix_file')

            wmhsppwf.connect(bullseye_wmparc, 'out_file', reorient2std_bulls, 'in_file')
            wmhsppwf.connect(reorient2std_bulls, 'out_file', bulls_to_flair, 'in_file')
            wmhsppwf.connect(n4biasfieldcorrect_fl, 'output_image', bulls_to_flair, 'reference')
            wmhsppwf.connect(t1fs_to_flair, 'out_matrix_file', bulls_to_flair, 'in_matrix_file')

            # vascular territories

            wmhsppwf.connect(maskout_fl, 'maskoutfile', reg_vter, 'moving_image')

            wmhsppwf.connect(maskout_fl, 'maskoutfile', vter_to_flair, 'reference_image')
            wmhsppwf.connect(reg_vter, 'inverse_composite_transform', vter_to_flair, 'transforms')

            wmhsppwf.connect(aseg_to_flair, 'out_file', tissue_mask, 'in_file')
            wmhsppwf.connect(tissue_mask, 'out_file', vter_tiss_parc, 'in1_file')
            wmhsppwf.connect(vter_to_flair, 'output_image', vter_tiss_parc, 'in2_file')

            wmhsppwf.connect(vter_tiss_parc, 'out_file', save_vols_vter, 'vter_file')
            wmhsppwf.connect(ensemble_segment, 'out_segmented_file', save_vols_vter, 'wmh_file')

            # outputs in FLAIR space
            wmhsppwf.connect(lobes_to_flair, 'out_file', datasinkout, '@lobes_flair')
            wmhsppwf.connect(shells_to_flair, 'out_file', datasinkout, '@shells_flair')
            wmhsppwf.connect(bulls_to_flair, 'out_file', datasinkout, '@bullseye_flair')

            # outputs in T1FS space
            wmhsppwf.connect(lobe_wmparc, 'out_file', datasinkout, '@lobes_freesurf')
            wmhsppwf.connect(depth_wmparc, 'out_file', datasinkout, '@shells_freesurf')
            wmhsppwf.connect(bullseye_wmparc, 'out_file', datasinkout, '@bullseye_freesurf')

            # outputs vascular territories
            wmhsppwf.connect(vter_tiss_parc, 'out_file', datasinkout, '@vascter')

            # save vols
            wmhsppwf.connect(save_vols_glob, 'json_file', merge_json_full, 'MRI_GLOBAL')
            wmhsppwf.connect(save_vols_lobes, 'json_file',  merge_json_full,'MRI_LOBE')
            wmhsppwf.connect(save_vols_shells, 'json_file',  merge_json_full,'MRI_SHELL')
            wmhsppwf.connect(save_vols_bulls, 'json_file',  merge_json_full,'MRI_BULL')
            wmhsppwf.connect(save_vols_vter, 'json_file',  merge_json_full,'MRI_VASCTER')
            wmhsppwf.connect(merge_json_full, 'json_file', datasinkout, '@volumes')

        else:

            wmhsppwf.connect(save_vols_glob, 'json_file', merge_json_basic, 'MRI_GLOBAL')
            wmhsppwf.connect(merge_json_basic, 'json_file', datasinkout, '@volumes')

    return wmhsppwf
    

