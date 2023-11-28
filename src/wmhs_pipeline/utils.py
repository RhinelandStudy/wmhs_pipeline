from nipype.interfaces.base import (
    traits,
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    File
)

#from nipype.utils.filemanip import copyfile
import os

# os.environ['ANTSPATH'] = '/groups/mri-rhinelandstudy/software/ants2.3'
# os.environ['PATH'] = os.path.join(os.environ['ANTSPATH'], 'bin') + ':' + os.environ['PATH']
# os.environ['FSLDIR'] = '/groups/mri-rhinelandstudy/software/fsl/fsl6.0.0'
# os.environ['PATH'] = os.path.join(os.environ['FSLDIR'], 'bin') + ':' + os.environ['PATH']
# os.environ['FREESURFER_HOME'] = '/groups/mri-rhinelandstudy/software/freesurfer/freesurfer6.0.0'
# os.environ['PATH'] = os.path.join(os.environ['FREESURFER_HOME'], 'bin') + ':' + os.environ['PATH']
# os.environ['SUBJECTS_DIR'] = os.path.join(os.environ['FREESURFER_HOME'], 'subjects')

################################################
########## WMH SEGMENTATION UTILS ##############
################################################

def get_first_file(file_list):
    return file_list[0]

def compute_mask(in_file):
    """
    create a brainmask from the input image
    """
    import nibabel as nib
    import numpy as np
    import os
    from scipy.ndimage.morphology import binary_dilation
    
    dilat_rad = [5, 5, 5]
    struct = np.ones(dilat_rad, dtype=np.bool)
    img_nib = nib.load(in_file)
    img = img_nib.get_data()
    
    mask = np.zeros(img.shape, dtype=np.bool)
    mask[img > 0] = True
    mask = binary_dilation(mask, struct)
    mask_nib = nib.Nifti1Image(mask, img_nib.affine, img_nib.header)
    nib.save(mask_nib, 'brainmask_FLAIR.nii.gz')
    return os.path.abspath('brainmask_FLAIR.nii.gz')

    
def convert_mgz(in_file):
    """
    convert mgz to nii.gz 
    """
    from nipype.interfaces.freesurfer.preprocess import MRIConvert
    import os.path as op

    fname, ext = op.splitext(op.basename(in_file))
    if ext == ".gz":
        return in_file
        #fname, ext2 = op.splitext(fname)
        #ext = ext2 + ext
    else:
        mc = MRIConvert()
        mc.inputs.in_file = in_file
        mc.inputs.out_type = 'niigz'
        mc.inputs.out_file = fname + '.nii.gz'
        res = mc.run()
        outputs = res.outputs.get()
        return outputs['out_file']


def maskout_image(mask_file, image_file):
    import nibabel as nib
    import numpy as np
    import os
    
    # read mask
    mask_nib = nib.load(mask_file)
    mask = mask_nib.get_data().astype(np.bool)
    
    # maskout images
    img_nib = nib.load(image_file)
    img = img_nib.get_data()
    img[~mask] = img.min()

    aux = nib.Nifti1Image(img, img_nib.affine, img_nib.header)
    fname, ext = os.path.splitext(os.path.basename(image_file))
    maskoutname = fname.replace('.nii','') + '_maskout.nii.gz'
    maskoutfile = os.path.join(os.getcwd(),maskoutname)
    nib.save(aux, maskoutfile)

    return maskoutfile

def normalize_image(mask_file, image_file):
    import nibabel as nib
    import numpy as np
    import os
    from scipy.stats.mstats import zscore
    
    # read mask
    mask_nib = nib.load(mask_file)
    mask = mask_nib.get_data().astype(np.bool)
    
    # normalize
    img_nib = nib.load(image_file)
    img = img_nib.get_data()
    
    img2 = np.zeros(img.shape, dtype=np.float32)
    img2[mask] = zscore(img[mask].astype(np.float32))
    aux = nib.Nifti1Image(img2, img_nib.affine, img_nib.header.set_data_dtype(np.float32))
    fname,ext = os.path.splitext(os.path.basename(image_file))
    normname = fname.replace('.nii','') + '_norm.nii.gz'
    normoutfile = os.path.join(os.getcwd(),normname)
    nib.save(aux, normoutfile)
    
    return normoutfile
    
    
def merge_probmaps(**probmaps_dict):

    import nibabel as nib
    import numpy as np
    import os
    from copy import copy

    ref_nib = None
    avg = None  # average probmap

    # average probmaps in across folders

    for key,probmap_file in probmaps_dict.items():

        prob_nib = nib.load(probmap_file)
        prob = prob_nib.get_data()

        if avg is None:
            avg = copy(prob) / float(len(probmaps_dict))
            ref_nib = copy(prob_nib)
        else:
            assert prob_nib.header.get_data_shape() == ref_nib.header.get_data_shape(), 'non-equal probmap shapes'
            avg += prob / float(len(probmaps_dict))

    # threshold average
    seg = np.zeros(avg.shape, dtype=np.uint8)
    seg[avg > 0.5] = 1

    # save
    seg_nib = nib.Nifti1Image(seg, ref_nib.affine, ref_nib.header)
    seg_nib.set_data_dtype(np.uint8)
    nib.save(seg_nib, 'WMH_FLAIR.nii.gz')

    return os.path.abspath('WMH_FLAIR.nii.gz')


def plot_qc_fig(wmh_file, flair_file, mask_file, out_file=None, N=35):

    import matplotlib.image as mpimg
    import nibabel as nib
    import numpy as np
    from copy import copy
    import os

    wmh_nib = nib.load(wmh_file)
    mask_nib = nib.load(mask_file)

    wmh_np = wmh_nib.get_data()
    mask_np = mask_nib.get_data()

    wmh_dist = wmh_np.sum(0).sum(0)  # sum of intensities in every slice of 'z' direction
    mask_dist = mask_np.sum(0).sum(0)

    total_wmh = float(wmh_dist.sum())
    if total_wmh > 0.:
        wmh_cum = wmh_dist.cumsum().astype(float) / total_wmh
        wmh_bounds = [np.argmax(wmh_cum > 1e-9), np.argmax(wmh_cum > (1. - 1e-9))] # if there's WMH, set bounds to ini and end of WMH mask in z-axis
    else:
        wmh_bounds = [len(wmh_dist) / 2] * 2  # if theres no WMH, then set bounds to half of z-axis

    total_mask = float(mask_dist.sum())
    assert total_mask > 0., 'mask image is zero everywhere'
    mask_cum = mask_dist.cumsum().astype(float) / total_mask
    mask_bounds = [np.argmax(mask_cum > 0.05), np.argmax(mask_cum > 0.95)]  # image bounds to 5% and 95% mask percentiles in z-axis

    # print('wmh bounds: (%d, %d)' % (wmh_bounds[0], wmh_bounds[1]))
    # print('mask bounds: (%d, %d)' % (mask_bounds[0], mask_bounds[1]))

    # set inferior and superior bounds half-way between the wmh and mask bounds (only if mask bounds extend beyond wmh bounds)
    inf_bound = np.mean([wmh_bounds[0], min(mask_bounds[0], wmh_bounds[0])]).astype(int)
    sup_bound = np.mean([wmh_bounds[1], max(mask_bounds[1], wmh_bounds[1])]).astype(int)

    # print('inf %d, sup %d' % (inf_bound, sup_bound))

    # enlarge bounds in case less than N slides
    residue = N - (sup_bound - inf_bound)
    margin = np.ceil(max(0., residue) / 2.).astype(int)
    img_bounds = [inf_bound - margin, sup_bound + margin]

    # print('img bounds: %d, %d' % (img_bounds[0], img_bounds[1]))

    # plot the image

    flair_np = nib.load(flair_file).get_data()
    flair_np = (flair_np - flair_np.min()) / (flair_np.max() - flair_np.min())  # normalize

    h = flair_np.shape[1]  # flair image dims
    w = flair_np.shape[0]

    height = h * 2  # QC image dims
    width = w * N

    qc_img = np.empty((height, width, 3))

    for i, idx in enumerate(np.linspace(img_bounds[0], img_bounds[1], N, dtype=int)):

        # flair slice
        slice_fl = flair_np[:, ::-1, idx].T

        # flair slices with overlaid wmh segment
        mask_idx = wmh_np[:, ::-1, idx].T
        slice_wmh1 = copy(slice_fl)
        slice_wmh1[mask_idx > 0] = 1.
        slice_wmh0 = copy(slice_fl)
        slice_wmh0[mask_idx > 0] = 0.

        # naked flair
        qc_img[:h, i*w:(i+1)*w, :] = copy(np.repeat(slice_fl[..., None], 3, 2))

        # flair with wmh overlay (activate red channel + desactivate green and blue channels in wmh pixels)
        qc_img[h:, i*w:(i+1)*w, :] = copy(np.concatenate((slice_wmh1[..., None], np.repeat(slice_wmh0[..., None], 2, 2)), axis=2))


    if out_file is None:
        out_file = 'QCsnapshot.png'

    mpimg.imsave(out_file, qc_img)

    return os.path.abspath(out_file)


def vols_to_json(wmh_file, global_file=None, lobes_file=None, shells_file=None, bulls_file=None, vter_file=None):
    """
    computes global wmh volume and stores in json file
    """
    import nibabel as nib
    import numpy as np
    from json import dump
    import os

    wmh_nib = nib.load(wmh_file)
    wmh = wmh_nib.get_data().astype(np.bool)
    mm3_vol = np.prod(np.asarray(wmh_nib.header.get_zooms()).astype(float))

    if global_file is not None:

        parc = nib.load(global_file).get_data()

        labels_list = [1, 2, 3, 4]
        wmh_names_list = ['MRI_GLOBAL_SUPRATENTORIAL_WMH_DM_R1', 'MRI_GLOBAL_SUPRATENTORIAL_GMH_DM_R1',
                          'MRI_GLOBAL_INFRATENTORIAL_WMH_DM_R1', 'MRI_GLOBAL_INFRATENTORIAL_GMH_DM_R1']
        reg_names_list = ['MRI_GLOBAL_SUPRATENTORIAL_WM_DM_R1', 'MRI_GLOBAL_SUPRATENTORIAL_GM_DM_R1',
                          'MRI_GLOBAL_INFRATENTORIAL_WM_DM_R1', 'MRI_GLOBAL_INFRATENTORIAL_GM_DM_R1']

        json_file = os.path.abspath('vols_global.json')

    elif lobes_file is not None:

        parc = nib.load(lobes_file).get_data()

        labels_list = [11, 21, 12, 22, 13, 23, 14, 24, 5]
        wmh_names_list = ['MRI_LOBE_FRONTAL_WMH_LH_DM_R1', 'MRI_LOBE_FRONTAL_WMH_RH_DM_R1',
                          'MRI_LOBE_OCCIPITAL_WMH_LH_DM_R1', 'MRI_LOBE_OCCIPITAL_WMH_RH_DM_R1',
                          'MRI_LOBE_TEMPORAL_WMH_LH_DM_R1', 'MRI_LOBE_TEMPORAL_WMH_RH_DM_R1',
                          'MRI_LOBE_PARIETAL_WMH_LH_DM_R1', 'MRI_LOBE_PARIETAL_WMH_RH_DM_R1',
                          'MRI_LOBE_BASGANGTHAL_GMH_DM_R1']
        reg_names_list = ['MRI_LOBE_FRONTAL_LH_DM_R1', 'MRI_LOBE_FRONTAL_RH_DM_R1',
                          'MRI_LOBE_OCCIPITAL_LH_DM_R1', 'MRI_LOBE_OCCIPITAL_RH_DM_R1',
                          'MRI_LOBE_TEMPORAL_LH_DM_R1', 'MRI_LOBE_TEMPORAL_RH_DM_R1',
                          'MRI_LOBE_PARIETAL_LH_DM_R1', 'MRI_LOBE_PARIETAL_RH_DM_R1',
                          'MRI_LOBE_BASGANGTHAL_DM_R1']

        json_file = os.path.abspath('vols_lobe.json')

    elif shells_file is not None:

        parc = nib.load(shells_file).get_data()

        labels_list = [1, 2, 3, 4]
        wmh_names_list = ['MRI_SHELL_PERIVENT_WGMH_DM_R1', 'MRI_SHELL_MEDPERIVENT_WGMH_DM_R1', 'MRI_SHELL_MEDDEEP_WGMH_DM_R1', 'MRI_SHELL_DEEP_WGMH_DM_R1']
        reg_names_list = ['MRI_SHELL_PERIVENT_DM_R1', 'MRI_SHELL_MEDPERIVENT_DM_R1', 'MRI_SHELL_MEDDEEP_DM_R1', 'MRI_SHELL_DEEP_DM_R1']

        json_file = os.path.abspath('vols_shell.json')

    elif bulls_file is not None:

        parc = nib.load(bulls_file).get_data()

        labels_list = [111, 211, 121, 221, 131, 231, 141, 241, 51,
                       112, 212, 122, 222, 132, 232, 142, 242, 52,
                       113, 213, 123, 223, 133, 233, 143, 243, 53,
                       114, 214, 124, 224, 134, 234, 144, 244, 54]

        wmh_names_list = ['MRI_BULL_FRONTAL_PERIVENT_WMH_LH_DM_R1', 'MRI_BULL_FRONTAL_PERIVENT_WMH_RH_DM_R1',
                          'MRI_BULL_OCCIPITAL_PERIVENT_WMH_LH_DM_R1', 'MRI_BULL_OCCIPITAL_PERIVENT_WMH_RH_DM_R1',
                          'MRI_BULL_TEMPORAL_PERIVENT_WMH_LH_DM_R1', 'MRI_BULL_TEMPORAL_PERIVENT_WMH_RH_DM_R1',
                          'MRI_BULL_PARIETAL_PERIVENT_WMH_LH_DM_R1', 'MRI_BULL_PARIETAL_PERIVENT_WMH_RH_DM_R1',
                          'MRI_BULL_BASGANGTHAL_PERIVENT_GMH_DM_R1',

                          'MRI_BULL_FRONTAL_MEDPERIVENT_WMH_LH_DM_R1', 'MRI_BULL_FRONTAL_MEDPERIVENT_WMH_RH_DM_R1',
                          'MRI_BULL_OCCIPITAL_MEDPERIVENT_WMH_LH_DM_R1', 'MRI_BULL_OCCIPITAL_MEDPERIVENT_WMH_RH_DM_R1',
                          'MRI_BULL_TEMPORAL_MEDPERIVENT_WMH_LH_DM_R1', 'MRI_BULL_TEMPORAL_MEDPERIVENT_WMH_RH_DM_R1',
                          'MRI_BULL_PARIETAL_MEDPERIVENT_WMH_LH_DM_R1', 'MRI_BULL_PARIETAL_MEDPERIVENT_WMH_RH_DM_R1',
                          'MRI_BULL_BASGANGTHAL_MEDPERIVENT_GMH_DM_R1',

                          'MRI_BULL_FRONTAL_MEDDEEP_WMH_LH_DM_R1', 'MRI_BULL_FRONTAL_MEDDEEP_WMH_RH_DM_R1',
                          'MRI_BULL_OCCIPITAL_MEDDEEP_WMH_LH_DM_R1', 'MRI_BULL_OCCIPITAL_MEDDEEP_WMH_RH_DM_R1',
                          'MRI_BULL_TEMPORAL_MEDDEEP_WMH_LH_DM_R1', 'MRI_BULL_TEMPORAL_MEDDEEP_WMH_RH_DM_R1',
                          'MRI_BULL_PARIETAL_MEDDEEP_WMH_LH_DM_R1', 'MRI_BULL_PARIETAL_MEDDEEP_WMH_RH_DM_R1',
                          'MRI_BULL_BASGANGTHAL_MEDDEEP_GMH_DM_R1',

                          'MRI_BULL_FRONTAL_DEEP_WMH_LH_DM_R1', 'MRI_BULL_FRONTAL_DEEP_WMH_RH_DM_R1',
                          'MRI_BULL_OCCIPITAL_DEEP_WMH_LH_DM_R1', 'MRI_BULL_OCCIPITAL_DEEP_WMH_RH_DM_R1',
                          'MRI_BULL_TEMPORAL_DEEP_WMH_LH_DM_R1', 'MRI_BULL_TEMPORAL_DEEP_WMH_RH_DM_R1',
                          'MRI_BULL_PARIETAL_DEEP_WMH_LH_DM_R1', 'MRI_BULL_PARIETAL_DEEP_WMH_RH_DM_R1',
                          'MRI_BULL_BASGANGTHAL_DEEP_GMH_DM_R1']

        reg_names_list = ['MRI_BULL_FRONTAL_PERIVENT_LH_DM_R1', 'MRI_BULL_FRONTAL_PERIVENT_RH_DM_R1',
                          'MRI_BULL_OCCIPITAL_PERIVENT_LH_DM_R1', 'MRI_BULL_OCCIPITAL_PERIVENT_RH_DM_R1',
                          'MRI_BULL_TEMPORAL_PERIVENT_LH_DM_R1', 'MRI_BULL_TEMPORAL_PERIVENT_RH_DM_R1',
                          'MRI_BULL_PARIETAL_PERIVENT_LH_DM_R1', 'MRI_BULL_PARIETAL_PERIVENT_RH_DM_R1',
                          'MRI_BULL_BASGANGTHAL_PERIVENT_DM_R1',

                          'MRI_BULL_FRONTAL_MEDPERIVENT_LH_DM_R1', 'MRI_BULL_FRONTAL_MEDPERIVENT_RH_DM_R1',
                          'MRI_BULL_OCCIPITAL_MEDPERIVENT_LH_DM_R1', 'MRI_BULL_OCCIPITAL_MEDPERIVENT_RH_DM_R1',
                          'MRI_BULL_TEMPORAL_MEDPERIVENT_LH_DM_R1', 'MRI_BULL_TEMPORAL_MEDPERIVENT_RH_DM_R1',
                          'MRI_BULL_PARIETAL_MEDPERIVENT_LH_DM_R1', 'MRI_BULL_PARIETAL_MEDPERIVENT_RH_DM_R1',
                          'MRI_BULL_BASGANGTHAL_MEDPERIVENT_DM_R1',

                          'MRI_BULL_FRONTAL_MEDDEEP_LH_DM_R1', 'MRI_BULL_FRONTAL_MEDDEEP_RH_DM_R1',
                          'MRI_BULL_OCCIPITAL_MEDDEEP_LH_DM_R1', 'MRI_BULL_OCCIPITAL_MEDDEEP_RH_DM_R1',
                          'MRI_BULL_TEMPORAL_MEDDEEP_LH_DM_R1', 'MRI_BULL_TEMPORAL_MEDDEEP_RH_DM_R1',
                          'MRI_BULL_PARIETAL_MEDDEEP_LH_DM_R1', 'MRI_BULL_PARIETAL_MEDDEEP_RH_DM_R1',
                          'MRI_BULL_BASGANGTHAL_MEDDEEP_DM_R1',

                          'MRI_BULL_FRONTAL_DEEP_LH_DM_R1', 'MRI_BULL_FRONTAL_DEEP_RH_DM_R1',
                          'MRI_BULL_OCCIPITAL_DEEP_LH_DM_R1', 'MRI_BULL_OCCIPITAL_DEEP_RH_DM_R1',
                          'MRI_BULL_TEMPORAL_DEEP_LH_DM_R1', 'MRI_BULL_TEMPORAL_DEEP_RH_DM_R1',
                          'MRI_BULL_PARIETAL_DEEP_LH_DM_R1', 'MRI_BULL_PARIETAL_DEEP_RH_DM_R1',
                          'MRI_BULL_BASGANGTHAL_DEEP_DM_R1']

        json_file = os.path.abspath('vols_bull.json')

    elif vter_file is not None:

        parc = nib.load(vter_file).get_data()

        labels_list = [ 102, 104, 106, 108, 110, 112, 114, 116, 118, 120,
                        202, 204, 206, 210, 212, 214, 216, 220]
        wmh_names_list = ['MRI_VASCTER_ACA_WMH_RH_DM_R1', 'MRI_VASCTER_MCA_WMH_RH_DM_R1', 'MRI_VASCTER_PCA_WMH_RH_DM_R1', 'MRI_VASCTER_PONS_WMH_RH_DM_R1', 'MRI_VASCTER_CEREBELLUM_WMH_RH_DM_R1',
                          'MRI_VASCTER_ACA_WMH_LH_DM_R1', 'MRI_VASCTER_MCA_WMH_LH_DM_R1', 'MRI_VASCTER_PCA_WMH_LH_DM_R1', 'MRI_VASCTER_PONS_WMH_LH_DM_R1', 'MRI_VASCTER_CEREBELLUM_WMH_LH_DM_R1',
                          'MRI_VASCTER_ACA_GMH_RH_DM_R1', 'MRI_VASCTER_MCA_GMH_RH_DM_R1', 'MRI_VASCTER_PCA_GMH_RH_DM_R1', 'MRI_VASCTER_CEREBELLUM_GMH_RH_DM_R1',
                          'MRI_VASCTER_ACA_GMH_LH_DM_R1', 'MRI_VASCTER_MCA_GMH_LH_DM_R1', 'MRI_VASCTER_PCA_GMH_LH_DM_R1', 'MRI_VASCTER_CEREBELLUM_GMH_LH_DM_R1']
        reg_names_list = ['MRI_VASCTER_ACA_WM_RH_DM_R1', 'MRI_VASCTER_MCA_WM_RH_DM_R1', 'MRI_VASCTER_PCA_WM_RH_DM_R1', 'MRI_VASCTER_PONS_WM_RH_DM_R1', 'MRI_VASCTER_CEREBELLUM_WM_RH_DM_R1',
                          'MRI_VASCTER_ACA_WM_LH_DM_R1', 'MRI_VASCTER_MCA_WM_LH_DM_R1', 'MRI_VASCTER_PCA_WM_LH_DM_R1', 'MRI_VASCTER_PONS_WM_LH_DM_R1', 'MRI_VASCTER_CEREBELLUM_WM_LH_DM_R1',
                          'MRI_VASCTER_ACA_GM_RH_DM_R1', 'MRI_VASCTER_MCA_GM_RH_DM_R1', 'MRI_VASCTER_PCA_GM_RH_DM_R1', 'MRI_VASCTER_CEREBELLUM_GM_RH_DM_R1',
                          'MRI_VASCTER_ACA_GM_LH_DM_R1', 'MRI_VASCTER_MCA_GM_LH_DM_R1', 'MRI_VASCTER_PCA_GM_LH_DM_R1', 'MRI_VASCTER_CEREBELLUM_GM_LH_DM_R1']

        json_file = os.path.abspath('vols_vascter.json')

    # Volumes
    dict = {}
    for label, wmh_name, reg_name in zip(labels_list, wmh_names_list, reg_names_list):
        mask = parc == label
        dict.update({wmh_name: float(wmh[mask].sum()) * mm3_vol})
        dict.update({reg_name: float(mask.sum()) * mm3_vol})

    with open(json_file, 'w') as f:
        dump(dict, f)

    return json_file


def merge_jsons(**jsons_dict):

    import os, json

    dict = {}

    for key, json_file in jsons_dict.items():

        with open(json_file) as f:

            s = f.read()
            d = json.loads(s)

            dict.update({key: d})

    out_file = 'vols_wmh.json'
    with open(out_file, 'w') as f:
        json.dump(dict, f)

    return os.path.abspath(out_file)



def create_deepmedic_channel_file(channel_file_path, is_pred=False):
    """
    create channel configuration file required by deepMedicRun
    """

    import os

    if is_pred:
        channel_file_path = 'pred'

    channel_config_file = 'testChannel.cfg'
    with open(channel_config_file, 'w') as fid:
        fid.write(channel_file_path + '\n')

    return os.path.abspath(channel_config_file)


def create_deepmedic_config_file(model_type, flair_channel_file, t1_channel_file, t2_channel_file, roi_channel_file, pred_channel_file):
    import os

    from wmhs_pipeline.configoptions import DM_MODEL_DIR
    # from configoptions import DM_MODEL_DIR

    test_config_file = 'testConfig.cfg'

    # this workaround to set the output path to the deepmedic run folder
    folder_for_output = os.path.join(os.path.abspath(os.path.join(os.getcwd(), os.pardir)), 'deepmedicrun_' + model_type)
    model_file_path = os.path.join(DM_MODEL_DIR, 'best_mod_compp_' + model_type + '.ckpt')

    channels_aux = '['
    for chan in [flair_channel_file, t1_channel_file, t2_channel_file]:
        if chan is not None:
            channels_aux += '"' + chan + '", '
    channels = channels_aux[:-2] + ']'

    with open(test_config_file, 'w') as fid:
        fid.write('sessionName = "deepmed_v1"' + '\n\n')
        fid.write('folderForOutput = "' + folder_for_output + '"\n\n')
        fid.write('cnnModelFilePath = "' + model_file_path + '"\n\n')
        fid.write('channels = ' + channels + '\n\n')
        fid.write('namesForPredictionsPerCase = "' + pred_channel_file + '"\n\n')
        fid.write('roiMasks = "' + roi_channel_file + '"\n\n')
        fid.write('batchsize = 10\n\n')
        fid.write('saveSegmentation = True\n')
        fid.write('saveProbMapsForEachClass = [True, True, True, True, True]\n')
        fid.write('saveIndividualFms = False\n')
        fid.write('saveAllFmsIn4DimImage = False\n')
        # fid.write('minMaxIndicesOfFmsToSaveFromEachLayerOfNormalPathway = []\n')
        # fid.write('minMaxIndicesOfFmsToSaveFromEachLayerOfSubsampledPathway = [[],[],[],[],[],[],[],[]]\n')
        # fid.write('minMaxIndicesOfFmsToSaveFromEachLayerOfFullyConnectedPathway = [[],[0,150],[]]\n')
        fid.write('padInputImagesBool = True\n')

    return os.path.abspath(test_config_file)


class DeepMedicInputSpec(CommandLineInputSpec):
    """
    interface for DeepMedic
    """
    model_config_file = File(exists=True, desc='deepMedic model config file.', argstr='-model %s', position=0,
                             mandatory=True)
    test_config_file = File(exists=True, desc='deepMedic test config file.', argstr='-test %s', position=1,
                            mandatory=True)
    # load_saved_model  = File(exists=True, desc='deepMedic saved model file.',   argstr='-load %s',  position=2, mandatory=True)
    device = traits.String(desc='device name', argstr='-dev %s', position=3, mandatory=True)
    use_gpu = traits.Bool(desc='set  the flag to use gpu')


class DeepMedicOutputSpec(TraitedSpec):
    # out_segmented_file = File(exists=True, desc='Output files from deepMedicRun')
    out_probmap_file = File(exists=True, desc='Output probability map from deepMedicRun')


class DeepMedic(CommandLine):
    _cmd = 'deepMedicRun'
    # _cmd = '/home/sanromag/CODE/external/deepmedic/deepMedicRun'
    input_spec = DeepMedicInputSpec
    output_spec = DeepMedicOutputSpec

    def __init__(self, **inputs):
        return super(DeepMedic, self).__init__(**inputs)

    def _format_arg(self, name, spec, value):
        if (name == 'load_saved_model'):
            # remove the .index extension here, the file will be existing file
            return spec.argstr % (self.inputs.load_saved_model.replace('.index', ''))

        return super(DeepMedic,
                     self)._format_arg(name, spec, value)

    def _run_interface(self, runtime):
        runtime = super(DeepMedic, self)._run_interface(runtime)

        # I need to comment this out to avoid a persistent error after deepmedic call finishes
        # if runtime.stderr:
        #     self.raise_exception(runtime)

        return runtime

    def _list_outputs(self):
        # outputs = self._outputs().get()
        # outputs['out_segmented_file'] = os.path.abspath(
        #     os.getcwd() + '/predictions/deepmed_v1/predictions/pred_Segm.nii.gz')
        outputs = self.output_spec().get()
        # outputs['out_segmented_file'] = os.path.abspath('predictions/deepmed_v1/predictions/pred_Segm.nii.gz')
        outputs['out_probmap_file'] = os.path.abspath('predictions/deepmed_v1/predictions/pred_ProbMapClass1.nii.gz')

        return outputs


################################################
####### BULLSEYE PARCELLATION UTILS ############
################################################


def filter_labels(in_file, include_superlist, fixed_id=None, map_pairs_list=None, out_file='filtered.nii.gz'):
    """filters-out labels not in the include-superset. Merges labels within superset. Transforms label-ids according to mappings (or fixed id)"""
    import nibabel as nib
    import numpy as np
    import os

    # read label file and create output
    in_nib = nib.load(in_file)
    in0 = in_nib.get_data()
    out0 = np.zeros(in0.shape, dtype=in0.dtype)

    # for each group of labels in subset assign them the same id (either 1st in subset or fixed-id, in case given)
    for labels_list in include_superlist:
        for label in labels_list:
            value = labels_list[0]
            if fixed_id is not None: value = fixed_id[0]
            out0[in0 == label] = value

    # transform label-ids in case mapping is specified
    if map_pairs_list is not None:
        out1 = np.copy(out0)
        for map_pair in map_pairs_list:
            out1[out0 == map_pair[0]] = map_pair[1]

    # save output
    out_final = out0 if not map_pairs_list else out1
    out_nib = nib.Nifti1Image(out_final, in_nib.affine, in_nib.header)
    nib.save(out_nib, out_file)

    return os.path.abspath(out_file)


def norm_dist_map(orig_file, dest_file):
    """compute normalized distance map given an origin and destination masks, resp."""
    import os
    import nibabel as nib
    import numpy as np
    from scipy.ndimage.morphology import distance_transform_edt

    orig_nib = nib.load(orig_file)
    dest_nib = nib.load(dest_file)

    orig = orig_nib.get_data()
    dest = dest_nib.get_data()

    dist_orig = distance_transform_edt(np.logical_not(orig.astype(np.bool)))
    dist_dest = distance_transform_edt(np.logical_not(dest.astype(np.bool)))

    # normalized distance (0 in origin to 1 in dest)
    ndist = dist_orig / (dist_orig + dist_dest)

    ndist_nib = nib.Nifti1Image(ndist.astype(np.float32), orig_nib.affine)
    nib.save(ndist_nib, 'ndist.nii.gz')

    return os.path.abspath('ndist.nii.gz')

def create_shells(ndist_file, n_shells=4, out_file = 'shells.nii.gz', mask_file=None):
    """creates specified number of shells given normalized distance map. When mask is given, output in mask == 0 is set to zero"""
    import os
    import nibabel as nib
    import numpy as np

    ndist_nib = nib.load(ndist_file)
    ndist = ndist_nib.get_data()

    # if mask is provided, use it to mask-out regions outside it
    if mask_file is not None:
        mask_nib = nib.load(mask_file)
        assert mask_nib.header.get_data_shape() == ndist_nib.header.get_data_shape(), "Different shapes of images"
        mask = mask_nib.get_data() > 0

    out = np.zeros(ndist.shape, dtype=np.int8)

    limits = np.linspace(0., 1., n_shells+1)
    for i in np.arange(n_shells)+1:
        # compute shell and assing increasing label-id
        mask2 = np.logical_and(ndist >= limits[i-1], ndist < limits[i])
        if mask_file is not None:  # maskout regions outside mask
            mask2 = np.logical_and(mask2, mask)
        out[mask2] = i
    out[np.isclose(ndist, 0.)] = 0  # need to assign zero to ventricles because of >= above

    aux_hdr = ndist_nib.header
    aux_hdr.set_data_dtype(np.int8)

    out_nib = nib.Nifti1Image(out, ndist_nib.affine, aux_hdr)
    nib.save(out_nib, out_file)

    return os.path.abspath(out_file)


def merge_labels(in1_file, in2_file, out_file='merged.nii.gz', intersect=False, num_digits=1):
    """merges labels from two input labelmaps, optionally computing intersection"""
    import os
    import nibabel as nib
    import numpy as np

    in1_nib = nib.load(in1_file)
    in2_nib = nib.load(in2_file)

    assert in1_nib.header.get_data_shape() == in2_nib.header.get_data_shape(), "Different shapes of images"

    in1 = in1_nib.get_data()
    in2 = in2_nib.get_data()

    out = None

    # if not intersection, simply include labels from 'in2' into 'in1'
    if not intersect:

        out = np.zeros(in1.shape, dtype=np.int8)

        out[:] = in1[:]
        mask = in2 > 0
        out[mask] = in2[mask]  # overwrite in1 where in2 > 0


        aux_hdr = in1_nib.header
        aux_hdr.set_data_dtype(np.int8)

    # if intersection, create new label-set as cartesian product of the two sets
    else:

        out = np.zeros(in1.shape, dtype=np.int32)

        u1_set = np.unique(in1.ravel())
        u2_set = np.unique(in2.ravel())

        for u1 in u1_set:
            if u1 == 0: continue
            str1 = '%d'
            mask1 = in1 == u1
            for u2 in u2_set:
                if u2 == 0: continue
                mask2 = in2 == u2
                mask3 = np.logical_and(mask1, mask2)
                if not np.any(mask3): continue
                str2 = '%%0%dd' % num_digits  # number of digits used to represent label id
                out[mask3] = int(str1 % u1 + str2 % u2)  # new label id by concatenating [u1, u2]

        aux_hdr = in1_nib.header
        aux_hdr.set_data_dtype(np.int32)

    out_nib = nib.Nifti1Image(out, in1_nib.affine, aux_hdr)
    nib.save(out_nib, out_file)

    return os.path.abspath(out_file)


def generate_wmparc(incl_file, ndist_file, label_file, incl_labels=None, verbose=False):
    """generates wmparc by propagating labels in 'label_file' down the gradient defined by distance map in 'ndist_file'.
    Labels are only propagated in regions where 'incl_file' > 0 (or 'incl_file' == incl_labels[i], if 'incl_labels is provided).
    """
    import os
    import nibabel as nib
    import numpy as np
    from scipy.ndimage.morphology import binary_dilation, generate_binary_structure, iterate_structure

    connectivity = generate_binary_structure(3, 2)

    # read images
    incl_nib = nib.load(incl_file)
    ndist_nib = nib.load(ndist_file)
    label_nib = nib.load(label_file)

    assert incl_nib.header.get_data_shape() == ndist_nib.header.get_data_shape() and \
           incl_nib.header.get_data_shape() == label_nib.header.get_data_shape(), "Different shapes of mask, ndist and label images"

    # create inclusion mask
    incl_mask = None
    incl_aux = incl_nib.get_data()
    if incl_labels is None:
        incl_mask = incl_aux > 0
    else:
        incl_mask = np.zeros(incl_nib.header.get_data_shape(), dtype=np.bool)
        for lab in incl_labels:
            incl_mask[incl_aux == lab] = True

    # get rest of numpy arrays
    ndist = ndist_nib.get_data()
    label = label_nib.get_data()

    # get DONE and processing masks
    DONE_mask = label > 0  # this is for using freesurfer wmparc
    proc_mask = np.logical_and(np.logical_and(ndist > 0., ndist < 1.), incl_mask)

    # setup the ouptut vol
    out = np.zeros(label.shape, dtype=label.dtype)

    # initialize labels in cortex
    out[DONE_mask] = label[DONE_mask]  # this is for using freesurfer wmparc

    # start with connectivity 1
    its_conn = 1

    # main loop
    while not np.all(DONE_mask[proc_mask]):

        if verbose:
            print('%0.1f done' % (100. * float(DONE_mask[proc_mask].sum()) / float(proc_mask.sum())))

        # loop to increase connectivity for non-reachable TO-DO points
        while True:

            # dilate the SOLVED area
            aux = binary_dilation(DONE_mask, iterate_structure(connectivity, its_conn))
            # next TO-DO: close to DONE, in the processing mask and not yet done
            TODO_mask = np.logical_and(np.logical_and(aux, proc_mask), np.logical_not(DONE_mask))

            if TODO_mask.sum() > 0:
                break

            if verbose:
                print('Non-reachable points. Increasing connectivity')

            its_conn += 1

        # sort TO-DO points by ndist
        Idx_TODO = np.argwhere(TODO_mask)
        Idx_ravel = np.ravel_multi_index(Idx_TODO.T, label.shape)
        I_sort = np.argsort(ndist.ravel()[Idx_ravel])

        # iterate along TO-DO points
        for idx in Idx_TODO[I_sort[::-1]]:

            max_dist = -1.

            # process each neighbor
            for off in np.argwhere(iterate_structure(connectivity, its_conn)) - its_conn:

                try:

                    # if it is not DONE then skip
                    if not DONE_mask[idx[0] + off[0], idx[1] + off[1], idx[2] + off[2]]:
                        continue

                    # if it is the largest distance (ie, largest gradient)
                    cur_dist = ndist[idx[0] + off[0], idx[1] + off[1], idx[2] + off[2]]
                    if cur_dist > max_dist:
                        out[idx[0], idx[1], idx[2]] = out[idx[0] + off[0], idx[1] + off[1], idx[2] + off[2]]
                        max_dist = cur_dist

                except:
                    print('something wrong with neighbor at: (%d, %d, %d)' % (
                    idx[0] + off[0], idx[1] + off[1], idx[2] + off[2]))
                    pass

            if max_dist < 0.: print("something went wrong with point: (%d, %d, %d)" % (idx[0], idx[1], idx[2]))

            # mark as solved and remove from visited
            DONE_mask[idx[0], idx[1], idx[2]] = True

    # # remove labels from cortex (old aparc version)
    # out[dest_mask] = 0

    print('Writing output labelmap')
    out_nib = nib.Nifti1Image(out, label_nib.affine, label_nib.header)
    nib.save(out_nib, 'wmparc.nii.gz')

    return os.path.abspath('wmparc.nii.gz')


class Annot2LabelInputSpec(CommandLineInputSpec):
    subject = traits.String(desc='subject id', argstr='--subject %s', position=0, mandatory=True)
    hemi = traits.Enum("rh", "lh", desc="hemisphere [rh | lh]", position=1, argstr="--hemi %s", mandatory=True)
    lobes = traits.Enum("lobes", desc='lobes type', argstr='--lobesStrict %s', position=2)
    in_annot = traits.File(desc='input annotation file', exists=True)

class Annot2LabelOutputSpec(TraitedSpec):
    out_annot_file = File(desc = "lobes annotation file", exists = True)

class Annot2Label(CommandLine):
    """wrapper for FreeSurfer command-line tool 'mri_annotation2label'"""
    input_spec = Annot2LabelInputSpec
    output_spec = Annot2LabelOutputSpec
    _cmd = os.path.join(os.environ['FREESURFER_HOME'], 'bin', 'mri_annotation2label')

    def _list_outputs(self):
            outputs = self.output_spec().get()
            outputs['out_annot_file'] = os.path.join(os.path.dirname(self.inputs.in_annot), self.inputs.hemi + ".lobes.annot")
            return outputs

    def _format_arg(self, name, spec, value):
        if(name=='subject'):
             # take only the last part of the subject path
             return spec.argstr % ( os.path.basename(os.path.normpath(self.inputs.subject)))

        return super(Annot2Label, self)._format_arg(name, spec, value)


class Aparc2AsegInputSpec(CommandLineInputSpec):
    subject = traits.String(desc='subject id', argstr='--s %s', position=0, mandatory=True)
    annot = traits.String(desc='name of annot file', argstr='--annot %s', position=1, mandatory=True)
    labelwm = traits.Bool(desc='percolate white matter', argstr='--labelwm', position=2)
    dmax = traits.Int(desc='depth to percolate', argstr='--wmparc-dmax %d', position=3)
    rip = traits.Bool(desc='rip unknown label', argstr='--rip-unknown', position=4)
    hypo = traits.Bool(desc='hypointensities as wm', argstr='--hypo-as-wm', position=5)
    out_file = traits.File(desc='output aseg file', argstr='--o %s', position=6)
    in_lobes_rh = traits.File(desc='input lobar file RH', exists=True)
    in_lobes_lh = traits.File(desc='input lobar file LH', exists=True)

class Aparc2AsegOutputSpec(TraitedSpec):
    out_file = File(desc = "lobes aseg file", exists = True)

class Aparc2Aseg(CommandLine):
    """wrapper for FreeSurfer command-line tool 'mri_aparc2aseg'"""
    input_spec = Aparc2AsegInputSpec
    output_spec = Aparc2AsegOutputSpec

    _cmd = os.path.join(os.environ['FREESURFER_HOME'], 'bin', 'mri_aparc2aseg')

    def _list_outputs(self):
            outputs = self.output_spec().get()
            outputs['out_file'] = os.path.abspath(self.inputs.out_file)
            return outputs

    def _format_arg(self, name, spec, value):
        if(name=='subject'):
             # take only the last part of the subject path
             return spec.argstr % ( os.path.basename(os.path.normpath(self.inputs.subject)))

        return super(Aparc2Aseg, self)._format_arg(name, spec, value)



