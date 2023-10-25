import os
import subprocess
import pkg_resources as p

import numpy as np
import nibabel as nb
import numpy.ma as ma
import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import matplotlib.cm as cm
from matplotlib import gridspec as mgs
from matplotlib.colors import ListedColormap

from nipype.interfaces import afni
from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.utility as util


def generate_qc_pages(qc_dir):
    """Generates the QC HTML files populated with the QC images that were
    created during the CPAC pipeline run.

    This function runs after the pipeline is over.

    Parameters
    ----------
    qc_dir : string
        path to qc directory

    Returns
    -------
    None

    """

    qc_dir = os.path.abspath(qc_dir)

    try:
        if not os.path.exists(qc_dir):
            os.makedirs(qc_dir)
    except IOError:
        print("\n\n[!] Could not create a directory for the QC dashboard. "
              "Please check write permissions.\n\nDirectory attempted:\n "
              "{0}".format(qc_dir))
        raise IOError

    files = []
    for root, _, fs in os.walk(qc_dir):
        root = root[len(qc_dir) + 1:]
        files += [os.path.join(root, f) for f in fs]

    with open(p.resource_filename('CPAC.qc', 'data/index.html'), 'rb') as f:
        qc_content = f.read()
        qc_content = qc_content.replace(
            b'/*CPAC*/``/*CPAC*/',
            ('`' + '\n'.join(files) + '`').encode()
        )
        with open(os.path.join(qc_dir, 'index.html'), 'wb') as f:
            f.write(qc_content)


def cal_snr_val(measure_file):
    """Calculate average snr value for snr image.

    Parameters
    ----------
    measure_file : string
        path to input nifti file

    Returns
    -------
    avg_snr_file : string
        a text file store average snr value

    """

    data = nb.load(measure_file).get_fdata()
    data_flat = data.flatten()
    data_no0 = data_flat[data_flat > 0]
    snr_val = ma.mean(data_no0)

    avg_snr_file = os.path.join(os.getcwd(), 'average_snr_file.txt')
    with open(avg_snr_file, 'w') as f:
        f.write(str(snr_val) + '\n')

    return avg_snr_file


def gen_plot_png(arr, measure, ex_vol=None):
    """Generate Motion FD Plot. Shows which volumes were dropped.

    Parameters
    ----------
    arr : list
        Frame wise Displacements

    measure : string
        Label of the Measure

    ex_vol : list
        Volumes excluded

    Returns
    -------
    png_name : string
            path to the generated plot png
    """

    matplotlib.rcParams.update({'font.size': 8})

    arr = np.loadtxt(arr)

    if ex_vol:
        try:
            ex_vol = np.genfromtxt(ex_vol, delimiter=',', dtype=int)
            ex_vol = ex_vol[ex_vol > 0]
        except:
            ex_vol = []
    else:
        ex_vol = []

    arr = arr[1:]
    del_el = [x for x in ex_vol if x < len(arr)]

    ex_vol = np.array(del_el)

    fig = plt.figure(figsize=(10, 6))
    plt.plot([i for i in range(len(arr))], arr, '-')
    fig.suptitle('%s plot with Mean %s = %0.4f' % (measure, measure,
                                                   arr.mean()))
    if measure == 'FD' and len(ex_vol) > 0:

        plt.scatter(ex_vol, arr[ex_vol], c="red", zorder=2)

        for x in ex_vol:
            plt.annotate('( %d , %0.3f)' % (x, arr[x]), xy=(x, arr[x]),
                            arrowprops=dict(facecolor='black', shrink=0.0))

    plt.xlabel('Volumes')
    plt.ylabel('%s' % measure)
    png_name = os.path.join(os.getcwd(), '%s_plot.png' % measure)
    fig.savefig(os.path.join(os.getcwd(), png_name))
    plt.close()
    matplotlib.rcdefaults()
    return png_name


def gen_carpet_plt(gm_mask, wm_mask, csf_mask, functional_to_standard, output):

    size = (950, 800)

    carpet_plot_path = os.path.join(os.getcwd(), output + '.png')

    func = nb.load(functional_to_standard).get_fdata()
    gm_voxels = func[nb.load(gm_mask).get_fdata().astype(bool)]
    wm_voxels = func[nb.load(wm_mask).get_fdata().astype(bool)]
    csf_voxels = func[nb.load(csf_mask).get_fdata().astype(bool)]
    del func

    data = np.concatenate((gm_voxels, wm_voxels, csf_voxels))
    seg = np.concatenate((
        np.ones(gm_voxels.shape[0]) * 1,
        np.ones(wm_voxels.shape[0]) * 2,
        np.ones(csf_voxels.shape[0]) * 3
    ))

    p_dec = 1 + data.shape[0] // size[0]
    if p_dec:
        data = data[::p_dec, :]
        seg = seg[::p_dec]

    t_dec = 1 + data.shape[1] // size[1]
    if t_dec:
        data = data[:, ::t_dec]

    interval = max((int(data.shape[-1] + 1) // 10, int(data.shape[-1] + 1) // 5, 1))
    xticks = list(range(0, data.shape[-1])[::interval])

    mycolors = ListedColormap(cm.get_cmap('tab10').colors[:4][::-1])

    gs = mgs.GridSpecFromSubplotSpec(1, 2, subplot_spec=mgs.GridSpec(1, 1)[0],
                                    width_ratios=[1, 100],
                                    wspace=0.0)
    ax0 = plt.subplot(gs[0])
    ax0.set_yticks([])
    ax0.set_xticks([])
    ax0.imshow(seg[:, np.newaxis], interpolation='none', aspect='auto',
            cmap=mycolors, vmin=1, vmax=4)
    ax0.grid(False)
    ax0.spines["left"].set_visible(False)
    ax0.spines["top"].set_visible(False)

    ax1 = plt.subplot(gs[1])
    ax1.imshow(data, interpolation='nearest', aspect='auto', cmap='gray')
    ax1.grid(False)
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    ax1.set_xticks(xticks)
    ax1.set_xlabel('time (frames)')
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    plt.savefig(carpet_plot_path, dpi=200, bbox_inches='tight')
    plt.close()

    return carpet_plot_path


def gen_motion_plt(motion_parameters):

    """
    Function to Generate Matplotlib plot for motion.
    Separate plots for Translation and Rotation are generated.

    Parameters
    ----------

    motion_parameters : string
                    Motion Parameters file

    Returns
    -------

    translation_plot : string
        path to translation plot

    rotation_plot : string
        path to rotation plot

    """

    rotation_plot = os.path.join(os.getcwd(), 'motion_rot_plot.png')
    translation_plot = os.path.join(os.getcwd(), 'motion_trans_plot.png')

    data = np.loadtxt(motion_parameters).T

    plt.gca().set_prop_cycle(color=['red', 'green', 'blue'])
    plt.plot(data[0])
    plt.plot(data[1])
    plt.plot(data[2])
    plt.legend(['roll', 'pitch', 'yaw'], loc='upper right')
    plt.ylabel('Rotation (degrees)')
    plt.xlabel('Volume')
    plt.savefig(rotation_plot)
    plt.close()

    plt.gca().set_prop_cycle(color=['red', 'green', 'blue'])
    plt.plot(data[3])
    plt.plot(data[4])
    plt.plot(data[5])
    plt.legend(['x', 'y', 'z'], loc='upper right')
    plt.ylabel('Translation (mm)')
    plt.xlabel('Volume')
    plt.savefig(translation_plot)
    plt.close()

    return translation_plot, rotation_plot


def gen_histogram(measure_file, measure):
    """Generates Histogram Image of intensities for a given input nifti file.

    Parameters
    ----------
    measure_file : string
        path to input nifti file

    measure : string
        Name of the measure label in the plot

    Returns
    -------
    hist_path : string
        Path to the generated histogram png

    """
    hist_path = None

    from CPAC.qc.utils import make_histogram
    import os
    m_ = measure
    if isinstance(measure_file, list):
        hist_path = []
        for file_ in measure_file:
            measure = m_
            if 'sca_roi' in measure.lower():
                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                if 'ROI_' in fname:
                    fname = fname.rsplit('ROI_')[1]
                elif 'roi_' in fname:
                    fname = fname.rsplit('roi_')[1]
                fname = 'sca_roi_' + fname.split('_')[0]
                measure = fname
            if 'sca_tempreg' in measure.lower():
                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                fname = fname.split('z_maps_roi_')[1]
                fname = 'sca_mult_regression_maps_roi_' + fname.split('_')[0]
                measure = fname
            if 'dr_tempreg' in measure.lower():
                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                for i in ['temp_reg_map_', 'tempreg_map_', 'tempreg_maps_', 'temp_reg_maps_']:
                    if i in fname:
                        try:
                            fname = fname.rsplit(i)[1]
                            break
                        except IndexError:
                            continue
                fname = 'dual_regression_map_'+ fname.split('_')[0]
                measure = fname
            if 'centrality' in measure.lower():
                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                type_, fname = fname.split('centrality_')
                fname = type_ + 'centrality_' + fname.split('_')[0]
                measure = fname

            hist_path.append(make_histogram(file_, measure))

    else:
        hist_path = make_histogram(measure_file, measure)

    return hist_path


def make_histogram(measure_file, measure):

    """
    Generates Histogram Image of intensities for a given input
    nifti file.

    Parameters
    ----------

    measure_file : string

                path to input nifti file

    measure : string

        Name of the measure label in the plot


    Returns
    -------

    hist_path : string

        Path to the generated histogram png

    """

    import matplotlib
    from matplotlib import pyplot as plt
    import numpy as np
    import nibabel as nb
    import os

    data = nb.load(measure_file).get_fdata()
    data_flat = data.flatten(order='F')
    y, binEdges = np.histogram(
        data_flat[
            (np.isfinite(data_flat)) & (data_flat != 0)
        ],
        bins=100
    )
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

    fig = plt.figure()
    fig.suptitle('%s intensity plot' % measure)
    plt.plot(bincenters, y, '-')
    plt.xlabel('intensity')
    plt.ylabel('# of voxels')

    png_name = os.path.join(os.getcwd(), '%s_hist_plot.png' % measure)
    fig.savefig(os.path.join(os.getcwd(), png_name))

    plt.close()
    hist_path = os.path.join(os.getcwd(), png_name)

    """
    ###
    hist_file = os.path.join(os.getcwd(), '%s_hist_path_file.txt' % measure)
    fl = open(hist_file, 'w')
    fl.write(str(measure_file) + '\n')
    fl.write(str(hist_path) + '\n')

    fl.close()
    """

    return hist_path


def drop_percent(measure_file, percent):
    """
    Zeros out voxels in measure files whose intensity doesnt fall in percent
    of voxel intensities

    Parameters
    ----------

    measure_file : string
                Input nifti file

    percent : percentage of the voxels to keep


    Returns
    -------

    modified_measure_file : string
                    measure_file with 1 - percent voxels zeroed out
    """

    import os
    import nibabel as nb
    import numpy as np

    img = nb.load(measure_file)
    data = img.get_fdata()

    max_val = np.percentile(data[data != 0.0], percent)
    data[data >= max_val] = 0.0

    save_img = nb.Nifti1Image(data, header=img.header, affine=img.affine)

    if '.nii.gz' in measure_file:
        ext = '.nii.gz'
    else:
        ext = '.nii'

    f_name = os.path.basename(os.path.splitext(os.path.splitext(measure_file)[0])[0])
    saved_name = '%s_%d_%s' % (f_name, percent, ext)
    save_img.to_filename(saved_name)

    modified_measure_file = os.path.join(os.getcwd(),
                                         saved_name)

    return modified_measure_file


def get_spacing(across, down, dimension):

    """
    Get Spacing in slices to be selected for montage
    display varying in given dimension

    Parameters
    ----------

    across : integer
        # images placed horizontally in montage

    down : integer
        # images stacked vertically in montage


    Returns
    -------

    space : integer
        # of images to skip before displaying next one

    """

    space = 10

    prod = (across*down*space)

    if prod > dimension:
        while(across*down*space) > dimension:
            space -= 1
    else:
        while(across*down*space) < dimension:
            space += 1

    return space


def determine_start_and_end(data, direction, percent):

    """
    Determine start slice and end slice in data file in
    given direction with at least threshold percent of voxels
    at start and end slices.

    Parameters
    ----------

    data : string
        input nifti file

    direction : string
        axial or sagittal

    percent : float
        percent(from total) of non zero voxels at starting and ending slice


    Returns
    -------

    start : integer
            Index of starting slice

    end : integer
            Index of the last slice

    """

    x, y, z = data.shape

    xx1 = 0
    xx2 = x - 1
    zz1 = 0
    zz2 = z - 1
    total_non_zero_voxels = len(np.nonzero(data.flatten())[0])
    thresh = percent * float(total_non_zero_voxels)
    start = None
    end = None

    if 'axial' in direction:

        while(zz2 > 0):

            d = len(np.nonzero(data[:, :, zz2].flatten())[0])
            if float(d) > thresh:
                break

            zz2 -= 1

        while(zz1 < zz2):

            d = len(np.nonzero(data[:, :, zz1].flatten())[0])
            if float(d) > thresh:
                break

            zz1 += 1

        start =  zz1
        end = zz2

    else:
        while(xx2 > 0):
            d = len(np.nonzero(data[xx2, :, :].flatten())[0])
            if float(d) > thresh:
                break

            xx2 -= 1

        while(xx1 < xx2):

            d = len(np.nonzero(data[xx1, :, :].flatten())[0])
            if float(d) > thresh:
                break

            xx1 += 1

        start = xx1
        end = xx2

    return start, end


def montage_axial(overlay, underlay, png_name, cbar_name):
    """Draws Montage using overlay on Anatomical brain in Axial Direction,
    calls make_montage_axial.

    Parameters
    ----------

    overlay : string
            Nifi file

    underlay : string
            Nifti for Anatomical Brain

    cbar_name : string
            name of the cbar

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """

    pngs = None
    if isinstance(overlay, list):
        pngs = []
        for ov in overlay:
            fname = os.path.basename(os.path.splitext(os.path.splitext(ov)[0])[0])
            pngs.append(make_montage_axial(ov, underlay,
                                           fname + '_' + png_name, cbar_name))
    else:
        pngs = make_montage_axial(overlay, underlay, png_name, cbar_name)

    png_name = pngs

    return png_name


def make_montage_axial(overlay, underlay, png_name, cbar_name):
    """
    Draws Montage using overlay on Anatomical brain in Axial Direction

    Parameters
    ----------

    overlay : string
            Nifi file

    underlay : string
            Nifti for Anatomical Brain

    cbar_name : string
            name of the cbar

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """
    import os
    import matplotlib
    matplotlib.rcParams.update({'font.size': 5})
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid1 import ImageGrid
    import matplotlib.pyplot as plt
    import nibabel as nb
    import numpy as np

    Y = nb.load(underlay).get_fdata()
    X = nb.load(overlay).get_fdata()
    X = X.astype(np.float32)
    Y = Y.astype(np.float32)

    if 'skull_vis' in png_name:
        X[X < 20.0] = 0.0
    if 'skull_vis' in png_name or \
            't1_edge_on_mean_func_in_t1' in png_name or \
                'MNI_edge_on_mean_func_mni' in png_name:
        max_ = np.nanmax(np.abs(X.flatten()))
        X[X != 0.0] = max_

    z1, z2 = determine_start_and_end(Y, 'axial', 0.0001)
    spacing = get_spacing(6, 3, z2 - z1)
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))

    if ('snr' in png_name) or ('reho' in png_name) or \
            ('vmhc' in png_name) or ('sca_' in png_name) or \
            ('alff' in png_name) or ('centrality' in png_name) or \
            ('dr_tempreg' in png_name):
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode="single", cbar_pad=0.2,
                         direction="row")
    else:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, direction="row")

    zz = z1
    for i in range(6*3):
        if zz >= z2:
            break
        try:
            im = grid[i].imshow(np.rot90(Y[:, :, zz]), cmap=cm.Greys_r)
        except IndexError as e:
            # TODO: send this to the logger instead
            print("\n[!] QC Interface: Had a problem with creating the "
                  "axial montage for {0}\n\nDetails:{1}. This error might occur because of a registration error encountered while using ANTs.\
                  Please refer to the png image located in your working directory for more insight."
                  "\n".format(png_name, e))
            pass
        zz += spacing

    x, y, z = X.shape
    X[X == 0.0] = np.nan
    max_ = np.nanmax(np.abs(X.flatten()))

    zz = z1
    im = None
    for i in range(6*3):
        if zz >= z2:
            break

        try:
            if cbar_name == 'red_to_blue':
                im = grid[i].imshow(np.rot90(X[:, :, zz]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=0, vmax=max_)
            elif cbar_name == 'green':
                im = grid[i].imshow(np.rot90(X[:, :, zz]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=0, vmax=max_)
            else:
                im = grid[i].imshow(np.rot90(X[:, :, zz]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=- max_, vmax=max_)
        except IndexError as e:
            # TODO: send this to the logger instead
            print("\n[!] QC Interface: Had a problem with creating the "
                  "axial montage for {0}\n\nDetails:{1}.This error might occur because of a registration error encountered while using ANTs.\
                   Please refer to the image located in your working directory for more insight"
                  "\n".format(png_name, e))
            pass

        grid[i].axes.get_xaxis().set_visible(False)
        grid[i].axes.get_yaxis().set_visible(False)
        zz += spacing

    cbar = grid.cbar_axes[0].colorbar(im)

    if 'snr' in png_name:
        cbar.ax.set_yticks(np.linspace(0, max_, 8))

    elif ('reho' in png_name) or ('vmhc' in png_name) or \
            ('sca_' in png_name) or ('alff' in png_name) or \
            ('centrality' in png_name) or ('dr_tempreg' in png_name):
        cbar.ax.set_yticks(np.linspace(-max_, max_, 8))

    plt.axis("off")
    png_name = os.path.join(os.getcwd(), png_name)
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.close()

    matplotlib.rcdefaults()

    return png_name


def montage_sagittal(overlay, underlay, png_name, cbar_name):
    """
    Draws Montage using overlay on Anatomical brain in Sagittal Direction
    calls make_montage_sagittal

    Parameters
    ----------

    overlay : string
            Nifi file

    underlay : string
            Nifti for Anatomical Brain

    cbar_name : string
            name of the cbar

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """

    pngs = None

    if isinstance(overlay, list):
        pngs = []
        for ov in overlay:
            fname = os.path.basename(os.path.splitext(os.path.splitext(ov)[0])[0])
            pngs.append(make_montage_sagittal(ov, underlay, fname + '_' + png_name, cbar_name))
    else:
        pngs = make_montage_sagittal(overlay, underlay, png_name, cbar_name)
    png_name = pngs

    return png_name


def make_montage_sagittal(overlay, underlay, png_name, cbar_name):
    """
    Draws Montage using overlay on Anatomical brain in Sagittal Direction

    Parameters
    ----------

    overlay : string
            Nifi file

    underlay : string
            Nifti for Anatomical Brain

    cbar_name : string
            name of the cbar

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """
    import matplotlib
    import os
    import numpy as np

    matplotlib.rcParams.update({'font.size': 5})

    try:
        from mpl_toolkits.axes_grid1 import ImageGrid
    except:
        from mpl_toolkits.axes_grid import ImageGrid
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import nibabel as nb

    from mpl_toolkits.axes_grid1 import ImageGrid

    from CPAC.qc.utils import determine_start_and_end, get_spacing

    Y = nb.load(underlay).get_fdata()
    X = nb.load(overlay).get_fdata()
    X = X.astype(np.float32)
    Y = Y.astype(np.float32)

    if 'skull_vis' in png_name:
        X[X < 20.0] = 0.0
    if 'skull_vis' in png_name or \
            't1_edge_on_mean_func_in_t1' in png_name or \
                'MNI_edge_on_mean_func_mni' in png_name:
        max_ = np.nanmax(np.abs(X.flatten()))
        X[X != 0.0] = max_

    x1, x2 = determine_start_and_end(Y, 'sagittal', 0.0001)
    spacing = get_spacing(6, 3, x2 - x1)
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))

    if ('snr' in png_name) or ('reho' in png_name) or \
            ('vmhc' in png_name) or ('sca_' in png_name) or \
            ('alff' in png_name) or ('centrality' in png_name) or \
            ('dr_tempreg' in png_name):
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode="single", cbar_pad=0.5,
                         direction="row")
    else:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode=None, direction="row")

    xx = x1
    for i in range(6*3):
        if xx >= x2:
            break

        try:
            im = grid[i].imshow(np.rot90(Y[xx, :, :]), cmap=cm.Greys_r)
        except IndexError as e:
            # TODO: send this to the logger instead
            print("\n[!] QC Interface: Had a problem with creating the "
                  "sagittal montage for {0}\n\nDetails:{1}.This error might occur because of a registration error encountered while using ANTs\
                   Please refer to the image located in your working directory for more insight"
                  "\n".format(png_name, e))
            pass

        grid[i].get_xaxis().set_visible(False)
        grid[i].get_yaxis().set_visible(False)
        xx += spacing

    x, y, z = X.shape
    X[X == 0.0] = np.nan
    max_ = np.nanmax(np.abs(X.flatten()))
    xx = x1
    for i in range(6*3):
        if xx >= x2:
            break
        im = None

        try:
            if cbar_name == 'red_to_blue':
                im = grid[i].imshow(np.rot90(X[xx, :, :]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=0, vmax=max_)
            elif cbar_name == 'green':
                im = grid[i].imshow(np.rot90(X[xx, :, :]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=0, vmax=max_)
            else:
                im = grid[i].imshow(np.rot90(X[xx, :, :]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=- max_, vmax=max_)
        except IndexError as e:
            # TODO: send this to the logger instead
            print("\n[!] QC Interface: Had a problem with creating the "
                  "sagittal montage for {0}\n\nDetails:{1}.This error might occur because of a registration error encountered while using ANTs.\
                   Please refer to the image located in your working directory for more insight"
                  "\n".format(png_name, e))
            pass

        xx += spacing

    try:
        cbar = grid.cbar_axes[0].colorbar(im)

        if 'snr' in png_name:
            cbar.ax.set_yticks(np.linspace(0, max_, 8))
        elif ('reho' in png_name) or ('vmhc' in png_name) or \
                ('sca_' in png_name) or ('alff' in png_name) or \
                ('centrality' in png_name) or ('dr_tempreg' in png_name):
            cbar.ax.set_yticks(np.linspace(-max_, max_, 8))

    except AttributeError as e:
        # TODO: send this to the logger instead
        print("\n[!] QC Interface: Had a problem with creating the "
              "sagittal montage for {0}\n\nDetails:{1}"
              "\n".format(png_name, e))
        pass

    plt.axis("off")
    png_name = os.path.join(os.getcwd(), png_name)
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.close()
    matplotlib.rcdefaults()

    return png_name


def montage_gm_wm_csf_axial(overlay_csf, overlay_wm, overlay_gm, underlay, png_name):

    """
    Draws Montage using GM WM and CSF overlays on Anatomical brain in Sagittal Direction

    Parameters
    ----------

    overlay_csf : string
            Nifi file CSF MAP

    overlay_wm : string
            Nifti file WM MAP

    overlay_gm : string
            Nifti file GM MAP

    underlay : string
            Nifti for Anatomical Brain

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid as ImageGrid
    import matplotlib.pyplot as plt
    import nibabel as nb
    import matplotlib.cm as cm

    Y = nb.load(underlay).get_fdata()
    z1, z2 = determine_start_and_end(Y, 'axial', 0.0001)
    spacing = get_spacing(6, 3, z2 - z1)
    X_csf = nb.load(overlay_csf).get_fdata()
    X_wm = nb.load(overlay_wm).get_fdata()
    X_gm = nb.load(overlay_gm).get_fdata()
    X_csf = X_csf.astype(np.float32)
    X_wm = X_wm.astype(np.float32)
    X_gm = X_gm.astype(np.float32)
    Y = Y.astype(np.float32)

    max_csf = np.nanmax(np.abs(X_csf.flatten()))
    X_csf[X_csf != 0.0] = max_csf
    max_wm = np.nanmax(np.abs(X_wm.flatten()))
    X_wm[X_wm != 0.0] = max_wm
    max_gm = np.nanmax(np.abs(X_gm.flatten()))
    X_gm[X_gm != 0.0] = max_gm
    fig = plt.figure(1)

    try:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                          aspect=True, cbar_mode=None, direction="row")
    except:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode=None, direction="row")

    zz = z1
    for i in range(6*3):
        if zz >= z2:
            break
        im = grid[i].imshow(np.rot90(Y[:, :, zz]), cmap=cm.Greys_r)
        zz += spacing

    x, y, z = X_csf.shape
    X_csf[X_csf == 0.0] = np.nan
    X_wm[X_wm == 0.0] = np.nan
    X_gm[X_gm == 0.0] = np.nan

    zz = z1
    im = None
    for i in range(6*3):
        if zz >= z2:
            break
        im = grid[i].imshow(np.rot90(X_csf[:, :, zz]), cmap=cm.get_cmap('green'), alpha=0.82, vmin=0, vmax=max_csf)
        im = grid[i].imshow(np.rot90(X_wm[:, :, zz]), cmap=cm.get_cmap('blue'), alpha=0.82, vmin=0, vmax=max_wm)
        im = grid[i].imshow(np.rot90(X_gm[:, :, zz]), cmap=cm.get_cmap('red'), alpha=0.82, vmin=0, vmax=max_gm)

        grid[i].axes.get_xaxis().set_visible(False)
        grid[i].axes.get_yaxis().set_visible(False)
        zz += spacing

    cbar = grid.cbar_axes[0].colorbar(im)

    plt.axis("off")
    png_name = os.path.join(os.getcwd(), png_name)
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.close()

    return png_name


def montage_gm_wm_csf_sagittal(overlay_csf, overlay_wm, overlay_gm, underlay, png_name):
    """
    Draws Montage using GM WM and CSF overlays on Anatomical brain in Sagittal Direction

    Parameters
    ----------

    overlay_csf : string
            Nifi file CSF MAP

    overlay_wm : string
            Nifti file WM MAP

    overlay_gm : string
            Nifti file GM MAP

    underlay : string
            Nifti for Anatomical Brain

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """

    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid as ImageGrid
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import nibabel as nb

    Y = nb.load(underlay).get_fdata()
    x1, x2 = determine_start_and_end(Y, 'sagittal', 0.0001)
    spacing = get_spacing(6, 3, x2 - x1)
    X_csf = nb.load(overlay_csf).get_fdata()
    X_wm = nb.load(overlay_wm).get_fdata()
    X_gm = nb.load(overlay_gm).get_fdata()
    X_csf = X_csf.astype(np.float32)
    X_wm = X_wm.astype(np.float32)
    X_gm = X_gm.astype(np.float32)
    Y = Y.astype(np.float32)

    max_csf = np.nanmax(np.abs(X_csf.flatten()))
    X_csf[X_csf != 0.0] = max_csf
    max_wm = np.nanmax(np.abs(X_wm.flatten()))
    X_wm[X_wm != 0.0] = max_wm
    max_gm = np.nanmax(np.abs(X_gm.flatten()))
    X_gm[X_gm != 0.0] = max_gm
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))

    try:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                          aspect=True, cbar_mode=None, direction="row")
    except:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode=None, direction="row")

    zz = x1
    for i in range(6*3):
        if zz >= x2:
            break
        im = grid[i].imshow(np.rot90(Y[zz, :, :]), cmap=cm.Greys_r)
        zz += spacing

    x, y, z = X_csf.shape
    X_csf[X_csf == 0.0] = np.nan
    X_wm[X_wm == 0.0] = np.nan
    X_gm[X_gm == 0.0] = np.nan

    zz = x1
    im = None
    for i in range(6*3):
        if zz >= x2:
            break

        im = grid[i].imshow(np.rot90(X_csf[zz, :, :]),
                            cmap=cm.get_cmap('green'), alpha=0.82, vmin=0,
                            vmax=max_csf)
        im = grid[i].imshow(np.rot90(X_wm[zz, :, :]),
                            cmap=cm.get_cmap('blue'), alpha=0.82, vmin=0,
                            vmax=max_wm)
        im = grid[i].imshow(np.rot90(X_gm[zz, :, :]),
                            cmap=cm.get_cmap('red'), alpha=0.82, vmin=0,
                            vmax=max_gm)

        grid[i].axes.get_xaxis().set_visible(False)
        grid[i].axes.get_yaxis().set_visible(False)
        zz += spacing

    cbar = grid.cbar_axes[0].colorbar(im)

    plt.axis("off")
    png_name = os.path.join(os.getcwd(), png_name)
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.close()

    return png_name


def register_pallete(colors_file, cbar_name):

    """
    Registers color pallete to matplotlib

    Parameters
    ----------

    colors_file : string
        file containing colors in hexadecimal formats in each line

    cbar_name : string
        Proposed name for the color bar


    Returns
    -------

    None

    """

    import matplotlib.colors as col
    import matplotlib.cm as cm

    with open(colors_file, 'r') as f:
        colors = [c.rstrip('\r\n') for c in reversed(f.readlines())]
        cmap3 = col.ListedColormap(colors, cbar_name)
        cm.register_cmap(cmap=cmap3)


def resample_1mm(file_):

    """
    Calls make_resample_1mm which resamples file to 1mm space

    Parameters
    ----------

    file_ : string
        path to the scan

    Returns
    -------

    new_fname : string
        path to 1mm resampled nifti file

    """
    new_fname = None

    if isinstance(file_, list):
        new_fname = []
        for f in file_:
            new_fname.append(make_resample_1mm(f))
    else:
        new_fname = make_resample_1mm(file_)

    return new_fname


def make_resample_1mm(file_):

    """
    Resamples input nifti file to 1mm space

    Parameters
    ----------

    file_ : string
        Input Nifti File

    Returns
    -------

    new_fname : string
            Input Nifti resampled to 1mm space
    """
    remainder, ext_ = os.path.splitext(file_)
    remainder, ext1_ = os.path.splitext(remainder)

    ext = ''.join([ext1_, ext_])

    new_fname = ''.join([remainder, '_1mm', ext])
    new_fname = os.path.join(os.getcwd(), os.path.basename(new_fname))
    cmd = " 3dresample -dxyz 1.0 1.0 1.0 -prefix %s " \
          "-inset %s " % (new_fname, file_)
    subprocess.getoutput(cmd)

    return new_fname

# own modules

# code
def dc(input1, input2):
    """
    Dice coefficient

    Computes the Dice coefficient (also known as Sorensen index) between the binary
    objects in two images.

    The metric is defined as

    .. math::

        DC=\frac{2|A\cap B|}{|A|+|B|}

    , where :math:`A` is the first and :math:`B` the second set of samples (here: binary objects).

    Parameters
    ----------
    input1 : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.

    input2 : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.

    Returns
    -------
    dc : float
        The Dice coefficient between the object(s) in ```input1``` and the
        object(s) in ```input2```. It ranges from 0 (no overlap) to 1 (perfect overlap).

    Notes
    -----
    This is a real metric.
    """
    input1 = numpy.atleast_1d(input1.astype(bool))
    input2 = numpy.atleast_1d(input2.astype(bool))

    intersection = numpy.count_nonzero(input1 & input2)

    size_i1 = numpy.count_nonzero(input1)
    size_i2 = numpy.count_nonzero(input2)

    try:
        dc = 2. * intersection / float(size_i1 + size_i2)
    except ZeroDivisionError:
        dc = 0.0

    return dc


def jc(input1, input2):
    """
    Jaccard coefficient

    Computes the Jaccard coefficient between the binary objects in two images.

    Parameters
    ----------
    input1 : array_like
            Input data containing objects. Can be any type but will be converted
            into binary: background where 0, object everywhere else.
    input2 : array_like
            Input data containing objects. Can be any type but will be converted
            into binary: background where 0, object everywhere else.

    Returns
    -------
    jc : float
        The Jaccard coefficient between the object(s) in `input1` and the
        object(s) in `input2`. It ranges from 0 (no overlap) to 1 (perfect overlap).

    Notes
    -----
    This is a real metric.
    """
    input1 = numpy.atleast_1d(input1.astype(bool))
    input2 = numpy.atleast_1d(input2.astype(bool))

    intersection = numpy.count_nonzero(input1 & input2)
    union = numpy.count_nonzero(input1 | input2)

    jc = float(intersection) / float(union)

    return jc

def crosscorr(input1,input2):
    """
    cross correlation
    computer compute cross correction bewteen input mask
    """

    input1 = numpy.atleast_1d(input1.astype(bool))
    input2 = numpy.atleast_1d(input2.astype(bool))

    from scipy.stats.stats import pearsonr
    cc=pearsonr(input1,input2)
    return cc

def coverage(input1,input2):
    """
    estimate the coverage between  two mask
    """
    input1 = numpy.atleast_1d(input1.astype(bool))
    input2 = numpy.atleast_1d(input2.astype(bool))

    intsec=numpy.count_nonzero(input1 & input2)
    if numpy.sum(input1)> numpy.sum(input2):
        smallv=numpy.sum(input2)
    else:
        smallv=numpy.sum(input1)
    cov=float(intsec)/float(smallv)
    return cov
