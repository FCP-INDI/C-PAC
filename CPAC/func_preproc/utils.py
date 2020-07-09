
import nibabel as nb
import subprocess

def add_afni_prefix(tpattern):
    if tpattern:
        if ".txt" in tpattern:
            tpattern = "@{0}".format(tpattern)
    return tpattern


def nullify(value, function=None):
    from traits.trait_base import Undefined
    if value is None:
        return Undefined
    if function:
        return function(value)
    return value


def chunk_ts(func_file):
    func_img = nb.load(func)
    trs = func_img.shape[3]
    chunk = trs/n_cpus
    TR_ranges = []

    for chunk_idx in range(0, n_cpus):
        if chunk_idx == n_cpus - 1:
            TR_ranges.append((chunk_idx*chunk, trs - 1))
        else:
            TR_ranges.append((chunk_idx*chunk, (chunk_idx+1)*chunk -1))
    return TR_ranges


def split_ts_chunks(func_file, tr_ranges):
    if '.nii' in func_file:
        ext = '.nii'
    if '.nii.gz' in func_file:
        ext = '.nii.gz'

    split_funcs = []
    for chunk_idx, tr_range in enumerate(tr_ranges):
        out_file = os.path.join(os.getcwd(), os.path.basename(func_file).replace(ext, "_{0}{1}".format(chunk_idx, ext)))
        in_file = "{0}[{1}..{2}]".format(func_file, tr_range[0], tr_range[1])

        cmd = ["3dcalc", "-a", in_file, "-expr", "a", "-prefix", out_file]

        retcode = subprocess.check_output(cmd)

        split_funcs.append(out_file)

    return split_funcs
