import scipy.signal
import nibabel as nb

def convolve(functional_file_path, function_file_paths):

    func_img = nb.load(functional_file_path)
    func_data = func_img.get_data()

    files = []
    for i, function_file_path in enumerate(function_file_paths):
        function = nb.load(function_file_path).get_data()
        convolved = scipy.signal.fftconvolve(func_data, function, axes=[-1], mode='same')
        function_img = nb.Nifti1Image(convolved, affine=func_img.affine)
        function_img.to_filename('./convolved_%d.nii.gz' % i)
        files += ['./convolved_%d.nii.gz' % i]
    
    return files