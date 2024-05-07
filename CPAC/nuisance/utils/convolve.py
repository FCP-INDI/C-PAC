import nibabel as nib
import scipy.signal


def convolve(functional_file_path, function_file_paths):
    func_img = nib.load(functional_file_path)
    func_data = func_img.get_fdata()

    files = []
    for i, function_file_path in enumerate(function_file_paths):
        function = nib.load(function_file_path).get_fdata()
        convolved = scipy.signal.fftconvolve(
            func_data, function, axes=[-1], mode="same"
        )
        function_img = nib.Nifti1Image(convolved, affine=func_img.affine)
        function_img.to_filename("./convolved_%d.nii.gz" % i)
        files += ["./convolved_%d.nii.gz" % i]

    return files
