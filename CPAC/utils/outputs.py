import pkg_resources as p
import pandas as pd


class Outputs():

    # Settle some things about the resource pool reference and the output directory
    reference_csv = p.resource_filename('CPAC', 'resources/cpac_outputs.tsv')

    try:
        reference = pd.read_csv(reference_csv, delimiter='\t', keep_default_na=False)
    except Exception as e:
        err = "\n[!] Could not access or read the cpac_outputs.tsv " \
              "resource file:\n{0}\n\nError details {1}\n".format(reference_csv, e)
        raise Exception(err)

    # all outputs
    any = list(reference.Resource)

    # extra outputs that we don't write to the output directory, unless the
    # user selects to do so
    debugging = list(
        reference[reference['Optional: Debugging'] == 'Yes']['Resource']
    )

    # functional data that are 4D time series, instead of derivatives
    functional_timeseries = list(
        reference[reference['4D Time Series'] == 'Yes']['Resource']
    )

    anat = list(reference[reference['Sub-Directory'] == 'anat']['Resource'])
    func = list(reference[reference['Sub-Directory'] == 'func']['Resource'])

    # outputs to send into smoothing, if smoothing is enabled, and
    # outputs to write out if the user selects to write non-smoothed outputs
    _template_filter = reference['Space'] == 'template'
    _epitemplate_filter = reference['Space'] == 'EPI template'
    _symtemplate_filter = reference['Space'] == 'symmetric template'
    _T1w_native_filter = reference['Space'] == 'T1w'
    _bold_native_filter = reference['Space'] == 'functional'
    _long_native_filter = reference['Space'] == 'longitudinal T1w'
    _nonsmoothed_filter = reference['To Smooth'] == 'Yes'
    _zstd_filter = reference['To z-std'] == 'Yes'
    _corr_filter = reference['Type'] == 'correlation'

    all_template_filter = (_template_filter | _epitemplate_filter |
                           _symtemplate_filter)
    all_native_filter = (_T1w_native_filter | _bold_native_filter |
                         _long_native_filter)

    native_nonsmooth = list(reference[all_native_filter &
                                      _nonsmoothed_filter]['Resource'])
    template_nonsmooth = list(reference[all_template_filter &
                                        _nonsmoothed_filter]['Resource'])

    to_smooth = list(reference[_nonsmoothed_filter]['Resource'])
    to_zstd = list(reference[_zstd_filter & ~_corr_filter]['Resource'])
    to_fisherz = list(reference[_zstd_filter & _corr_filter]['Resource'])

    # don't write these, unless the user selects to write native-space outputs
    native_smooth = list(reference[~all_template_filter & ~_nonsmoothed_filter]['Resource'])

    # ever used??? contains template-space, smoothed, both raw and z-scored
    template_smooth = list(reference[all_template_filter & ~_nonsmoothed_filter]['Resource'])

    _bold_filter = reference['Type'] == 'bold'
    _ts_filter = reference['4D Time Series'] == 'Yes'
    bold_ts = list(reference[_bold_filter & _bold_native_filter & _ts_filter]['Resource'])

    # outputs to send into z-scoring, if z-scoring is enabled, and
    # outputs to write out if user selects to write non-z-scored outputs
    native_raw = list(reference[all_native_filter &
                      (reference['To z-std'] == 'Yes')]['Resource'])

    template_raw = list(reference[all_template_filter &
                        (reference['To z-std'] == 'Yes')]['Resource'])

    def _is_cifti(_file_key):
         return _file_key.upper().startswith('CIFTI ')
    ciftis = reference[reference.File.map(_is_cifti)][['Resource', 'File']]
    ciftis = {cifti.Resource: cifti.File.split(' ')[-1] for cifti in ciftis.itertuples() if ' ' in cifti.File}

    def _is_gifti(_file_key):
         return _file_key.upper().startswith('GIFTI ')
    giftis = reference[reference.File.map(_is_gifti)][['Resource', 'File']]
    giftis = {gifti.Resource: gifti.File.split(' ')[-1] for gifti in giftis.itertuples() if ' ' in gifti.File}
