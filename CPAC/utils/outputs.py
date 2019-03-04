import pkg_resources as p
import pandas as pd


class Outputs():

    # Settle some things about the resource pool reference and the output directory
    reference_csv = p.resource_filename('CPAC', 'resources/cpac_outputs.csv')

    try:
        reference = pd.read_csv(reference_csv)
    except Exception as e:
        err = "\n[!] Could not access or read the cpac_outputs.csv " \
              "resource file:\n{0}\n\nError details {1}\n".format(reference_csv, e)
        raise Exception(err)

    # all outputs
    any = list(reference.Resource)

    # outputs marked as optional in the matrix file, but we want them to be
    # written out no matter what for a specific reason
    override_optional = list(
        reference[reference['Override optional'] == 'yes']['Resource'])

    # extra outputs that we don't write to the output directory, unless the
    # user selects to do so
    debugging = list(
        reference[reference['Optional outputs: Debugging outputs'] == 'yes']['Resource'])

    # outputs to write out if the user selects to write all the functional
    # resources and files CPAC generates
    extra_functional = list(
        reference[reference['Optional outputs: Extra functionals'] == 'yes']['Resource'])

    # outputs to send into smoothing, if smoothing is enabled, and
    # outputs to write out if the user selects to write non-smoothed outputs
    # "_mult" is for items requiring mapnodes
    _template_filter = reference['Space'] == 'template'
    _native_filter = reference['Optional outputs: Native space'] == 'yes'

    _nonsmoothed_filter = reference['Optional outputs: Non-smoothed'] == 'yes'
    _multiple_filter = reference['Multiple outputs'] == 'yes'
    _derivative_filter = reference['Derivative'] == 'yes'

    native_nonsmooth = list(reference[_native_filter & _nonsmoothed_filter & ~_multiple_filter]['Resource'])
    native_nonsmooth_mult = list(reference[_native_filter & _nonsmoothed_filter & _multiple_filter]['Resource'])
    template_nonsmooth = list(reference[_template_filter & _nonsmoothed_filter & ~_multiple_filter]['Resource'])
    template_nonsmooth_mult = list(reference[_template_filter & _nonsmoothed_filter & _multiple_filter]['Resource'])

    # don't write these, unless the user selects to write native-space outputs
    native_smooth = list(reference[~_template_filter & ~_nonsmoothed_filter & _derivative_filter]['Resource'])

    # ever used??? contains template-space, smoothed, both raw and z-scored
    template_smooth = list(reference[_template_filter & ~_nonsmoothed_filter & _derivative_filter]['Resource'])

    # outputs to send into z-scoring, if z-scoring is enabled, and
    # outputs to write out if user selects to write non-z-scored outputs
    # "_mult" is for items requiring mapnodes
    template_raw = list(
        reference[_template_filter & ~_multiple_filter & 
        (reference['Optional outputs: Raw scores'] == 'yes')]['Resource'])
    template_raw_mult = list(reference[_template_filter & _multiple_filter &
        (reference['Optional outputs: Raw scores'] == 'yes')]['Resource'])

    # outputs to send into the average calculation nodes
    # "_mult" is for items requiring mapnodes
    average = list(reference[~_multiple_filter & (reference['Calculate averages'] == 'yes')]['Resource'])
    average_mult = list(reference[_multiple_filter & (reference['Calculate averages'] == 'yes')]['Resource'])

    # outputs to link for QC pages
    qc = list(reference[_derivative_filter & _template_filter]['Resource'])
