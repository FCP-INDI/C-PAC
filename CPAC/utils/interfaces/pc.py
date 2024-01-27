from nipype import logging
from nipype.interfaces.afni.base import AFNICommand, AFNICommandInputSpec
from nipype.interfaces.base import File, TraitedSpec, traits

iflogger = logging.getLogger("nipype.interface")


class PCInputSpec(AFNICommandInputSpec):
    in_file = File(
        desc="input file to 3dpc",
        argstr="%s",
        position=-1,
        mandatory=True,
        exists=True,
        copyfile=False,
    )
    out_file = File(
        name_template="%s_pcs",
        desc="output image file name",
        argstr="-prefix %s",
        name_source="in_file",
    )
    mask = File(desc="input mask", argstr="-mask %s", exists=True)
    pcs_file = File(
        name_template="%s_vec.1D",
        desc="output image file name",
        name_source="out_file",
        keep_extension=True,
    )
    pcs = traits.Int(
        desc="number of components to save in the output;"
        "it cant be more than the number of input bricks",
        argstr="-pcsave %s",
    )


class PCOutputSpec(TraitedSpec):
    out_file = File(desc="mask file", exists=True)
    pcs_file = File(desc="pcs", exists=True)


class PC(AFNICommand):
    _cmd = "3dpc"
    input_spec = PCInputSpec
    output_spec = PCOutputSpec
