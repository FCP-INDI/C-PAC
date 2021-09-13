"""FSL Nipype interfaces with customized functionality"""
import os
from nipype.interfaces.base import isdefined
from nipype.interfaces.fsl.utils import Merge as fslMerge


class Merge(fslMerge):
    """Use relative paths for input files"""
    def _format_arg(self, name, spec, value):
        if name == "tr":
            if self.inputs.dimension != "t":
                raise ValueError("When TR is specified, dimension must be t")
            return spec.argstr % value
        if name == "dimension":
            if isdefined(self.inputs.tr):
                return "-tr"
            return spec.argstr % value
        # Begin custom code
        # -----------------
        if name == "in_files":
            return ' '.join([
                os.path.relpath(in_file) for in_file in value
            ])
        # ---------------
        # End custom code
        return super(Merge, self)._format_arg(name, spec, value)
