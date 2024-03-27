# Copyright (C) 2018-2024  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
from logging import basicConfig, INFO

from CPAC.utils.monitoring.custom_logging import getLogger

logger = getLogger("CPAC.pipeline.test")
basicConfig(format="%(message)s", level=INFO)


def run_gather_outputs_func(pipeline_out_dir):
    from CPAC.pipeline import cpac_group_runner as cgr

    df_dct = cgr.gather_outputs(
        pipeline_out_dir, ["functional_to_standard"], None, False, False, get_func=True
    )
    logger.info(df_dct)
