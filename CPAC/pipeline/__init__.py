"""The C-PAC pipeline and its underlying infrastructure

Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with C-PAC. If not, see <https://www.gnu.org/licenses/>."""
import os
import pkg_resources as p

ALL_PIPELINE_CONFIGS = os.listdir(
    p.resource_filename("CPAC", os.path.join("resources", "configs")))
ALL_PIPELINE_CONFIGS = [x.split('_')[2].replace('.yml', '') for
                        x in ALL_PIPELINE_CONFIGS if 'pipeline_config' in x]
ALL_PIPELINE_CONFIGS.sort()
AVAILABLE_PIPELINE_CONFIGS = [preconfig for preconfig in ALL_PIPELINE_CONFIGS
                              if preconfig != 'benchmark-ANTS' and
                              not preconfig.startswith('regtest-')]


__all__ = ['ALL_PIPELINE_CONFIGS', 'AVAILABLE_PIPELINE_CONFIGS']
