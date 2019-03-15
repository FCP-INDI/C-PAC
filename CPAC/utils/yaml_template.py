from __future__ import absolute_import

import re
import yaml

noalias_dumper = yaml.dumper.SafeDumper
noalias_dumper.ignore_aliases = lambda self, data: True


def create_yaml_from_template(d, template):

    """Save dictionary to a YAML file, keeping the structure
    (such as first level comments and ordering) from the template
    
    It may not be fully robust to YAML structures, but it works
    for C-PAC config files!
    """

    output = ""

    with open(template, 'r') as dtf:
        d_default = yaml.load(dtf)

    empty_lines = 0
    with open(template, 'r') as f:
        for line in f:

            # keep empty lines
            if re.match(r'^$', line.strip()):
                if empty_lines < 2:
                    output += "\n"
                empty_lines += 1
                continue

            # keep comments
            if re.match(r'^#', line):
                output += line
                empty_lines = 0
                continue

            # keep fields and add values
            # not robust for any YAML, focusing on C-PAC config files
            key_group = re.match(r'^([a-zA-Z_-]+)\s*:', line)
            if key_group:
                key = key_group.group(1)
                if key in d and d[key] is not None:
                    pass
                elif key in d_default and d_default[key] is not None:
                    d[key] = d_default[key]
                else:
                    output += key + ":\n"
                    continue
                    
                empty_lines = 0
                default_flow_style = False

                # Flow style for non-dict list or empty list
                if type(d[key]) is list and \
                    (
                        (len(d[key]) > 0 and type(d[key][0]) is not dict) or \
                        (len(d[key]) == 0)
                    ):

                    default_flow_style = True

                output += yaml.dump(
                    { key: d[key] } ,
                    default_flow_style=default_flow_style,
                    Dumper=noalias_dumper
                ).strip("{}\n\r") + "\n"
                

    return output