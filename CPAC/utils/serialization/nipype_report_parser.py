from abc import ABC
from os import PathLike
from typing import Union, Dict

import docutils.frontend
import docutils.nodes
import docutils.parsers.rst
import docutils.utils


class ReportSection(ABC):
    EXECUTION_INPUTS = 'execution-inputs'
    EXECUTION_OUTPUTS = 'execution-outputs'
    EXECUTION_INFO = 'execution-info'


def _parse_rst(text: str) -> docutils.nodes.document:
    """
    Parse ReStructured Text
    from:
    https://stackoverflow.com/questions/12883428/how-to-parse-restructuredtext-in-python
    """
    parser = docutils.parsers.rst.Parser()
    components = (docutils.parsers.rst.Parser,)
    settings = docutils.frontend.OptionParser(components=components).get_default_values()
    document = docutils.utils.new_document('<rst-doc>', settings=settings)
    parser.parse(text, document)
    return document


def read_report_rst(filename: Union[str, PathLike]) -> Dict[str, Dict[str, str]]:
    """
    Read NiPype report.rst data.

    Parameters
    ----------
    filename : Filepath to report.rst

    Returns
    -------
    Nested dictionary of sections and key-value pairs.
    """

    with open(filename, encoding='utf8') as file:
        rst = file.read()

    try:
        doc = _parse_rst(rst).asdom()
    except:  # noqa
        print(f'Could not parse RST file: "{filename}"')
        # ToDo: Report failure to parse with filepath?
        return {}

    extract_sections = (
        ReportSection.EXECUTION_INPUTS,
        ReportSection.EXECUTION_OUTPUTS,
        ReportSection.EXECUTION_INFO
    )

    out_dict = {}

    for section in doc.getElementsByTagName('section'):
        section_title = section.getAttribute('ids')
        if section_title in extract_sections:
            if len(section.childNodes) < 2:
                continue
            out_dict[section_title] = {}
            for bullet in section.childNodes[1].childNodes:
                if len(bullet.childNodes) == 0 or len(bullet.firstChild.childNodes) == 0:
                    continue
                item_text = bullet.firstChild.firstChild.nodeValue
                if not isinstance(item_text, str) or ' : ' not in item_text:
                    continue
                key, val = item_text.split(' : ', maxsplit=1)
                out_dict[section_title][key] = val

    return out_dict
