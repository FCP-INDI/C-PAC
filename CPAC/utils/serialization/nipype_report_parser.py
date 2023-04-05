import re
from abc import ABC
from os import PathLike
from typing import List, Tuple, Union, Dict


class ReportSection(ABC):
    EXECUTION_INPUTS = 'Execution Inputs'
    EXECUTION_OUTPUTS = 'Execution Outputs'
    EXECUTION_INFO = 'Runtime info'
    ENVIRONMENT = 'Environment'
    ORIGINAL_INPUTS = 'Original Inputs'


rx_star_item = r"^\*\s(\S+)\s:\s(.*)$"


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
    tokens: List[Tuple[str, str]] = []
    with open(filename, encoding='utf8') as file:

        # Lexer

        line = ''
        skip = 0

        while True:
            last_line = line
            line = file.readline()
            if skip > 0:
                skip -= 1
                continue
            if not line:
                break
            line = line[:-1]

            # Headings
            # More efficient with this order of expressions
            # noinspection PyChainedComparisons
            if len(line) > 3 and line[0] in ('=', '-', '~') and line.count(line[0]) == len(line):
                tokens.pop()
                if len(tokens) > 2:
                    for _ in range(2):
                        tokens.pop()
                tokens.append(('header' + line[0], last_line))
                skip = 2
                continue

            # Key value list
            match_star = re.search(rx_star_item, line)
            if match_star is not None:
                tokens.append(('key*', match_star.group(1)))
                tokens.append(('val*', match_star.group(2)))
                continue

            # Append to previous item
            if len(tokens) > 0 and tokens[-1][0] == 'val*':
                tokens[-1] = (tokens[-1][0], tokens[-1][1] + '\n' + line)
                continue

            tokens.append(('text', line))

    # Parser

    document = {}
    section = {}
    key = ''
    # val = ''

    for tok_name, tok_value in tokens:
        if tok_name.startswith('header'):
            section = {}
            document[tok_value] = section
            continue
        if tok_name.startswith('key'):
            key = tok_value
            continue
        if tok_name.startswith('val'):
            val = tok_value
            section[key] = val
            continue

    return document
