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
                tokens.append(('header' + line[0], last_line))
                skip = 2
                continue

            # Key value list
            match_star = re.search(rx_star_item, line)
            if match_star is not None:
                tokens.append(('key*', match_star.group(1)))
                tokens.append(('val*', match_star.group(2)))
                continue

            tokens.append(('text', line))

    # remove last three text tokens before header
    tokens2 = []
    for i in range(len(tokens)):
        tok_name, tok_value = tokens[i]
        if tok_name.startswith('header') and len(tokens2) > 0:
            tokens2.pop()
            for _ in range(2):
                if len(tokens2) == 0:
                    break
                if tokens[-1][0].startswith('text'):
                    tokens2.pop()
        tokens2.append(tokens[i])
    tokens = tokens2

    # merge text tokens into preceding value tokens
    tokens2 = []
    for i in range(len(tokens)):
        tok_name, tok_value = tokens[i]
        if i > 0 and tok_name.startswith('text') and tokens[i-1][0].startswith('val'):
            tokens[i - 1] = (tokens[i-1][0], tokens[i-1][1] + '\n' + tok_value)
        tokens2.append(tokens[i])
    tokens = tokens2

    # Parser

    document = {}
    section = {}
    key = ''

    for tok_name, tok_value in tokens:
        if tok_name.startswith('header'):
            section = {}
            document[tok_value] = section
            continue
        if tok_name.startswith('key'):
            key = tok_value
            continue
        if tok_name.startswith('val'):
            section[key] = tok_value
            continue

    return document
