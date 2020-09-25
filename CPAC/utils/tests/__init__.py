from git import Repo
from os import path


def get_BIDS_examples_dir(parent_dir=None):
    """
    Function to return the path to bids-examples, cloning if necessary

    Parameters
    ----------
    parent_dir: str or None

    Returns
    -------
    example_dir: str
    """
    example_dir = path.abspath(path.join(
        parent_dir if parent_dir is not None else path.dirname(__file__),
        'bids-examples'
    ))
    if not path.exists(example_dir):
        Repo.clone_from(
            "https://github.com/bids-standard/bids-examples.git",
            example_dir
        )
    return(example_dir)
