"""Helper functions to parse BED files."""
import pathlib

from .bed import Bed


def bed_objects_from_file(path):
    """Load a BED file and return a list of Bed objects""
    Args:
        path: Path to the BED file.
    """
    path = validate_file(path)
    return [Bed(line) for line in path.open()]


def validate_file(path):
    "Check if the path is a valid BED file and return a pathlib.Path object."
    path = pathlib.Path(path)

    # check if path is a valid file
    if not path.is_file():
        raise Exception(f'{path} is not a valid path')

    # check if path is a bed file
    # TODO this is just prelimenary check
    #if not path.suffix == '.bed':
    #    raise Exception(f'{path} is not a BED file')

    return path
