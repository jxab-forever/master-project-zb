"""Read Julia reference matrices from the full and compact distributions."""
import gzip
import hashlib
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import highspy


def load_reference_mps(solver, path):
    """Load a matrix and return the SHA-256 of its original, uncompressed bytes."""
    path = Path(path).resolve()
    compressed = path.with_suffix(path.suffix + '.gz')
    previous = Path.cwd()
    try:
        if path.is_file():
            os.chdir(path.parent)
            status = solver.readModel(path.name)
            with path.open('rb') as stream:
                checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
        else:
            # Some Windows HiGHS builds cannot read gzip directly. Give the
            # native reader a temporary MPS with an ASCII relative filename.
            with TemporaryDirectory(prefix='stream_julia_', dir=path.parent) as directory:
                unpacked = Path(directory) / 'reference.mps'
                with gzip.open(compressed, 'rb') as source, unpacked.open('wb') as target:
                    shutil.copyfileobj(source, target)
                try:
                    os.chdir(directory)
                    status = solver.readModel(unpacked.name)
                finally:
                    os.chdir(previous)
                with unpacked.open('rb') as stream:
                    checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
    finally:
        os.chdir(previous)
    if status == highspy.HighsStatus.kError:
        raise RuntimeError(f'Could not read Julia matrix: {path}')
    return checksum
