#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Download of a genome (fasta) and its annotation (GFF) from NCBI.

The download is performed with the NCBI *datasets* command line tool if it is
available, otherwise with the NCBI datasets REST API (pure Python), which does
not require any external software. Both methods retrieve the very same zip
archive.

The files are renamed to fulfill the pipeline convention, that is a directory
named after the accession, containing <accession>.fa and <accession>.gff
"""

import os
import shutil

# only fixed commands are executed, without a shell
import subprocess  # nosec B404
import tempfile
import zipfile

from sequana_pipetools import logger

__all__ = ["download_genome"]

NCBI_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"


def _download_with_datasets(accession, archive):
    """Download the data package using the NCBI *datasets* executable"""
    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        accession,
        "--include",
        "genome,gff3",
        "--filename",
        archive,
        "--no-progressbar",
    ]
    try:
        # the command is a fixed list of arguments; no shell is involved
        subprocess.run(cmd, check=True, capture_output=True, text=True, shell=False)  # nosec B603
    except subprocess.CalledProcessError as err:
        raise RuntimeError((err.stderr or err.stdout or "").strip() or "the 'datasets' command failed") from err


def _download_with_api(accession, archive):
    """Download the data package using the NCBI datasets REST API

    This is the fallback used when the *datasets* executable is not available.
    Only requests (a sequana dependency) is required here.
    """
    import requests

    url = f"{NCBI_API}/genome/accession/{accession}/download"
    params = {"include_annotation_type": ["GENOME_FASTA", "GENOME_GFF"]}

    try:
        response = requests.get(url, params=params, stream=True, timeout=(30, 300))
        response.raise_for_status()
        with open(archive, "wb") as fout:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                fout.write(chunk)
    except requests.RequestException as err:
        raise RuntimeError(str(err)) from err

    if not zipfile.is_zipfile(archive):
        # on errors, the API returns either a JSON message or a truncated archive
        # (e.g. with an unknown accession) rather than a valid zip archive
        with open(archive, "rb") as fin:
            content = fin.read(500)
        try:
            message = content.decode().strip()
        except UnicodeDecodeError:
            message = "invalid archive returned by NCBI. The accession may be incorrect"
        raise RuntimeError(message)


def _safe_extract(zipdata, target):
    """Extract a zip archive, refusing members pointing outside of *target*

    The archive comes from NCBI but an unvalidated extraction would allow a
    corrupted or malicious archive to write anywhere on the file system.
    """
    target = os.path.realpath(target)
    for member in zipdata.namelist():
        destination = os.path.realpath(os.path.join(target, member))
        if destination != target and not destination.startswith(target + os.sep):
            raise RuntimeError(f"Unexpected path ({member}) found in the NCBI archive. Extraction aborted.")
    zipdata.extractall(target)


def download_genome(accession, outdir=".", force=False):
    """Download a genome and its annotation from NCBI using *datasets*.

    :param str accession: a NCBI assembly accession (e.g. GCF_000146045.2)
    :param str outdir: where the <accession> directory is created
    :param bool force: download again even though the files already exist
    :return: the absolute path of the genome directory

    The resulting directory contains the <accession>.fa and <accession>.gff
    files expected by the pipeline.
    """
    genome_directory = os.path.abspath(os.path.join(outdir, accession))
    fasta = os.path.join(genome_directory, f"{accession}.fa")
    gff = os.path.join(genome_directory, f"{accession}.gff")

    if not force and os.path.exists(fasta) and os.path.exists(gff):
        logger.info(f"Found existing {fasta} and {gff}. Using them (use --force-genome-download to download again)")
        return genome_directory

    # the 'datasets' executable is used when available; the REST API is used as a
    # fallback so that no external software is required
    if shutil.which("datasets"):
        methods = [("datasets executable", _download_with_datasets), ("datasets API", _download_with_api)]
    else:
        methods = [("datasets API", _download_with_api)]

    with tempfile.TemporaryDirectory() as tmpdir:
        archive = os.path.join(tmpdir, "ncbi_dataset.zip")

        for i, (name, method) in enumerate(methods):
            logger.info(f"Downloading {accession} from NCBI using the {name}. This may take a while.")
            try:
                method(accession, archive)
                break
            except RuntimeError as err:
                if i + 1 == len(methods):
                    raise RuntimeError(f"The download of {accession} failed ({name}):\n{err}") from err
                logger.warning(f"Download with the {name} failed ({err}). Trying with the {methods[i + 1][0]}.")

        with zipfile.ZipFile(archive) as zipdata:
            _safe_extract(zipdata, tmpdir)

        data_directory = os.path.join(tmpdir, "ncbi_dataset", "data", accession)
        if not os.path.exists(data_directory):
            raise RuntimeError(f"No data found for the accession {accession}. Please check its validity on NCBI.")

        # sorted so that the selected files do not depend on the file system ordering
        filenames = sorted(os.listdir(data_directory))
        fasta_files = [x for x in filenames if x.endswith((".fna", ".fa", ".fasta"))]
        gff_files = [x for x in filenames if x.endswith(".gff")]

        if not fasta_files:
            raise RuntimeError(f"No genome sequence (fasta) found for the accession {accession}.")
        if not gff_files:
            raise RuntimeError(
                f"No annotation (GFF) found for the accession {accession}. The rnaseq pipeline requires an "
                "annotated genome. Please choose an annotated assembly or provide your own files with "
                "--genome-directory"
            )

        os.makedirs(genome_directory, exist_ok=True)
        shutil.move(os.path.join(data_directory, fasta_files[0]), fasta)
        shutil.move(os.path.join(data_directory, gff_files[0]), gff)

    logger.info(f"Genome and annotation saved in {genome_directory}")
    return genome_directory
