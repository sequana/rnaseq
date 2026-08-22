"""Tests of the NCBI genome/annotation download"""

import os
import shutil
import zipfile

import pytest

from sequana_pipelines.rnaseq import download as download_module
from sequana_pipelines.rnaseq.download import download_genome

ACCESSION = "GCF_000000000.1"


def create_fake_datasets(directory, fasta=True, gff=True, returncode=0):
    """Create an executable named 'datasets' that mimics the NCBI tool"""
    bindir = directory / "bin"
    bindir.mkdir()

    data = f"ncbi_dataset/data/{ACCESSION}"
    archive = directory / "reference.zip"
    with zipfile.ZipFile(archive, "w") as zipdata:
        if fasta:
            zipdata.writestr(f"{data}/{ACCESSION}_ASM1_genomic.fna", ">chr1\nACGT\n")
        if gff:
            zipdata.writestr(f"{data}/genomic.gff", "##gff-version 3\n")

    executable = bindir / "datasets"
    executable.write_text(
        "#!/bin/bash\n"
        f"if [ {returncode} -ne 0 ]; then echo 'Error: invalid accession' >&2; exit {returncode}; fi\n"
        'while [ "$#" -gt 0 ]; do\n'
        '    if [ "$1" = "--filename" ]; then shift; cp ' + str(archive) + ' "$1"; fi\n'
        "    shift\n"
        "done\n"
    )
    executable.chmod(0o755)
    return bindir


def create_fake_api(monkeypatch, archive, calls):
    """Replace the API downloader by a local copy of a fake archive"""

    def fake_api(accession, target):
        calls.append(accession)
        shutil.copy(archive, target)

    monkeypatch.setattr(download_module, "_download_with_api", fake_api)


def test_download_genome(tmp_path, monkeypatch):
    """The genome and its annotation are downloaded and renamed"""
    bindir = create_fake_datasets(tmp_path)
    monkeypatch.setenv("PATH", str(bindir), prepend=os.pathsep)

    genome_directory = download_genome(ACCESSION, outdir=str(tmp_path))

    assert genome_directory == str(tmp_path / ACCESSION)
    assert os.path.exists(f"{genome_directory}/{ACCESSION}.fa")
    assert os.path.exists(f"{genome_directory}/{ACCESSION}.gff")


def test_download_genome_reuse(tmp_path, monkeypatch):
    """Existing files are re-used instead of being downloaded again"""
    # no 'datasets' executable available: existing files must be re-used
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))

    genome_directory = tmp_path / ACCESSION
    genome_directory.mkdir()
    (genome_directory / f"{ACCESSION}.fa").write_text(">chr1\nACGT\n")
    (genome_directory / f"{ACCESSION}.gff").write_text("##gff-version 3\n")

    assert download_genome(ACCESSION, outdir=str(tmp_path)) == str(genome_directory)


def test_download_genome_no_datasets(tmp_path, monkeypatch):
    """Without the 'datasets' executable, the REST API is used"""
    # no 'datasets' executable: the REST API must be used instead
    create_fake_datasets(tmp_path)
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    calls = []
    create_fake_api(monkeypatch, tmp_path / "reference.zip", calls)

    genome_directory = download_genome(ACCESSION, outdir=str(tmp_path))

    assert calls == [ACCESSION]
    assert os.path.exists(f"{genome_directory}/{ACCESSION}.fa")
    assert os.path.exists(f"{genome_directory}/{ACCESSION}.gff")


def test_download_genome_fallback_on_error(tmp_path, monkeypatch):
    """When 'datasets' fails, the REST API takes over"""
    # the 'datasets' executable fails: the REST API must take over
    bindir = create_fake_datasets(tmp_path, returncode=1)
    monkeypatch.setenv("PATH", str(bindir), prepend=os.pathsep)
    calls = []
    create_fake_api(monkeypatch, tmp_path / "reference.zip", calls)

    genome_directory = download_genome(ACCESSION, outdir=str(tmp_path))

    assert calls == [ACCESSION]
    assert os.path.exists(f"{genome_directory}/{ACCESSION}.fa")


def test_download_genome_error(tmp_path, monkeypatch):
    """When all methods fail, nothing is left behind"""
    # both methods fail
    bindir = create_fake_datasets(tmp_path, returncode=1)
    monkeypatch.setenv("PATH", str(bindir), prepend=os.pathsep)

    def fake_api(accession, target):
        raise RuntimeError("404 Client Error")

    monkeypatch.setattr(download_module, "_download_with_api", fake_api)

    with pytest.raises(RuntimeError, match="failed"):
        download_genome(ACCESSION, outdir=str(tmp_path))

    # nothing must be left behind
    assert not os.path.exists(tmp_path / ACCESSION)


def test_download_with_api(tmp_path, monkeypatch):
    """A valid archive downloaded from the API is kept as it is"""
    create_fake_datasets(tmp_path)
    content = (tmp_path / "reference.zip").read_bytes()

    class FakeResponse:
        """Minimal stand-in for a requests response"""

        def raise_for_status(self):
            """No HTTP error"""

        def iter_content(self, chunk_size=1):
            """Yield the whole payload at once"""
            yield content

    import requests

    monkeypatch.setattr(requests, "get", lambda *args, **kwargs: FakeResponse())

    archive = str(tmp_path / "downloaded.zip")
    download_module._download_with_api(ACCESSION, archive)
    assert zipfile.is_zipfile(archive)


def test_download_with_api_not_a_zip(tmp_path, monkeypatch):
    """A JSON error message returned by the API is reported"""

    class FakeResponse:
        """Minimal stand-in for a requests response"""

        def raise_for_status(self):
            """No HTTP error"""

        def iter_content(self, chunk_size=1):
            """Yield the whole payload at once"""
            yield b'{"message": "invalid accession"}'

    import requests

    monkeypatch.setattr(requests, "get", lambda *args, **kwargs: FakeResponse())

    with pytest.raises(RuntimeError, match="invalid accession"):
        download_module._download_with_api(ACCESSION, str(tmp_path / "downloaded.zip"))


def test_download_genome_no_annotation(tmp_path, monkeypatch):
    """An assembly without annotation is rejected"""
    bindir = create_fake_datasets(tmp_path, gff=False)
    monkeypatch.setenv("PATH", str(bindir), prepend=os.pathsep)

    with pytest.raises(RuntimeError, match="No annotation"):
        download_genome(ACCESSION, outdir=str(tmp_path))


def test_download_with_api_truncated_archive(tmp_path, monkeypatch):
    """A truncated archive returned by the API is reported"""

    # with an unknown accession, NCBI returns a truncated (binary) archive
    class FakeResponse:
        """Minimal stand-in for a requests response"""

        def raise_for_status(self):
            """No HTTP error"""

        def iter_content(self, chunk_size=1):
            """Yield the whole payload at once"""
            yield b"PK\x03\x04-\x00\x08\x00\x08\x00\xff\xff\xff\xff"

    import requests

    monkeypatch.setattr(requests, "get", lambda *args, **kwargs: FakeResponse())

    with pytest.raises(RuntimeError, match="accession may be incorrect"):
        download_module._download_with_api(ACCESSION, str(tmp_path / "downloaded.zip"))


def test_safe_extract_rejects_path_traversal(tmp_path):
    """Archive members pointing outside of the target directory are refused"""
    archive = tmp_path / "malicious.zip"
    with zipfile.ZipFile(archive, "w") as zipdata:
        zipdata.writestr("../escaped.fna", ">chr1\nACGT\n")

    target = tmp_path / "extraction"
    target.mkdir()
    with zipfile.ZipFile(archive) as zipdata:
        with pytest.raises(RuntimeError, match="Unexpected path"):
            download_module._safe_extract(zipdata, str(target))

    assert not os.path.exists(tmp_path / "escaped.fna")
