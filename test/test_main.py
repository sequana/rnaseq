import os
import subprocess
import sys
import tempfile

from click.testing import CliRunner

from sequana_pipelines.rnaseq.main import main

from . import test_dir

sharedir = f"{test_dir}/data"
saccer3 = f"{test_dir}/data/Saccer3/"
conta = f"{test_dir}/data/Saccer3/Saccer3_rRNA.fa"


# fast
def test_standalone_subprocess():
    directory = tempfile.TemporaryDirectory()
    cmd = """sequana_rnaseq --input-directory {} --working-directory {} """.format(sharedir, directory.name)
    subprocess.call(cmd.split())


# slow
def test_standalone_script():
    directory = tempfile.TemporaryDirectory()

    runner = CliRunner()
    results = runner.invoke(
        main,
        [
            "--input-directory",
            sharedir,
            "--genome-directory",
            saccer3,
            "--force",
            "--aligner-choice",
            "bowtie2",
            "--feature-counts-feature-type",
            "gene,tRNA",
            "--working-directory",
            directory.name,
            "--rRNA-feature",
            "rRNA_gene",
        ],
    )  # ideally should be rRNA but current
    assert results.exit_code == 0


def test_standalone_script_contaminant():
    directory = tempfile.TemporaryDirectory()
    runner = CliRunner()
    results = runner.invoke(
        main,
        [
            "--input-directory",
            sharedir,
            "--genome-directory",
            saccer3,
            "--force",
            "--aligner-choice",
            "bowtie2",
            "--feature-counts-feature-type",
            "gene",
            "--contaminant-file",
            conta,
            "--working-directory",
            directory.name,
        ],
    )
    assert results.exit_code == 0


# fast
def test_version():
    cmd = "sequana_rnaseq --version"
    subprocess.call(cmd.split())


# fast
def test_standalone_script_wrong_feature():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.rnaseq.main as m

    sys.argv = [
        "test",
        "--input-directory",
        sharedir,
        "--genome-directory",
        saccer3,
        "--force",
        "--aligner-choice",
        "bowtie2",
        "--feature-counts-feature-type",
        "dummy",
        "--working-directory",
        directory.name,
        "--rRNA-feature",
        "rRNA_gene",
    ]  # ideally should be rRNA but current
    try:
        m.main()
        assert False
    except:
        assert True


# fast
def test_standalone_script_wrong_reference():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.rnaseq.main as m

    sys.argv = [
        "test",
        "--input-directory",
        sharedir,
        "--genome-directory",
        "dummy",
        "--force",
        "--aligner-choice",
        "bowtie2",
        "--working-directory",
        directory.name,
        "--rRNA-feature",
        "rRNA_gene",
    ]  # ideally should be rRNA but current
    try:
        m.main()
        assert False
    except:
        assert True


# fast
def test_standalone_script_wrong_triming():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.rnaseq.main as m

    sys.argv = [
        "test",
        "--input-directory",
        sharedir,
        "--genome-directory",
        saccer3,
        "--force",
        "--aligner-choice",
        "bowtie2",
        "--software-choice",
        "dummy",
        "--working-directory",
        directory.name,
        "--rRNA-feature",
        "rRNA_gene",
    ]  # ideally should be rRNA but current
    try:
        m.main()
        assert False
    except SystemExit:
        assert True


# slow
def test_full():

    with tempfile.TemporaryDirectory() as directory:
        wk = directory

        cmd = f"sequana_rnaseq --input-directory {sharedir} --genome-directory {saccer3} --aligner-choice bowtie2 --working-directory {wk} --force --rRNA-feature rRNA_gene"
        subprocess.call(cmd.split())

        cmd = "bash rnaseq.sh"

        stat = subprocess.call(cmd.split(), cwd=wk)

        assert os.path.exists(wk + "/summary.html")
        assert os.path.exists(wk + "/multiqc/multiqc_report.html")


# slow
def test_full_star():

    with tempfile.TemporaryDirectory() as directory:
        wk = directory

        cmd = f"sequana_rnaseq --input-directory {sharedir} --genome-directory {saccer3} --aligner-choice star --working-directory {wk} --force --rRNA-feature rRNA_gene"
        subprocess.call(cmd.split())

        cmd = "snakemake -s rnaseq.rules --wrapper-prefix https://raw.githubusercontent.com/sequana/sequana-wrappers/  -p --cores 2 "

        stat = subprocess.call(cmd.split(), cwd=wk)

        assert os.path.exists(wk + "/summary.html")
        assert os.path.exists(wk + "/multiqc/multiqc_report.html")


# slow
def __test_full_salmon():

    with tempfile.TemporaryDirectory() as directory:
        wk = directory

        cmd = f"sequana_rnaseq --input-directory {sharedir} --genome-directory {saccer3} --aligner-choice salmon --working-directory {wk} --force"
        subprocess.call(cmd.split())

        cmd = "snakemake -s rnaseq.rules --wrapper-prefix https://raw.githubusercontent.com/sequana/sequana-wrappers/  -p --cores 2 "

        stat = subprocess.call(cmd.split(), cwd=wk)

        assert os.path.exists(wk + "/summary.html")
        assert os.path.exists(wk + "/multiqc/multiqc_report.html")


# fast
def test_genome_options_mutually_exclusive():
    """--genome-directory and --genome-accession cannot be used together"""
    with tempfile.TemporaryDirectory() as directory:
        runner = CliRunner()
        results = runner.invoke(
            main,
            [
                "--input-directory",
                sharedir,
                "--genome-directory",
                saccer3,
                "--genome-accession",
                "GCF_000146045.2",
                "--force",
                "--aligner-choice",
                "bowtie2",
                "--working-directory",
                directory,
            ],
        )
        assert results.exit_code == 1


# fast
def test_genome_options_required():
    """One of --genome-directory or --genome-accession is mandatory"""
    with tempfile.TemporaryDirectory() as directory:
        runner = CliRunner()
        results = runner.invoke(
            main,
            [
                "--input-directory",
                sharedir,
                "--force",
                "--aligner-choice",
                "bowtie2",
                "--working-directory",
                directory,
            ],
        )
        assert results.exit_code == 1


# fast
def test_genome_accession(monkeypatch):
    """The downloaded genome directory ends up in the configuration file"""
    # the download itself is tested in test_download.py; here we only check that
    # the downloaded directory is the one stored in the configuration file
    accession = "GCF_000000000.1"

    def fake_download(acc, outdir=".", force=False):
        assert acc == accession
        return os.path.abspath(saccer3)

    from sequana_pipelines.rnaseq import download

    monkeypatch.setattr(download, "download_genome", fake_download)

    with tempfile.TemporaryDirectory() as directory:
        runner = CliRunner()
        results = runner.invoke(
            main,
            [
                "--input-directory",
                sharedir,
                "--genome-accession",
                accession,
                "--force",
                "--aligner-choice",
                "bowtie2",
                "--working-directory",
                directory,
            ],
        )
        assert results.exit_code == 0

        with open(f"{directory}/.sequana/config.yaml") as fin:
            assert os.path.abspath(saccer3) in fin.read()
