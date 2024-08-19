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

        cmd = "snakemake -s rnaseq.rules --wrapper-prefix https://raw.githubusercontent.com/sequana/sequana-wrappers/  -p --cores 2 "

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
