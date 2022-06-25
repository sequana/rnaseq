import easydev
import os
import tempfile
import subprocess
import sys

from . import test_dir


sharedir = f"{test_dir}/data"
saccer3 = f"{test_dir}/data/Saccer3/"


def test_standalone_subprocess():
    directory = tempfile.TemporaryDirectory()
    cmd = """sequana_rnaseq --input-directory {} --working-directory {} """.format(
        sharedir, directory.name)
    subprocess.call(cmd.split())


def test_standalone_script():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.rnaseq.main as m
    sys.argv = ["test", "--input-directory", sharedir, "--genome-directory",
        saccer3, "--force", "--aligner", "bowtie2", 
             "--feature-counts-feature-type", 'gene,tRNA',
        "--rRNA-feature", "rRNA_gene"]   # ideally should be rRNA but current
    m.main()

def test_standalone_script_contaminant():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.rnaseq.main as m
    sys.argv = ["test", "--input-directory", sharedir, "--genome-directory",
        saccer3, "--force", "--aligner", "bowtie2", 
             "--feature-counts-feature-type", 'gene,tRNA', 
            "--contaminant-file", "test.fa",
        "--rRNA-feature", "rRNA_gene"]   # ideally should be rRNA but current
    m.main()

def test_version():
    cmd = "sequana_rnaseq --version"
    subprocess.call(cmd.split())


def test_standalone_script_wrong_feature():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.rnaseq.main as m
    sys.argv = ["test", "--input-directory", sharedir, "--genome-directory",
        saccer3, "--force", "--aligner", "bowtie2", 
             "--feature-counts-feature-type", 'dummy', 
        "--rRNA-feature", "rRNA_gene"]   # ideally should be rRNA but current
    try:
        m.main()
        assert False
    except:
        assert True

def test_standalone_script_wrong_reference():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.rnaseq.main as m
    sys.argv = ["test", "--input-directory", sharedir, "--genome-directory",
        "dummy", "--force", "--aligner", "bowtie2", 
        "--rRNA-feature", "rRNA_gene"]   # ideally should be rRNA but current
    try:
        m.main()
        assert False
    except:
        assert True

def test_standalone_script_wrong_triming():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.rnaseq.main as m
    sys.argv = ["test", "--input-directory", sharedir, "--genome-directory",
        saccer3, "--force", "--aligner", "bowtie2", "--software-choice", "dummy",
        "--rRNA-feature", "rRNA_gene"]   # ideally should be rRNA but current
    try:
        m.main()
        assert False
    except SystemExit:
        assert True



def test_full():

    with tempfile.TemporaryDirectory() as directory:
        wk = directory

        cmd = f"sequana_rnaseq --input-directory {sharedir} --genome-directory {saccer3} --aligner bowtie2 --working-directory {wk} --force"
        subprocess.call(cmd.split())


        cmd = "snakemake -s rnaseq.rules --wrapper-prefix https://raw.githubusercontent.com/sequana/sequana-wrappers/  -p --cores 2 "

        stat = subprocess.call(cmd.split(), cwd=wk)

        assert os.path.exists(wk + "/summary.html")
        assert os.path.exists(wk + "/multiqc/multiqc_report.html")

