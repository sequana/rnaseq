import easydev
import os
import tempfile
import subprocess
import sys


sequana_path = easydev.get_package_location('sequana_rnaseq')
sharedir = os.sep.join([sequana_path , "sequana", 'data'])


def test_standalone_subprocess():
    directory = tempfile.TemporaryDirectory()
    #cmd = """sequana_pipelines_rnaseq --fastq-directory {} --output-directory {}""".format(
    #    sharedir, directory.name)
    #subprocess.call(cmd.split())


def test_standalone_script():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.rnaseq.main as m
    #sys.argv = ["test", "--fastq-directory", sharedir, "--output-directory", directory.name]
    #m.main()
