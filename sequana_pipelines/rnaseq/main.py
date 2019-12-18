import sys
import os
import argparse

from sequana.pipelines_common import *
from sequana.snaketools import Module
from sequana import logger
logger.level = "INFO"

col = Colors()

NAME = "rnaseq"
m = Module(NAME)
m.is_executable()


class Options(argparse.ArgumentParser):
    def __init__(self, prog=NAME):
        usage = col.purple(
            """This script prepares the sequana pipeline rnaseq layout to
            include the Snakemake pipeline and its configuration file ready to
            use.

            In practice, it copies the config file and the pipeline into a
            directory (rnaseq) together with an executable script

            ::

                sequana_pipelines_rnaseq --input-directory PATH_TO_DATA

            You must provide a path to a genome of reference (with its GFF file).

            We save the genomes in the genome directory where lots of genomes
            can be stored. So, you also need to provide the name of a genome.

            For instance, if you work on  Saccer3, type::

                sequana_pipelines_rnaseq --input-directory \
                    --genome-directory GENOMESPATH/Saccer3

            The pipeline will search for the files Saccer3.fa and Saccer3.gff in
            the directory GENOMESPATH/Saccer3.


        """
        )
        super(Options, self).__init__(usage=usage, prog=prog, description="",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        # add a new group of options to the parser
        so = SlurmOptions()
        so.add_options(self)

        # add a snakemake group of options to the parser
        so = SnakemakeOptions(working_directory=NAME)
        so.add_options(self)

        so = InputOptions()
        so.add_options(self)

        so = GeneralOptions()
        so.add_options(self)

        pipeline_group = self.add_argument_group("pipeline_general")
        pipeline_group.add_argument("--genome-directory", dest="genome_directory",
            default=".", required=True)
        pipeline_group.add_argument("--aligner", dest="aligner", required=True,
            choices=['bowtie2', 'bowtie1', 'STAR', 'star'],
            help= "a mapper in bowtie, bowtie2, STAR")

        # cutadapt related
        so = CutadaptOptions()
        so.add_options(self)


        # fastq_screen
        pipeline_group = self.add_argument_group("pipeline_fastq_screen")
        pipeline_group.add_argument("--do-fastq-screen", action="store_true",
            help="run fastq_screen ")
        pipeline_group.add_argument("--fastq-screen-conf",
            default="fastq_screen.conf", type=str, 
            help="""a valid fastqc_screen.conf file. See fastq_screen
documentation for details. In a nutsheel, add a line for each genome you want to
search for in your input data. Each line is 'DATABASE name path BOWTIE2'. The
path includes the path to the genome + its prefix name.  """)


        pipeline_group = self.add_argument_group("pipeline_others")
        pipeline_group.add_argument('--do-igvtools', action="store_true")
        pipeline_group.add_argument('--do-bam-coverage', action="store_true")
        pipeline_group.add_argument('--do-mark-duplicates', action="store_true")
        pipeline_group.add_argument('--do-rnaseq-qc', action="store_true")




def main(args=None):

    if args is None:
        args = sys.argv

    if "--version" in sys.argv:
        print_version(NAME)
        sys.exit(0)

    options = Options(NAME).parse_args(args[1:])

    manager = PipelineManager(options, NAME)

    # create the beginning of the command and the working directory
    manager.setup()

    # fill the config file with input parameters
    cfg = manager.config.config

    # general
    cfg.general.genome_directory = os.path.abspath(options.genome_directory)
    cfg.general.aligner = options.aligner

    # cutadapt
    cfg.cutadapt.do = options.skip_cutadapt
    cfg.cutadapt.fwd = options.cutadapt_fwd
    cfg.cutadapt.rev = options.cutadapt_rev
    cfg.cutadapt.quality = options.cutadapt_quality
    cfg.cutadapt.tool_choice = options.cutadapt_tool_choice

    cfg.input_directory = os.path.abspath(options.input_directory)
    cfg.input_pattern = options.input_pattern

    cfg.igvtools.do = options.do_igvtools
    cfg.coverage.do = options.do_bam_coverage
    cfg.mark_duplicates.do = options.do_mark_duplicates
    cfg.RNAseQC.do = options.do_rnaseq_qc


    cfg.fastq_screen.do = options.do_fastq_screen
    # get absolute path
    fastq_screen_conf = os.path.abspath(options.fastq_screen_conf)
    cfg.fastq_screen.conf = fastq_screen_conf


    # copy the fastq_screen.conf input or default file
    import shutil
    import sequana_pipelines.rnaseq
    if options.fastq_screen_conf:
        assert os.path.exists(fastq_screen_conf)
        shutil.copy(fastq_screen_conf, manager.workdir)
        cfg.fastq_screen.do = True
    else:
        # Just copy the default file
        shutil.copy(os.path.join(sequana_pipelines.rnaseq.__path__[0] ,
            "fastq_screen.conf"), manager.workdir)
    

    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown()

if __name__ == "__main__":
    main()
