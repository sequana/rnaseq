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

            For a local run, use :

                sequana_pipelines_rnaseq --input-directory PATH_TO_DATA --run-mode local

            For a run on a SLURM cluster:

                sequana_pipelines_rnaseq --input-directory PATH_TO_DATA --run-mode slurm

        """
        )
        super(Options, self).__init__(usage=usage, prog=prog, description="")

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

        pipeline_group = self.add_argument_group("pipeline_genome")
        pipeline_group.add_argument("--genome-do", dest="genome_do", 
            action="store_true")
        pipeline_group.add_argument("--genome-directory", dest="genome_directory", 
            default=".", required=True)
        pipeline_group.add_argument("--genome-name", dest="genome_name", 
            required=True)

        pipeline_group = self.add_argument_group("pipeline_cutadapt")
        pipeline_group.add_argument("--cutadapt-fwd", dest="cutadapt_fwd", 
            default="")
        pipeline_group.add_argument("--cutadapt-rev", dest="cutadapt_rev", 
            default="")
        pipeline_group.add_argument("--cutadapt-quality", dest="cutadapt_quality", 
            default=30, type=int)
        
"""
cutadapt:
 37     do: true
 38     tool_choice: cutadapt
 39     adapter_choice: ''
 40     design_file: ''
 43     m: 20
 44     mode: b
 45     options: -O 6 --trim-n
 47     threads: 4

fastq_screen:
 49     do: false
 50     conf: fastq_screen.conf
 51     options: --subset 200000  --aligner bowtie2
 52     pf2_report: false

 53 kraken:
 54     do: false
 55     database_directory: ''
 56     threads: 4

 57 bowtie1_mapping_rna:
 58     do: true
 59     options: ''
 60     threads: 4

 61 star_mapping:
 62     do: false
 63     options: --outFilterMismatchNoverLmax 0.05 --seedSearchStartLmax 20
 64     threads: 4

 65 bowtie1_mapping_ref:
 66     do: true
 67     options: --chunkmbs 400 -m 1
 68     threads: 4

 69 bowtie2_mapping:
 70     do: false
 71     options: ''
 72     threads: 4

 73 feature_counts:
 74     do: true
 75     options: -t gene -g ID -s 1
 76     threads: 2

coverage:
 78     do: false
 79     binSize: 10
 80     genomeSize: 2150570000
 81     extendReads: 65
 82     minFragmentLength: 0
 83     maxFragmentLength: 0
 84     threads: 4
 85     chromSize_file: hg19.chrom.sizes

 86 mark_duplicates:
 87     do: false
 88     remove: 'False'
 89     tmpdir: /tmp/
 90     threads: 4

 91 RNAseQC:
 92     do: false
 93     gtf_file: /path/to/your/directory/annotation/gtf/hg19.gtf
 94     BWArRNA_file: /path/to/your/directory/human_rRNA/human_all_rRNA.fasta
 95     options: ''

 96 multiqc:
 97     options: -f -x *_init_* -x *left_kept_reads* -x *report_rnaseq* -e htseq
-e slamdunk
 98     output-directory: multiqc
 99     config_file: multiqc_config.yaml

100 SARTools:
101     do: false
102     design: path/to/file
103     projectName: RNAseq
104     author: NAME
105     rawDir: featureCounts/
106     featuresToRemove: 'NULL'
107     varInt: group
108     condRef: WT
109     batch: 'NULL'
110     fitType: parametric
111     cooksCutoff: 'TRUE'
112     independentFiltering: 'TRUE'
113     alpha: '0.05'
114     pAdjustMethod: BH
115     typeTrans: VST
116     locfunc: median
"""

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
    # EXAMPLE TOREPLACE WITH YOUR NEEDS
    #cfg.TODO = os.path.abspath(options.working_directory)
    #cfg.YOURSECTION.TODO = options.TODO

    cfg.genome.genome_do = options.genome_do
    cfg.genome.genome_directory = os.path.abspath(options.genome_directory)
    cfg.genome.name = options.genome_name
    cfg.genome.fasta_file = options.genome_name + ".fa"
    cfg.genome.gff_file = options.genome_name + ".gff"

    cfg.cutadapt.fwd = options.cutadapt_fwd
    cfg.cutadapt.rev = options.cutadapt_rev
    cfg.cutadapt.quality = options.cutadapt_quality
    
    cfg.input_directory = os.path.abspath(options.input_directory)
    cfg.input_pattern = options.input_pattern


    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown()


if __name__ == "__main__":
    main()
