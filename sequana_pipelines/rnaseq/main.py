#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os
import shutil
import subprocess
import sys

import click_completion
import rich_click as click
from sequana_pipetools import SequanaManager
from sequana_pipetools.options import *

click_completion.init()
NAME = "rnaseq"


help = init_click(
    NAME,
    groups={
        "Pipeline Specific": [
            "--aligner",
            "--contaminant-file",
            "--do-igvtools",
            "--do-bam-coverage",
            "--do-mark-duplicates",
            "--do-rnaseqc",
            "--do-rseqc",
            "--genome-directory",
            "--rnaseqc-gtf-file",
            "--rRNA-feature",
            "--rseqc-bed-file",
            "--skip-rRNA",
            "--skip-gff-check",
            "--trimming-quality",
        ],
    },
)


@click.command(context_settings=help)
@include_options_from(ClickInputOptions)
@include_options_from(ClickSnakemakeOptions, working_directory=NAME)
@include_options_from(ClickSlurmOptions)
@include_options_from(ClickGeneralOptions)
@include_options_from(ClickTrimmingOptions)
@include_options_from(ClickFeatureCountsOptions)
@click.option(
    "--genome-directory",
    "genome_directory",
    default=".",
    show_default=True,
    type=click.Path(dir_okay=True, file_okay=False),
    required=True,
)
@click.option(
    "--aligner",
    "aligner",
    required=True,
    type=click.Choice(["bowtie2", "bowtie1", "star", "salmon"]),
    help="a mapper in bowtie, bowtie2, star",
)
@click.option(
    "--rRNA-feature",
    "rRNA",
    help="""Feature name corresponding to the rRNA to be identified in
the input GFF/GTF files. Must exist and be valid. If you do not have any,
you may skip this step using --skip-rRNA or provide a fasta file using --contaminant-file""",
)
@click.option(
    "--skip-rRNA",
    "skip_rRNA",
    is_flag=True,
    help="""skip the mapping on rRNA feature. ignored if --contaminant-file is provided""",
)
@click.option(
    "--contaminant-file",
    default=None,
    show_default=True,
    help="""A fasta file. If used, the rRNA-feature is not used
This option is useful if you have a dedicated list of rRNA feature or a dedicated
fasta file to search for contaminants""",
)
@click.option(
    "--skip-gff-check",
    is_flag=True,
    default=False,
    show_default=True,
    help="""By default we check the coherence between the input
GFF file and related options (e.g. --feature_counts_feature_type and
--feature_counts_attribute options). This may take time e.g. for mouse or human.
Using this option skips the sanity checks""",
)
@click.option(
    "--do-igvtools",
    is_flag=True,
    help="""if set, this will compute TDF files that can be imported in
IGV browser. TDF file allows to quickly visualise the coverage of the mapped
reads.""",
)
@click.option(
    "--do-bam-coverage",
    is_flag=True,
    help="Similar to --do-igvtools using bigwig",
)
@click.option(
    "--do-mark-duplicates",
    is_flag=True,
    help="""Mark duplicates. To be used e.g. with QCs""",
)
@click.option("--do-rnaseqc", is_flag=True, help="do RNA-seq QC using RNAseQC v2")
@click.option(
    "--rnaseqc-gtf-file",
    help="""The GTF file to be used for RNAseQC. Without a valid GTF,
    RNAseqQC will not work. You may try sequana gff-to-gtf application.""",
)
@click.option(
    "--do-rseqc",
    is_flag=True,
    help="""do RNA-seq QC using RseQC. This will need a BED file
corresponding to your GFF file. For prokaryotes, the BED file is created on the
fly.""",
)
@click.option("--rseqc-bed-file", help="""The rseQC input bed file.""")
def main(**options):

    if options["from_project"]:
        click.echo("--from-project Not yet implemented")
        sys.exit(1)
    # the real stuff is here
    manager = SequanaManager(options, NAME)
    manager.setup()

    # aliases
    options = manager.options
    cfg = manager.config.config

    from sequana_pipetools import logger

    logger.setLevel(options.level)

    manager.fill_data_options()
    # --------------------------------------------------------- general
    cfg.general.genome_directory = os.path.abspath(options.genome_directory)
    cfg.general.aligner = options.aligner

    # genome name = cfg.genome.genome_directory
    genome_name = cfg.general.genome_directory.rsplit("/", 1)[1]
    prefix = cfg.general.genome_directory
    fasta = cfg.general.genome_directory + f"/{genome_name}.fa"
    if os.path.exists(fasta) is False:
        logger.critical(
            """Could not find {}. You must have the genome sequence in fasta with the extension .fa named after the genome directory.""".format(
                fasta
            )
        )
        sys.exit()

    # mutually exclusive options
    if options.contaminant_file:
        cfg.general.contaminant_file = os.path.abspath(options.contaminant_file)
        logger.warning("You are using a custom FASTA --contaminant_file so --rRNA-feature will be ignored")
        cfg.general.rRNA_feature = None
    elif options.skip_rRNA:
        cfg.general.rRNA_feature = None
    else:
        cfg.general.rRNA_feature = options.rRNA

    # --------------------------------------------------------- trimming
    cfg.trimming.software_choice = options.trimming_software_choice
    cfg.trimming.do = not options.disable_trimming
    qual = options.trimming_quality

    if options.trimming_software_choice in ["cutadapt", "atropos"]:
        cfg.cutadapt.tool_choice = options.trimming_software_choice
        cfg.cutadapt.fwd = options.trimming_adapter_read1
        cfg.cutadapt.rev = options.trimming_adapter_read2
        cfg.cutadapt.m = options.trimming_minimum_length
        cfg.cutadapt.mode = options.trimming_cutadapt_mode
        cfg.cutadapt.options = options.trimming_cutadapt_options  # trim Ns -O 6
        cfg.cutadapt.quality = 30 if qual == -1 else qual
    else:
        cfg.fastp.minimum_length = options.trimming_minimum_length
        cfg.fastp.quality = 15 if qual == -1 else qual
        cfg.fastp.adapters = ""
        if options.trimming_adapter_read1:
            cfg.fastp.adapters += f"--adapter_sequence {options.trimming_adapter_read1}"
        if options.trimming_adapter_read2:
            cfg.fastp.adapters += f"--adapter_sequence_r2 {options.trimming_adapter_read2}"

        cfg.fastp.options = " --cut_tail "
        cfg.fastp.disable_quality_filtering = False
        cfg.fastp.disable_adapter_trimming = False

    # ----------------------------------------------------- feature counts
    cfg.feature_counts.options = options.feature_counts_options
    cfg.feature_counts.strandness = options.feature_counts_strandness
    cfg.feature_counts.attribute = options.feature_counts_attribute
    cfg.feature_counts.feature = options.feature_counts_feature_type
    cfg.feature_counts.extra_attributes = options.feature_counts_extra_attributes

    # ------------------------------------------------------ optional
    cfg.igvtools.do = options.do_igvtools
    cfg.bam_coverage.do = options.do_bam_coverage
    cfg.mark_duplicates.do = False
    if options.do_mark_duplicates:
        cfg.mark_duplicates.do = True

    # -------------------------------------------------------- RNAseqQC
    cfg.rnaseqc.do = options.do_rnaseqc

    if options.do_rnaseqc:
        if options.rnaseqc_gtf_file is None:
            logger.info(
                "You asked for RNA_seqc QC assessements but no GTF"
                " file provided; Please use --rnaseqc-gtf-file option. Switching off in your"
                " config file and continuing. You may use 'sequana gff2gtf input.gff' to create"
                " the gtf file"
            )
            cfg.rnaseqc.do = False
        if options.aligner in ["salmon"]:
            logger.info(
                "WARNING"
                "You asked for RNA_seqc QC assessements but no"
                " BAM will be generated by the salmon aligner. Switching off this option. "
            )
            cfg.rnaseqc.do = False

    cfg.rnaseqc.gtf_file = options.rnaseqc_gtf_file

    cfg.rseqc.do = options.do_rseqc
    cfg.rseqc.bed_file = options.rseqc_bed_file

    # -------------------------------------------------------- RNAdiff

    # import sequana_pipelines.rnaseq

    # SANITY CHECKS
    # -------------------------------------- do we find rRNA feature in the GFF ?
    # if we do not build a custom feature_counts set of options, no need to
    # check carefully the GFF; if users knows what he is doing; no need to
    # check the GFF either
    if options.skip_gff_check is False and "," not in cfg.feature_counts.feature:
        logger.info("Checking your input GFF file and rRNA feature if provided")

        from sequana.gff3 import GFF3

        genome_directory = os.path.abspath(cfg.general.genome_directory)
        genome_name = genome_directory.rsplit("/", 1)[1]
        prefix_name = genome_directory + "/" + genome_name
        gff_file = prefix_name + ".gff"

        gff = GFF3(gff_file)
        df_gff = gff.df  # This takes one minute on eukaryotes. No need to
        valid_features = gff.features  # about 3 seconds
        valid_attributes = gff.attributes  # about 10 seconds

        # first check the rRNA feature
        if cfg["general"]["rRNA_feature"] and cfg["general"]["rRNA_feature"] not in valid_features:

            logger.error(
                "rRNA feature not found in the input GFF ({})".format(gff_file)
                + " This is probably an error. Please check the GFF content and /or"
                " change the feature name with --rRNA-feature based on the content"
                " of your GFF. Valid features are: {}".format(valid_features)
            )
            sys.exit()

        # then, check the main feature
        fc_type = cfg.feature_counts.feature
        fc_attr = cfg.feature_counts.attribute

        logger.info("Checking your input GFF file and feature counts options.")
        logger.info(f"You chose '{fc_type}' feature and '{fc_attr}' attribute")
        # if only one feature (99% of the projet)
        if "," not in fc_type:
            fc_types = [fc_type]

        for fc_type in fc_types:
            S = sum(df_gff["genetic_type"] == fc_type)
            if S == 0:
                logger.error(
                    "Found 0 entries for feature '{}'. Please choose a valid feature from: {}".format(
                        fc_type, valid_features
                    )
                )
                sys.exit()
            else:
                logger.info(f"Found {S} '{fc_type}' entries")

            # now we check the attribute:
            dd = df_gff.query("genetic_type==@fc_type")
            attributes = [y for x in dd.attributes for y in x.keys()]
            S = attributes.count(fc_attr)
            if S == 0:
                uniq_attributes = set(attributes)
                logger.error(
                    f"Found 0 entries for attribute '{fc_attr}'. Please choose a valid attribute from: {uniq_attributes}"
                )
                sys.exit()
            else:
                unique = set([x[fc_attr] for k, x in dd.attributes.items() if fc_attr in x])
                logger.info(f"Found {S} '{fc_attr}' entries for the attribute [{len(unique)} unique entries]")

            if S != len(unique):
                logger.warning("Attribute non-unique. Feature counts should handle it")

            if options.feature_counts_extra_attributes:
                for extra_attr in cfg.feature_counts.extra_attributes.split(","):
                    if extra_attr not in set(attributes):
                        logger.error("{extra_attr} not found in the GFF attributes. Try one of {set(attributes)}")
                        sys.exit()

    # need to move the custom file into the working directoty
    if "," in cfg.feature_counts.feature:
        logger.info("Building a custom GFF file (custom.gff) using Sequana. Please wait")
        genome_directory = os.path.abspath(cfg.general.genome_directory)
        genome_name = genome_directory.rsplit("/", 1)[1]
        prefix_name = genome_directory + "/" + genome_name

        from sequana import GFF3

        gff = GFF3(prefix_name + ".gff")
        fc_types = cfg.feature_counts.feature.strip().split(",")
        gff.save_gff_filtered(features=fc_types, filename="custom.gff")
        cfg.general.custom_gff = "custom.gff"

    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown()

    try:  # option added in latest version
        if cfg.general.custom_gff:
            shutil.copy(cfg.general.custom_gff, options.workdir)
    except:
        pass


if __name__ == "__main__":
    main()
