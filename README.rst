This is is the **rnaseq** pipeline from the `Sequana <https://sequana.readthedocs.org>`_ projet

:Overview: RNASeq analysis from raw data to feature counts
:Input: A set of Fastq Files and genome reference and annotation.
:Output: MultiQC reports and feature Counts
:Status: Production
:Citation(sequana): Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, JOSS DOI doi:10.21105/joss.00352
:Citation(pipeline): 
    .. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4047837.svg
       :target: https://doi.org/10.5281/zenodo.4047837

Installation
~~~~~~~~~~~~

You must install Sequana first::

    pip install sequana

Then, just install this package::

    pip install sequana_rnaseq


Usage
~~~~~

::

    sequana_pipelines_rnaseq --help
    sequana_pipelines_rnaseq --input-directory DATAPATH --genome-directory genome --aligner star

This creates a directory with the pipeline and configuration file. You will then need 
to execute the pipeline::

    cd rnaseq
    sh rnaseq.sh  # for a local run

This launch a snakemake pipeline. If you are familiar with snakemake, you can 
retrieve the pipeline itself and its configuration files and then execute the pipeline yourself with specific parameters::

    snakemake -s rnaseq.rules -c config.yaml --cores 4 --stats stats.txt

Or use `sequanix <https://sequana.readthedocs.io/en/master/sequanix.html>`_ interface.

Requirements
~~~~~~~~~~~~

This pipelines requires the following executable(s):

- bowtie
- bowtie2
- STAR
- featureCounts (subread package)
- picard
- multiqc

More may be needed depending on the configuration file options. For instance,
you may use fastq_screen, in which case you need to install it and configure it. 

.. image:: https://raw.githubusercontent.com/sequana/sequana_rnaseq/master/sequana_pipelines/rnaseq/dag.png


Details
~~~~~~~~~

This pipeline runs **rnaseq** in parallel on the input fastq files (paired or not). 
A brief sequana summary report is also produced.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the `latest documented configuration file <https://raw.githubusercontent.com/sequana/sequana_rnaseq/master/sequana_pipelines/rnaseq/config.yaml>`_
to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file. 


.. warning:: the RNAseQC rule is switch off and is not currently functional in
   version 0.9.X

Changelog
~~~~~~~~~

========= ====================================================================
Version   Description
========= ====================================================================
0.9.18    * Fix typo (regression bug) + add tolerance in schema + generic 
            title in multiqc_config. (oct 2020)
0.9.17    * add the *tolerance* parameter in the feature_counts rule as a user
            parameter (config and pipeline). 
0.9.16    * Best feature_counts is now saved into rnadiff/feature_counts 
            directory and rnadiff scripts have been updated accordingly
          * the most probable feature count option is now computed more
            effectivily and incorporated inside the Snakemake pipeline (not in
            the onsuccess) so that multiqc picks the best one (not the 3 
            results)
          * the target.txt file can be generated inside the pipeline if user
            fill the rnadiff/conditions section in the config file
          * indexing options are filled automatically when calling
            sequana_rnaseq based on the presence/absence of the index 
            of the aligner being used.
          * salmon now integrated and feature counts created (still WIP in
            sequana)
0.9.15    * FastQC on raw data skipped by default (FastQC
            for processed data is still available)
          * Added paired options (-p) for featureCounts
          * Switch back markduplicates to False for now.
0.9.14    * Use only R1 with bowtie1
          * set the memory requirements for mark_duplicates in cluster_config
            file
          * Set temporary directory for mark_duplicates to be local ./tmp
0.9.13    * set mark_duplicate to true by default
          * use new sequana pipeline manager
          * export all features counts in a single file
          * custom HTML report
          * faster --help calls
          * --from-project option added
0.9.12    * include salmon tool as an alternative to star/bowtie2
          * include rnadiff directory with required input for Differential
            analysis
0.9.11    * Automatic guessing of the strandness of the experiment
0.9.10    * Fix multiqc for RNAseQC rule
0.9.9     * Fix RNAseQC rule, which is now available. 
          * Fix ability to use existing rRNA file as input
0.9.8     * Fix indexing for bowtie1 to not be done if aligner is different
          * add new options: --feature-counts-options and --do-rnaseq-qc,
            --rRNA-feature
          * Based on the input GFF, we now check the validity of the rRNA
            feature and feature counts options to check whether the feature 
            exists in the GFF
          * schema is now used to check the config file values
          * add a data test for testing and documentation
0.9.7     * fix typo found in version 0.9.6
0.9.6     * Fixed empty read tag in the configuration file
          * Possiblity to switch off cutadapt section
          * Fixing bowtie2 rule in sequana and update the pipeline accordingly
          * Include a schema file
          * output-directory parameter renamed into output_directory (multiqc 
            section)
          * handle stdout correctly in fastqc, bowtie1, bowtie2 rules
0.9.5     * Fixed https://github.com/sequana/sequana/issues/571
          * More cutadapt commands and sanity checks
          * Fixed bowtie2 options import in rnaseq.rules
0.9.4  
0.9.3     if a fastq_screen.conf is provided, we switch the fastqc_screen 
          section ON automatically
0.9.0     **Major refactorisation.**

          * remove sartools, kraken rules. 
          * Indexing is now optional and can be set in the configuration.
          * Configuration file is simplified  with a general section to enter
            the genome location and aligner. 
          * Fixed rules in  sequana (0.8.0) that were not up-to-date with
            several executables used in the  pipeline including picard,
            fastq_screen, etc. See Sequana Changelog for details with respect
            to rules changes. 
          * Copying the feature counts in main directory  ready to use for 
            a differential analysis.
========= ====================================================================
