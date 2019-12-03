This is is the **rnaseq** pipeline from the `Sequana <https://sequana.readthedocs.org>`_ projet

:Overview: TODO 
:Input: A set of Fastq Files and genome reference and annotation.
:Output: MultiQC reports and feature Counts 
:Status: Production
:Citation: Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, JOSS DOI doi:10.21105/joss.00352


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
    sequana_pipelines_rnaseq --input-directory DATAPATH --run-mode local
    sequana_pipelines_rnaseq --input-directory DATAPATH --run-mode slurm

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
- featureCounts

.. image:: https://raw.githubusercontent.com/sequana/sequana_rnaseq/master/sequana_pipelines/rnaseq/dag.png


Details
~~~~~~~~~

This pipeline runs **rnaseq** in parallel on the input fastq files (paired or not). 
A brief sequana summary report is also produced.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the `latest documented configuration file <https://raw.githubusercontent.com/sequana/sequana_rnaseq/master/sequana_pipelines/rnaseq/config.yaml>`_
to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file. 

