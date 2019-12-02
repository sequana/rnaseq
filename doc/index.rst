Sequana rnaseq pipeline documentation
#####################################################

|version|, |today|, status: production

The **rnaseq** pipeline is a `Sequana <https://github.com/sequana/sequana>`_ pipeline. You can find the source code 
on  `https://github.com/sequana/sequana_rnaseq <https://github.com/sequana/sequana_rnaseq/>`_. Would you have issues
about the code, usage or lack of information, please fill a report 
on `**Sequana** itself <https://github.com/sequana/sequana/issues>`_ indicating the pipeline name (We centralized all 
pipelines issues on **Sequana** repository only so as to be more responsive).

If you use **Sequana**, please do not forget to cite us:

    Cokelaer et al, (2017), 'Sequana': a Set of Snakemake NGS pipelines, Journal of
    Open Source Software, 2(16), 352, `JOSS DOI doi:10.21105/joss.00352 <http://www.doi2bib.org/bib/10.21105/joss.00352>`_

.. contents::
   :depth: 2

What is Sequana ?
=====================

**Sequana** is a versatile tool that provides

#. A Python library dedicated to NGS analysis (e.g., tools to visualise standard NGS formats).
#. A set of Pipelines dedicated to NGS in the form of Snakefiles
#. Standalone applications
    #. sequana_coverage ease the
       extraction of genomic regions of interest and genome coverage information
    #. sequana_taxonomy performs a quick
       taxonomy of your FastQ. This requires dedicated databases to be downloaded.
    #. Sequanix, a GUI for Snakemake workflows (hence Sequana pipelines as well)

To join the project, please let us know on `github <https://github.com/sequana/sequana/issues/306>`_.

For more information, please see `github <https://sequana.readthedocs.io>`_.

The Sequana rnaseq pipeline
==============================================

.. include:: ../README.rst


Rules
=====

.. snakemake:: fastqc
.. snakemake:: multiqc

