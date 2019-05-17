PGPY
====

PgPy is a python library **under development** designed for population genomic analysis. PgPy is written using python and interacting with vcf files with [pysam library](http://pysam.readthedocs.io/en/latest/api.html), allowing to quickly iterate through whole genomic data. The current release is developped under python 3.5 and no support is provided for 2.x version. 

Installing PgPy
---------------

Using pip

```bash
pip install git+https://github.com/jsgounot/PgPy.git
```

Or download / clone the github

```bash
git clone https://github.com/jsgounot/PgPy.git
cd PgPy
python setup.py install --user
```

Input dataset
-------------

The main purpose of this library is to work with a merged vcf file based on multiple samples sequencing data. The final vcf file must be tabulated using [tabix](http://www.htslib.org/doc/tabix.html). Merging multiple vcfs into on single vcf can be done using [vcftools](http://vcftools.sourceforge.net/)'s `vcfmerge` function. Since pysam works well with compressed file, you should use bgzip from tabix as well at the end. If you want to work with snpEff results, do not forget to annotate your merged vcf files during the process.


