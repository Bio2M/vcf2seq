# vcf2seq


## Aim

like seqTailor, give a VCF file, it return genomic sequence.

## Installation

```
pip install vcf2ses
```


## usage

```
usage: vcf2seq.py [-h] -g genome [-s SIZE] [-t {alt,ref,both}] [-b BLANK] [-a ADD_COLUMNS [ADD_COLUMNS ...]] [-o OUTPUT] [-v] vcf

================================================================================
Like seqtailor, give a VCF file, it return genomic sequence (default length: 31)

Nota:
- When a insertion is larger than '--size' option, only first '--size' nucleotides are outputed.
- header ID are formated like "<chr>_<position>_<ref>_<alt>".

VCF format specifications: https://github.com/samtools/hts-specs/blob/master/VCFv4.4.pdf
================================================================================

positional arguments:
  vcf                   vcf file (mandatory)

options:
  -h, --help            show this help message and exit
  -g genome, --genome genome
                        genome as fasta file (mandatory)
  -s SIZE, --size SIZE  size of the output sequence (defalt: 31)
  -t {alt,ref,both}, --type {alt,ref,both}
                        alt, ref, or both output? (default: alt)
  -b BLANK, --blank BLANK
                        Missing nucleotide character, default is dot (.)
  -a ADD_COLUMNS [ADD_COLUMNS ...], --add-columns ADD_COLUMNS [ADD_COLUMNS ...]
                        Add one or more columns to header (ex: '-a 3 AA' will add columns 3 and 27). The first column is '1' (or 'A')
  -o OUTPUT, --output OUTPUT
                        Output file (default: <input_file>-vcf2seq.fa)
  -v, --version         show program's version number and exit
```