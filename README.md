# vcf2seq


## Aim

Similar to seqtailor [PMID:31045209] : reads a VCF file, outputs a genomic sequence (default length: 31)

Unlike seqtailor, all sequences will have the same length. Moreover, it is possible to have an absence character (by default the dot ` .` ) for indels.

- When a insertion is larger than ``--size`` parameter, only first ``--size`` nucleotides are outputed.
- Sequence headers are formated as "<chr>_<position>_<ref>_<alt>".

VCF format specifications: https://github.com/samtools/hts-specs/blob/master/VCFv4.4.pdf


## Installation

```
pip install vcf2seq
```


## usage

```
usage: vcf2seq.py [-h] -g genome [-s SIZE] [-t {alt,ref,both}] [-b BLANK] [-a ADD_COLUMNS [ADD_COLUMNS ...]] [-o OUTPUT] [-v] vcf


positional arguments:
  vcf                   vcf file (mandatory)

options:
  -h, --help            show this help message and exit
  -g, --genome genome   genome as fasta file (mandatory)
  -s, --size SIZE       size of the output sequence (default: 31)
  -t, --type {alt,ref,both}
                        alt, ref, or both output? (default: alt)
  -b, --blank BLANK     Nucleotide absence character, default is dot (.)
  -a, --add-columns ADD_COLUMNS [ADD_COLUMNS ...]
                        Add one or more columns to header (ex: '-a 3 AA' will 
                        add columns 3 and 27). The first column is '1' (or 'A')
  -d, --delimiter DELIMITER
                        with -a/--add-columns and a fasta format output, 
                        specifies a delimiter (default: '|')
  -o, --output OUTPUT   Output file (default: <input_file>-vcf2seq.fa/tsv)
  -f, --output-format {fa,tsv}
                        Output file format (default: fa)
  -v, --version         show program's version number and exit
```
