#!/usr/bin/env python3

"""
Similar to seqtailor [PMID:31045209] : reads a VCF file, outputs a genomic sequence
(default length: 31)

Unlike seqtailor, all sequences will have the same length. Moreover, it is possible to have an
absence character (by default the dot ` .` ) for indels.

- When a insertion is larger than `--size` parameter, only first `--size` nucleotides are outputed.
- Sequence headers are formated as "<chr>_<position>_<ref>_<alt>".

VCF format specifications: https://github.com/samtools/hts-specs/blob/master/VCFv4.4.pdf
"""

import sys
import os
import argparse
import ascii
import pyfaidx

import info


def main():
    """ Function doc """
    args = usage()
    try:
        chr_dict = pyfaidx.Fasta(args.genome) # if fai file doesn't exists, it will be automatically created
    except pyfaidx.FastaNotFoundError as err:
        sys.exit(f"FastaNotFoundError: {err}")
    except OSError as err:
        sys.exit(f"\n{COL.RED}WriteError: directory {os.path.dirname(args.genome)!r} may not be "
                  "writable.\nIf you can't change the rights, you can create a symlink and target "
                  f"it. For example:\n  ln -s {args.genome} $HOME\n{COL.END}")
    vcf_ctrl(args, chr_dict)
    results, warnings = compute(args, chr_dict)
    write(args, results, warnings)


def vcf_ctrl(args, chr_dict):
    with open(args.vcf.name) as fh:
        for row in fh:
            if row.startswith('#'):
                continue
            chr, pos, id, ref, alt, *rest = row.rstrip('\n').split('\t')

            ### Check some commonly issues
            if not pos.isdigit():
                sys.exit(f"{COL.RED}ErrorVcfFormat: second column is the position. It must be a "
                         f"digit (found: {pos!r}).\n"
                          "A commonly issue is that headers are not commented by a '#' ")
            if chr not in chr_dict:
                sys.exit(f"{COL.RED}ErrorChr: Chromosomes are not named in the same way in the "
                          "query and the genome file. Below the first chromosome found: \n"
                         f" your query: {chr}\n"
                         f" genome: {next(iter(chr_dict.keys()))}\n"
                         f"Please, correct your request (or modify the file '{args.genome}.fai').")
            break


def compute(args, chr_dict):
    res_ref = []
    res_alt = []
    warnings = []
    num_row = 0
    cols_id = ascii.get_index(args.add_columns)    # columns chars are converted as index, ex: AA -> 27
    for variant in args.vcf:
        num_row += 1
        if not variant.rstrip('\n') or variant.startswith('#'):
            continue
        fields = variant.rstrip('\n').split('\t')
        chr, position, id, ref, alts = fields[:5]

        ### check if --add-columns is compatible with number of columns
        if args.add_columns and max(cols_id) > len(fields):
            sys.exit(f"\n{COL.RED}Error: vcf file has {len(fields)} columns, but you asked for "
                  f"{max(args.add_columns)}.{COL.END}\n")

        alts = alts.split(',')
        for alt in alts:

            header = f"{chr}_{position}_{ref}_{alt}"

            ### WARNING: event bigger than kmer size
            if max(len(ref), len(alt)) > args.size :
                warnings.append(f"Warning: large alteration ({chr}_{position}_{ref}_{alt}), truncated in output.")


            ### ERROR: REF base is not valid
            NUC = ["A", "T", "C", "G", args.blank]
            for nuc in (ref[0], alt[0]):
                if nuc not in NUC:
                    sys.exit(f"{COL.RED}Error: REF base {nuc!r} is not valid "
                            f"(line {num_row}).\n       You might add the '--blank {nuc}' "
                            "option or check your VCF file.")

            #####################################################################################
            #                Some explanations on variable naming                               #
            #                                                                                   #
            #  l = length                                                                       #
            #  ps = position start                                                              #
            #  pe = position end                                                                #
            #                                                                                   #
            #                  ps_ref2: position of the first base of REF                       #
            #                  |                                                                #
            #        l_ref1    | l_ref2     l_ref3                                              #
            #   |--------------|---------|--------------|                                       #
            #        l_alt1     l_alt2       l_alt3                                             #
            #   |------------|-------------|------------|                                       #
            #  ps_alt1                                 pe_alt3                                  #
            #                                                                                   #
            #####################################################################################

            ### define some corrections
            corr_ref = 0
            corr_alt = 0
            if not args.size&1:                                 # k is pair
                if len(ref)&1 and ref != args.blank: corr_ref += 1  # corr_ref + 1 if REF length is unpair
                if len(alt)&1 and alt != args.blank: corr_alt += 1  # corr_alt + 1 if ALT length is unpair
            else:                                               # k is unpair
                if not len(ref)&1: corr_ref += 1                    # corr_ref + 1 if REF length is pair
                if not len(alt)&1: corr_alt += 1                    # corr_alt + 1 if ALT length is pair
                if ref == args.blank: corr_ref += 1                 # missing value for REF
                if alt == args.blank: corr_alt += 1                 # missing value for ALT

            ## define REF kmer
            l_ref2  = 0 if ref == args.blank else len(ref)
            l_ref1  = (args.size - l_ref2) // 2
            l_ref3  = l_ref1 + corr_ref
            ps_ref2 = int(position)-1                               # -1 for pyfaidx
            ps_ref1 = ps_ref2 - l_ref1
            pe_ref3 = ps_ref2 + l_ref2 + l_ref3
            ref_seq = chr_dict[chr][ps_ref1:pe_ref3]

            ## define ALT kmer
            l_alt2 = 0 if alt == args.blank else len(alt)
            l_alt1 = (args.size - l_alt2) // 2
            l_alt3 = (args.size - l_alt2) // 2 + corr_alt
            ps_alt2 = ps_ref2 - (l_alt2 - l_ref2) // 2
            ps_alt1 = ps_ref2 - l_alt1
            ps_alt3 = ps_ref2 + l_ref2
            pe_alt3 = ps_alt3 + l_alt3
            seq_alt1 = chr_dict[chr][ps_alt1:ps_ref2]
            alt = alt if alt != args.blank else ""
            seq_alt3 = chr_dict[chr][ps_alt3:pe_alt3]
            alt_seq = f"{seq_alt1}{alt}{seq_alt3}"

            ### WARNING: REF bases must be the same as the calculated position
            seq_ref2 = chr_dict[chr][ps_ref2:ps_ref2+l_ref2]
            if l_ref2 and not ref == seq_ref2:
                warnings.append("Warning: mismatch between REF and genome "
                                f"at line {num_row} ({chr}:{ps_ref2+1}).\n"
                                f"    - REF on the vcf file: {ref!r}\n"
                                f"    - Found on the genome: '{seq_ref2}'\n"
                                "    Please check if the given genome is appropriate.\n")
            col_sep = ' ' if args.output_format == 'fa' else '\t'
            col_txt =  col_sep.join([fields[num-1] for num in cols_id])
            if len(ref_seq) == args.size == len(alt_seq):
                res_ref.append(f"{header}_ref{col_sep}{col_txt}\n{ref_seq}")
                res_alt.append(f"{header}_alt{col_sep}{col_txt}\n{alt_seq}")
            elif len(alt_seq) > args.size:
                warnings.append(f"Warning: ALT length ({len(alt_seq)} bp) larger than sequence ({args.size} bp) at line {num_row}, ignored.")
            else:
                warnings.append(f"Warning: sequence size not correct at line {num_row}, ignored ({len(alt_seq)} != {args.size}).")

    if args.type == 'alt':
        res = res_alt
    elif args.type == 'ref':
        res = res_ref
    else:
        res = list()
        for i,_ in enumerate(res_alt):
            res.append(res_ref[i])
            res.append(res_alt[i])
    return res, warnings



def write(args, results, warnings):
    ### RESULTS
    ext = args.output_format
    ## define output file
    if not args.output:
        name, _ = os.path.splitext(os.path.basename(args.vcf.name))
        args.output = f"{name}-vcf2seq.{ext}"
    ## write results in file
    with open(args.output, 'w') as fh:
        if not results:
            return
        if ext == 'fa':
            fh.write(">" + '\n>'.join([a for a in results]) + "\n")
        else:
            for row in results:
                header, seq = row.split('\n')
                fh.write(f"{seq}\t")
                fh.write(f"{header}\n")


    ### WARNINGS
    if warnings:
        print(COL.PURPLE)
        print(*warnings, sep='\n')
        print(COL.END)


class COL:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def usage():
    doc_sep = '=' * min(80, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument("vcf",
                        help="vcf file (mandatory)",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("-g", "--genome",
                        help="genome as fasta file (mandatory)",
                        metavar="genome",
                        required=True,
                       )
    parser.add_argument('-s', '--size',
                        type=int,
                        help="size of the output sequence (default: 31)",
                        default=31,
                       )
    parser.add_argument("-t", "--type",
                        type=str,
                        choices=['alt', 'ref', 'both'],
                        default='alt',
                        help="alt, ref, or both output? (default: alt)"
                        )
    parser.add_argument("-b", "--blank",
                        type=str,
                        help="Missing nucleotide character, default is dot (.)",
                        default='.',
                        )
    parser.add_argument("-a", "--add-columns",
                        help="Add one or more columns to header (ex: '-a 3 AA' will add columns "
                             "3 and 27). The first column is '1' (or 'A')",
                        nargs= '+',
                        )
    parser.add_argument("-o", "--output",
                        type=str,
                        help=f"Output file (default: <input_file>-{info.APPNAME}.fa/tsv)",
                        )
    parser.add_argument("-f", "--output-format",
                        choices=['fa', 'tsv'],
                        default='fa',
                        help=f"Output file format (default: fa)",
                        )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{info.VERSION}",
                       )
    ### Go to "usage()" without arguments or stdin
    if len(sys.argv) == 1 and sys.stdin.isatty():
        parser.print_help()
        sys.exit()
    return parser.parse_args()


if __name__ == "__main__":
    main()
