#!/usr/bin/env python3

"""
Similar to seqtailor [PMID:31045209] : reads a VCF file, outputs a genomic
sequence (default length: 31)

Unlike seqtailor, all sequences will have the same length. Moreover, it is
possible to have an absence character (by default the dot ` .` ) for indels.

- When a insertion is larger than `--size` parameter, only first `--size`
  nucleotides are outputed.
- Sequence headers are formated as "<chr>_<position>_<ref>_<alt>".

VCF format spec: https://github.com/samtools/hts-specs/blob/master/VCFv4.4.pdf
"""

import sys
import os
import argparse
import ascii
import pyfaidx
import shutil

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
    # ~ vcf_ok, vcf_msg = input_ctrl(args, chr_dict)
    # ~ if not vcf_ok:
        # ~ sys.exit(f"{COL.RED}{vcf_msg}")
    resp = compute(args, chr_dict)
    output(args, resp)


def _input_ok(args, rows, resp, chr_dict, cols_id):
    for row in rows:
        if row.startswith('#'):
            continue
        try:
            fields = row.rstrip('\n').split('\t')
            nfields = len(fields)
            chr, pos, id, ref, alt = fields[:5]
        except ValueError:
            resp["error"] = ("ErrorVcfFormat: not enough columns for a vcf (expected at least 5).")
            resp["is_ok"] = False
            return False

        ### check if --add-columns is compatible with number of columns
        if args.add_columns and max(cols_id) > nfields:
            resp["error"] = (f"VCF file has {nfields} columns, but you asked for "
                  f"{max(cols_id)}.")
            resp["is_ok"] = False
            return False

        ### Check some commonly issues
        if not pos.isdigit():
            resp["error"] = (f"ErrorVcfFormat: second column is the position. It must be a "
                     f"digit (found: {pos!r}).\n"
                      "A commonly issue is that the header is not commented by a '#' ")
            resp["is_ok"] = False
            return False
        if chr not in chr_dict:
            resp["error"] = ("ErrorChr: Chromosomes are not named in the same way in the "
                      "query and the genome file. Below the first chromosome found: \n"
                     f" your query: {chr}\n"
                     f" genome: {next(iter(chr_dict.keys()))}\n"
                     f"Please, correct your request (or modify the file '{args.genome}.fai').")
            resp["is_ok"] = False
            return False
        break
    return True, "ok"


def compute(args, chr_dict):
    ### object to return
    resp = {
        "is_ok": True,
        "result": [],
        "warning": [],
        "error": None
        }

    ### convert input as list
    if isinstance(args.input, str):
        rows = args.input.splitlines()
    else:
        rows = args.input.read().splitlines()

    ### define generic variables
    res_ref = []
    res_alt = []
    valid_nuc = ["A", "T", "C", "G", args.blank]
    cols_id = ascii.get_index(args.add_columns)    # columns chars are converted as index, ex: AA -> 27

    ### check input syntax
    if not _input_ok(args, rows, resp, chr_dict, cols_id):
        return resp

    ### starts computing
    for i,row in enumerate(rows):
        if not row or row.startswith('#'):
            continue
        fields = row.split('\t')
        chr, position, id, ref, alts = fields[:5]

        alts = alts.split(',')
        for alt in alts:

            header = f"{chr}:{position}_{ref}_{alt}"
            tsv_cols =  '\t' + '\t'.join([chr, position, ref, alt]) if args.output_format == 'tsv' else ''

            ### WARNING: event bigger than kmer size
            if len(ref) > args.size :
                resp["warning"].append(f" REF deletion larger than {args.size} at line {i+1}, truncated in output ({len(ref)} pb).")

            ### ERROR: REF/ALT base is not valid
            bad_nuc = [ a for a in (alt[0], ref[0]) if a not in valid_nuc]
            if bad_nuc:
                bad_nuc = bad_nuc[0]
                resp["warning"].append(f"The base {bad_nuc!r} is not valid at line {i+1}, ignored.\n"
                        f"    You might add the '-b/--blank {bad_nuc}' option or check your VCF file."
                        )
                continue

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

            try:
                ## define REF kmer
                l_ref2  = 0 if ref == args.blank else len(ref)
                l_ref1  = (args.size - l_ref2) // 2
                l_ref3  = l_ref1 + corr_ref
                ps_ref2 = int(position)-1                               # -1 for pyfaidx
                ps_ref1 = ps_ref2 - l_ref1
                pe_ref3 = ps_ref2 + l_ref2 + l_ref3
                ref_seq = str(chr_dict[chr][ps_ref1:pe_ref3])

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
            except:
                resp["warning"].append(f"Warning: something went wrong at line {num_row}, ignored.")
                break

            ### WARNING: REF bases must be the same as the calculated position
            seq_ref2 = chr_dict[chr][ps_ref2:ps_ref2+l_ref2]
            if l_ref2 and not ref == seq_ref2:
                resp["warning"].append(" Mismatch between REF and genome "
                                f"at line {i+1} (chr{chr}:{ps_ref2+1}).\n"
                                f"   - REF in the vcf file: {ref!r}\n"
                                f"   - Found in the genome: '{seq_ref2}'\n"
                                "   Please check if the given genome is appropriate.")
            col_sep = ' ' if args.output_format == 'fa' else '\t'

            ### Append results in lists
            if len(ref_seq) == args.size == len(alt_seq):
                ### append additional selected columns to the header
                added_cols = f"{col_sep}{col_sep.join([fields[num-1] for num in cols_id])}" if cols_id else ''
                ### append to list according of output format
                if args.output_format == "tsv":
                    res_ref.append(f"{ref_seq}{col_sep}{header}_ref{tsv_cols}{col_sep}ref{added_cols}")
                    res_alt.append(f"{alt_seq}{col_sep}{header}_alt{tsv_cols}{col_sep}alt{added_cols}")
                else:
                    res_ref.append(f">{header}_ref{added_cols}")
                    res_ref.append(ref_seq)
                    res_alt.append(f">{header}_alt{added_cols}")
                    res_alt.append(alt_seq)
            elif len(alt_seq) > args.size:
                resp["warning"].append(f" ALT insertion larger than {args.size} at line {i+1}, truncated in output ({len(alt)} bp).")
            else:
                resp["warning"].append(f" Sequence size not correct at line {i+1}, ignored"
                                f"({len(alt_seq)} != {args.size}).")


    # ~ res = list()
    if args.output_format == 'tsv':
        str_cols = '\t' + "col_{}".format('\tcol_'.join(args.add_columns)) if args.add_columns else ''
        resp["result"].append(f"sequence\tid\tchr\tposition\tREF\tALT\ttype{str_cols}")

    if args.type == 'alt':
        resp["result"] += res_alt
    elif args.type == 'ref':
        resp["result"] += res_ref
    else:
        if args.output_format == 'fa':
            for i in range(0, len(res_alt), 2):
                resp["result"] += [res_ref[i], res_ref[i+1]]
                resp["result"] += [res_alt[i], res_alt[i+1]]
        else:
            for i,_ in enumerate(res_alt):
                resp["result"].append(res_ref[i])
                resp["result"].append(res_alt[i])
    return resp



def output(args, resp):
    ### OUTPUT RESULTS
    ext = args.output_format
    ## define output file
    if not args.output:
        name, _ = os.path.splitext(os.path.basename(args.input.name))
        args.output = f"{name}-vcf2seq.{ext}"

    if resp["is_ok"]:
        ## write results in file
        if resp["result"]:
            with open(args.output, 'w') as fh:
                for result in resp["result"]:
                    fh.write(f"{result}\n")
            print(f"\n🧬 {args.output} succefully created.\n")
        ### WARNINGS
        if resp["warning"]:
            print(f"{COL.PURPLE}⚠️  Warnings:")
            for warning in resp["warning"]:
                print(f"  - {warning}")
            print(COL.END)
    else:
        print(f"\n☠️  {COL.RED}{resp['error']}\n")


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
    doc_sep = '=' * min(80, shutil.get_terminal_size()[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument("input",
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
    ### Go to "usage()" without arguments
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return parser.parse_args()


if __name__ == "__main__":
    main()
