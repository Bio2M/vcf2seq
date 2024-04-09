#!/usr/bin/env python3

"""
input: characters list
output: position of each characters in ascii uppercase alphabet (1 based)

ex:
>>> ascii_index.py A Z AA 44
>>> [1, 26, 27, 44]
"""

import sys
import argparse
import string


__appname__   = "ascii"
__shortdesc__ = "Give column from sheetspread (ex: 'A BC'), it return indexes (ex: '1 55')"
__licence__   = "none"
__version__   = "0.1.1"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"


def main():
    """ Function doc """
    values = sys.argv[1:]
    indexes = get_index(values)
    print(*indexes)


def get_index(values):
    if not values:
        return []
    ### if values is a string
    if isinstance(values, str):
        values = [values]

    indexes = []
    for chars in values:
        ### valeurs must be in ASCII Alphabet
        if not chars.isascii():
            return(f"Error: character {chars} is not in ASCII alphabet")

        if chars.isdecimal():
            indexes.append(int(chars))
        else:
            idx = 0
            for i,char in enumerate(chars[::-1]):
                try:
                    idx += 26**i * (string.ascii_uppercase.index(char.upper()) + 1)
                except ValueError as err:
                    sys.exit(f"Error: {err}")
            indexes.append(idx)

    return indexes


def usage():
    """
    Help function with argument parser.
    https://docs.python.org/3/howto/argparse.html?highlight=argparse
    """
    doc_sep = '=' * min(72, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    ### OPTION
    parser.add_argument("file",
                        help="file or stdin",
                        type=argparse.FileType('r'),
                        nargs='?' if sys.stdin else 1,
                        default=sys.stdin,
                        metavar=('file1 file2...'),
                       )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{__version__}",
                       )
    ### Go to "usage()" without arguments or stdin
    if len(sys.argv) == 1 and sys.stdin.isatty():
        parser.print_help()
        sys.exit()
    return parser.parse_args()


if __name__ == "__main__":
    main()
