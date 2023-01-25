import sys
import argparse
from viir.__init__ import __version__


class Params(object):

    def __init__(self, program_name):
        self.program_name = program_name

    def set_options(self):
        if self.program_name == 'viir':
            parser = self.viir_options()

        if len(sys.argv) == 1:
            args = parser.parse_args(['-h'])
        else:
            args = parser.parse_args()
        return args

    def viir_options(self):
        parser = argparse.ArgumentParser(description='ViiR version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = 'viir -l <FASTQ_LIST> -o <OUT_DIR> [-t <INT>]'

        # set options
        parser.add_argument('-l',
                            '--fastq-list',
                            action='store',
                            required=True,
                            type=str,
                            help='Fastq list.',
                            metavar='')

        parser.add_argument('-o',
                            '--out',
                            action='store',
                            required=True,
                            type=str,
                            help=('Output directory. Specified name must not\n'
                                  'exist.'),
                            metavar='')

        parser.add_argument('-t',
                            '--threads',
                            action='store',
                            default=16,
                            type=int,
                            help='Number of threads. [16]',
                            metavar='')

        parser.add_argument('-a',
                            '--adapter',
                            action='store',
                            default="Default_adapter",
                            type=str,
                            help=("FASTA of adapter sequences. If you don't\n"
                                  'specify this option, the defaul adapter set\n'
                                  'will be used.'),
                            metavar='')

        parser.add_argument('--pfam',
                            action='store',
                            default="Default_list",
                            type=str,
                            help=("List of Pfam IDs. If you don't specify\n"
                                  'this option, the defaul list will be used.'),
                            metavar='')

        parser.add_argument('--SS-lib-type',
                            action='store',
                            default='No',
                            type=str,
                            help=('Type of strand specific library (No/FR/RF). [No]'),
                            metavar='')

        parser.add_argument('--blastndb',
                            action='store',
                            default='Default_db',
                            type=str,
                            help=('FASTA to annotate your trinity assembly.'),
                            metavar='')

        parser.add_argument('--pvalue',
                            action='store',
                            default=0.01,
                            type=float,
                            help='Threshold of pvalue in DESeq2. [0.01]',
                            metavar='')

        parser.add_argument('--max-memory',
                            action='store',
                            default='32G',
                            type=str,
                            help=('Max memory used in Trinity. [32G]'),
                            metavar='')

        # set version
        parser.add_argument('-v',
                            '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    