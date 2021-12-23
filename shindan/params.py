import os
import sys
import argparse
from multiprocessing import Pool
import multiprocessing as multi
from shindan.utils import time_stamp
from shindan.__init__ import __version__


class Params(object):

    def __init__(self, program_name):
        self.program_name = program_name

    def set_options(self):
        if self.program_name == 'shindan':
            parser = self.shindan_options()

        if len(sys.argv) == 1:
            args = parser.parse_args(['-h'])
        else:
            args = parser.parse_args()
        return args

    def shindan_options(self):
        parser = argparse.ArgumentParser(description='Vid-kit version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('shindan -l <FASTQ_LIST> \n'
                        '       -o <OUT_DIR> \n'
                        '       [-t <INT>]')

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
                            default=2,
                            type=int,
                            help='Number of threads.',
                            metavar='')

        # set version
        parser.add_argument('-v',
                            '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    