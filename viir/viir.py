#!/usr/bin/env python3

from viir.params import Params


pm = Params('viir')
args = pm.set_options()

import os
import sys
import subprocess as sbp


class ViiR(object):

    def __init__(self, args):
        self.args = args

    def run(self):
        os.mkdir(self.args.out)
        cmd1 = 'wget https://raw.githubusercontent.com/YuSugihara/ViiR/master/run_viir.sh -O ./{}/run_viir.sh'.format(self.args.out)
        cmd2 = 'bash -e ./{0}/run_viir.sh {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(self.args.out, \
                                                                                      self.args.fastq_list, \
                                                                                      self.args.pvalue, \
                                                                                      self.args.pfam, \
                                                                                      self.args.threads, \
                                                                                      self.args.max_memory, \
                                                                                      self.args.SS_lib_type, \
                                                                                      self.args.adapter, \
                                                                                      self.args.blastdb)

        print(cmd1, file=sys.stderr, flush=True)
        sbp.run(cmd1,
                stdout=sys.stdout,
                stderr=sys.stderr,
                shell=True,
                check=True)

        print(cmd2, file=sys.stderr, flush=True)
        sbp.run(cmd2,
                stdout=sys.stdout,
                stderr=sys.stderr,
                shell=True,
                check=True)



def main():
    ViiR(args).run()

if __name__ == '__main__':
    main()



