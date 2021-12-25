#!/usr/bin/env python3

from viir.params import Params


pm = Params('viir')
args = pm.set_options()


import subprocess as sbp


class ViiR(object):

    def __init__(self, args):
        self.args = args

    def run():

        cmd1 = 'wget https://raw.githubusercontent.com/YuSugihara/ViiR/master/run_viir.sh'

        cmd2 = 'sh ./run_viir.sh {} {} {} {} {} {} {} {}'.format(self.args.out, \
                                                                 self.args.fastq_list, \
                                                                 self.args.pvalue, \
                                                                 self.args.pfam, \
                                                                 self.args.threads, \
                                                                 self.args.max_memory, \
                                                                 self.args.SS_lib_type, \
                                                                 self.args.adapter)

        sbp.run(cmd1,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        sbp.run(cmd2,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)



def main():
    ViiR(args).run()

if __name__ == '__main__':
    main()



