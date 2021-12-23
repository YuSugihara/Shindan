#!/usr/bin/env python3

import subprocess as sbp
from shindan.params import Params


pm = Params('shindan')
args = pm.set_options()

class Shindan(object):

    def __init__(self, args):
        self.args = args

    def run():

        cmd1 = 'wget https://raw.githubusercontent.com/YuSugihara/Shindan/master/run_shindan.sh'

        cmd2 = 'sh ./run_shindan.sh {} {} {} {} {} {} {} {}'.format(self.args.out, \
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
    Shindan(args).run()

if __name__ == '__main__':
    main()



