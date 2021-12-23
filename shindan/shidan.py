#!/usr/bin/env python3

from shindan.params import Params


pm = Params('shindan')
args = pm.set_options()


import os
import sys
import glob
import subprocess as sbp

class Shindan(object):

    def __init__(self, args):
        self.args = args

    def run():



def main():
    Shindan(args).run()

if __name__ == '__main__':
    main()



