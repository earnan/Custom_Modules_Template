#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   M_dna2acid.py
#         Author:   yujie
#    Description:   M_dna2acid.py
#        Version:   1.0
#           Time:   2022/04/18 15:22:18
#  Last Modified:   2022/04/18 15:22:18
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
#from icecream import ic
import argparse
import linecache
import os
import re
import time


parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   根据第n套密码子表翻译DNA')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input', type=str,
                      metavar='[input]', help='str/file', required=False)
optional.add_argument('-n', '--number', type=int,
                      metavar='[密码子表]', help='1-11,默认叶绿体11', default=11, required=False)
optional.add_argument("-l", "--lenth", action="store_true", help="-l触发true")
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def input_format2str(input):
    if os.path.isfile(input):
        with open(input, 'r') as f:
            seq = ''
            for line in f:
                if not line.startswith('>'):
                    seq += line.strip('\n')
        if args.lenth:
            seq = len(seq)
    else:
        seq = input
        if args.lenth:
            seq = len(seq)
    return seq


seq = input_format2str(args.input)
if type(seq) == type(''):
    coding_dna = Seq(seq)
    print(coding_dna.translate(table=args.number))
else:
    if type(seq) == type(1):
        print(seq)
