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

# 两种方法,自己写,调包


def trans2acid(codon):  # 翻译成氨基酸
    """
    The Bacterial and Plant Plastid Code (11):
    Stnd    AAs = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    This    AAs = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Starts      = ---M---------------M------------MMMM---------------M------------
    Base1       = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2       = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
    Base3       = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
    """
    genetic_code_number = 11
    acid = ''
    code_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                  'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', }
    acid = code_table[codon]
    return acid


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
