

import argparse
from Bio import SeqIO
import os
import re

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   read_fasta_to_dic')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[dir]', help='输入gbk所在目录', type=str, default="F:\\ref_tre\\gene\\feature\\Mm_G1.gene.seq", required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def read_fasta_to_dic(infasta):
    f = open(infasta, 'r')
    # fo = open(output_file, 'w')
    seq_id = ''
    seq_index = []
    dict_seq = {}
    # rres = []

    for line in f:
        # 如果是">"开头的，就创建一个key键
        if line.startswith('>'):
            seq_id = line.strip('\n')  # ID为键
            seq_index.append(line.replace(
                "\n", "").replace(">", ""))  # 顺便创建索引的列表
            dict_seq[seq_id] = ''  # 有key无value
        # 如果不是">"开头的，在这个key键下添加value值
        else:
            dict_seq[seq_id] += line.strip('\n')

    print(len(dict_seq))  # 包含说明行序列行的字典结构
    return dict_seq


if __name__ == '__main__':
    dict_seq = read_fasta_to_dic(
        'F:\\ref_tre\\gene\\feature\\Mm_G1.gene.seq')


"""
def readfasta(input_file, output_file):
    fa = open(input_file, 'r')
    fo = open(output_file, 'w')
    seqid = ''
    index = []
    seq = {}
    rres = []

    for line in fa:
        # 如果是">"开头的，就创建一个key键
        if line.startswith('>'):
            seqid = line.strip('\n')  # ID为键
            index.append(line.replace("\n", "").replace(">", ""))  # 顺便创建索引的列表
            seq[seqid] = ''  # 有key无value
        # 如果不是">"开头的，在这个key键下添加value值
        else:
            seq[seqid] += line.strip('\n')

    print(seq)  # 包含说明行序列行的字典结构
    # 把value值提取出，写入输出文件
    for value in seq.values():
        rres.append(value)
        fo.write(value + '\n')
    # 返回 ['TAGACCA','ATACA','GATTACA']
    return rres


# seq_list = readfasta('14_1sample_in.txt', '14_1output.txt')
seq_list = readfasta(
    'F:\\ref_tre\\gene\\feature\\Mm_G1.gene.seq', 'F:\\ref_tre\\gene\\feature\\Mm_G1.gene.out.txt')
"""
