#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   M_read_fasta_to_dic.py
#         Author:   yujie
#    Description:   M_read_fasta_to_dic.py
#        Version:   1.0
#           Time:   2022/03/11 17:45:37
#  Last Modified:   2022/03/11 17:45:37
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
import argparse
from Bio import SeqIO
import os
import re

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   read_fasta_to_dic')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--infile',
                      metavar='[file]', help='file-path', type=str, default="F:\\ref_tre\\gene\\feature\\Mm_G1.gene.seq", required=False)
optional.add_argument('-o', '--outfile',
                      metavar='[file]', help='file-path', type=str, default="F:\\ref_tre\\gene\\feature\\Mm_G1_seq_log", required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def read_fasta_to_dic1(infasta):  # 最简单,针对普通fasta文件 >物种名
    with open(infasta, 'r') as f:
        seq_id = ''
        id_index = []
        dict_seq = {}
        for line in f:
            # 如果是">"开头的，就创建一个key键
            if line.startswith('>'):
                seq_id = line.strip('\n')  # ID为键!!!!!!!!注意这里的键 有 >号
                id_index.append(line.replace(
                    "\n", "").replace(">", ""))  # 顺便创建索引的列表
                dict_seq[seq_id] = ''  # 有key无value
            # 如果不是">"开头的，在这个key键下添加value值
            else:
                dict_seq[seq_id] += line.strip('\n')
    print(len(dict_seq))  # 包含说明行序列行的字典结构
    return dict_seq


def read_fasta_to_dic2(infasta, outfile):
    with open(infasta, 'r') as fi, open(outfile, 'w') as fo:
        seq_id = ''
        id_index = []
        dict_seq = {}
        list_seq = []
        for line in fi:
            # 如果是">"开头的，就创建一个key键
            if line.startswith('>'):
                seq_id = line.strip('\n')  # ID为键
                id_index.append(line.replace(
                    "\n", "").replace(">", ""))  # 顺便创建索引的列表
                dict_seq[seq_id] = ''  # 有key无value
            # 如果不是">"开头的，在这个key键下添加value值
            else:
                dict_seq[seq_id] += line.strip('\n')
        print(len(dict_seq))
        # 把value值提取出，写入输出文件
        for value in dict_seq.values():
            list_seq.append(value)
            fo.write(value + '\n')
        # 返回 ['TAGACCA','ATACA','GATTACA']
    return dict_seq, list_seq


def read_fasta_to_dic3(infasta):  # 适用于带详细位置cds的fa文件
    with open(infasta, 'r') as f:
        seq_id = ''  # 基因名
        dict_seq = {}  # 基因名-序列
        dict_len = {}  # 基因名-长度
        dict_pos = {}  # 基因名-位置,sort()
        d_pos = {}  # 基因名-位置,未排序
        for line in f:
            list_gene_pos = []  # 某基因对应的位置组成的列表 sort
            l_gene_pos = []
            list_n = [0]  # 计算同名基因是第几个
            if line.startswith('>'):
                seq_id = line.strip('\n').split()[2].split('=')[
                    1].strip(']')  # 基因名
                if seq_id in dict_seq.keys():
                    list_n.append(0)
                    seq_id = seq_id+'-'+str(len(list_n))  # 基因名+1 ycf1-2形式
                list_tmp = re.findall(
                    r'\d+', line.strip('\n').split()[1].lstrip(
                        '[').rstrip(']'))  # 位置打散成一个个起点或终点
                [l_gene_pos.append(int(i))
                 for i in list_tmp]  # 转换成数字,放进l_gene_pos,未排序
                d_pos[seq_id] = l_gene_pos

                l_gene_pos.sort()
                [list_gene_pos.append(i) for i in l_gene_pos]
                dict_pos[seq_id] = list_gene_pos
                dict_seq[seq_id] = ''
                dict_len[seq_id] = ''
            else:
                dict_seq[seq_id] += line.strip('\n')
                dict_len[seq_id] += str(len(line.strip('\n')))
    print('{0} Item Total: {1} {2} {3}'.format(os.path.basename(infasta),
                                               len(dict_seq), len(dict_len), len(dict_pos)))
    return dict_seq, dict_len, dict_pos, d_pos


if __name__ == '__main__':
    dict_seq = read_fasta_to_dic1(args.infile)
    dict_seq, list_seq = read_fasta_to_dic2(args.infile, args.outfile)
    dict_seq, dict_len, dict_pos, d_pos = read_fasta_to_dic3(args.infile)
