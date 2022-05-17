#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   dynamic_generating_variable.py
#         Author:   yujie
#    Description:   dynamic_generating_variable.py
#        Version:   1.0
#           Time:   2022/03/14 16:09:44
#  Last Modified:   2022/03/14 16:09:44
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
def create_var(inlist):
    createvar = locals()  # 核心是locals()这个内置函数
    listTemp = inlist
    for i, s in enumerate(listTemp):  # i是listtemp的角标,s是值,这个内置函数是特殊的迭代器
        createvar['list_'+str(s)] = []
    #print(createvar['list_9'], createvar['list_10'])
    return [createvar['list_'+str(s)] for s in listTemp]


if __name__ == '__main__':
    inlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    for i in create_var(inlist):
        print(i)

"""
def test_list_pre():
    createVar = locals()  # 核心是locals()这个内置函数
    listTemp = range(1, 11)
    for i, s in enumerate(listTemp):  # i是listtemp的角标,s是值,这个内置函数是特殊的迭代器
        createVar['a'+str(s)] = s

    print(createVar['a9'], createVar['a10'])
    prepare_list = locals()
    for i in range(16):
        prepare_list['list_' + str(i)] = []
        prepare_list['list_' + str(i)].append(('我是第' + str(i)) + '个list')
    print(prepare_list['list_0'])
    print(prepare_list['list_1'])
    print(prepare_list['list_2'])
    print(prepare_list['list_15'])
"""
