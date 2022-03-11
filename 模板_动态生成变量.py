createVar = locals()  # 核心是locals()这个内置函数
listTemp = range(1, 11)
for i, s in enumerate(listTemp):  # i是listtemp的角标,s是值,这个内置函数是特殊的迭代器
    createVar['a'+str(s)] = s

print(createVar['a9'], createVar['a10'])


def test_list_pre():
    prepare_list = locals()
    for i in range(16):
        prepare_list['list_' + str(i)] = []
        prepare_list['list_' + str(i)].append(('我是第' + str(i)) + '个list')
    print(prepare_list['list_0'])
    print(prepare_list['list_1'])
    print(prepare_list['list_2'])
    print(prepare_list['list_15'])


if __name__ == '__main__':
    test_list_pre()
