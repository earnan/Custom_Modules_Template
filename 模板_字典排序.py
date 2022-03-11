def readfasta(input_file, output_file):
    fa = open(input_file, 'r')
    fo = open(output_file, 'w')
    d = {}
    for line in fa:
        id = int(line.split()[1].split('-')[0])
        d[id] = line  # 有key无value
    print(d)  # 包含说明行序列行的字典结构
    return d


d = {151: 'atp 151-188:+', 10: 'ef1 10-50:+', 51: 'alt 51-61:-;75-148:-'}
# 以下为字典按键排序的两种方法
for i in sorted(d):
    print(i, d[i], end="        ")

# x[0] 改为 x[1]则是按值排序
d_order = sorted(d.items(), key=lambda x: x[0], reverse=False)
print('\n')
print(d_order)
