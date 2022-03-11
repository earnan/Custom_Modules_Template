"""把fasta文件存储成俩列表，说明行 序列行各为一个列表"""


f = open('rosalind_gc.txt', 'r')

lines = f.readlines()

f.close()


def readfasta(lines):
    """读入并处理FASTA文件的函数"""

    seq = []

    index = []

    seqplast = ""

    numlines = 0

    for i in lines:

        if '>' in i:  # 说明行

            index.append(i.replace("\n", "").replace(">", ""))

            seq.append(seqplast.replace("\n", ""))  # 将已有半边序列存入seq列表末尾

            seqplast = ""  # 序列片段又重置为空

            numlines += 1

        else:  # 序列行

            seqplast = seqplast + i.replace("\n", "")  # 把分行的序列拼接成一个字符串

            numlines += 1

        if numlines == len(lines):

            seq.append(seqplast.replace("\n", ""))

    seq = seq[1:]
    print(index)
    print(seq)
    return index, seq


def countGC(list):

    numGC = list.count('G') + list.count('C')

    perGC = float(numGC) / len(list)

    return perGC * 100


(index, seq) = readfasta(lines)  # 接收序列名和序列
result = []

for i in seq:

    result.append(countGC(i))

maxID = index[result.index(max(result))].replace('>', "").replace("\n", "")

seqGC = max(result)

print(maxID)

print(round(seqGC, 6))

with open("05fasta_out.txt", 'w') as f2:
    f2.write(str(index)+str(seq))
