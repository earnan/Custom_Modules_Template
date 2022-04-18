def judgment_input_type(string):
    if os.path.isfile(string):
        abs_path = os.path.abspath(string)
        abs_dir = os.path.dirname(abs_path)
        input_file = open(string, 'r')
        for line in input_file:
            if len(line.strip('\n')):  # and (not line.startswith('>')):
                seq = line.strip('\n')
    else:  # type(string) == type("a"):
        seq = string
        abs_dir = os.getcwd()
    return seq, abs_dir


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
