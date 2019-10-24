def read_TS(filename):     # 读取文件
    inf = open(filename, 'r')
    states_line = inf.readline()
    states_list = states_line.split()
    initial_states = inf.readline().split()
    transitions_list = []
    formula_list = []
    while True:
        line = inf.readline()
        l = line.split()
        if not l:
            break
        transitions_list.append((l[1], l[0], l[2]))
    while True:
        line = inf.readline()
        if not line:
            inf.close()
            break
        formula_list.append(line)
    return states_list, transitions_list, initial_states, formula_list

if __name__ == '__main__':
    print(read_TS('LTL.txt'))           # debug用
    a, b, c, d = read_TS('LTL.txt')     # debug用