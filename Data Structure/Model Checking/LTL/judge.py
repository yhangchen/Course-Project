from LTL_Formula import LTL
from TransitionSystem import TransitionSystem
from graph import persistence_checking
from Data_transform import LTL2NBA, product
from input import read_TS
from Subsets import subsets

#初始化
filename = 'LTL.txt'
input_data = read_TS(filename)
TS = TransitionSystem.input(input_data[0], input_data[1], input_data[2])
TS = TS.complete()#假设无终状态。
ltls = []#LTL公式列表。
for i in input_data[3]:
    ltl = LTL.make_tree(i)
    if ltl is None:#空行。
        continue
    ltl.store(i)#储存原始形式。
    ltls.append(ltl)

def model_checking_ltl(TS, ltl):#对一个LTL公式检验。
    if ltl == LTL('true'):
        return True
    AP = TS.AP()
    NBA = LTL2NBA(ltl.negative(), AP)#生成NBA。
    new_TS = product(TS, NBA)[0]#构造乘积TS。
    F = product(TS, NBA)[1]#构造乘积TS上的eventually always语句。
    return persistence_checking(new_TS, F)#利用两次深搜判断。

def model_checking_ltls():#合理输出。
    for ltl in ltls:
        print(model_checking_ltl(TS, ltl), end=' ')
        print(ltl.content)


if __name__ == '__main__':
    model_checking_ltls()



