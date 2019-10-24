from Subsets import subsets
from LTL_Formula import LTL
from TransitionSystem import TransitionSystem
from Data_transform import LTL2NBA, product
def persistence_checking(TS, F): #检验问题。
    R = set(); U = []; T = set(); V = []; cycle_found = False #初始化没有圈。
    def cycle_check(s):#找圈。内部深搜。
        nonlocal cycle_found, T, V
        T = set(); V = []#V用来记录路径，T用来记录经过的状态，子函数内部初始化为空。
        cycle_found = False
        V.append(s)
        T.add(s)
        while True:
            if not V or cycle_found:
                break
            s_ = V[-1]
            if s in TS.post(s_):#存在圈
                cycle_found = True
                V.append(s)
            else:
                new = TS.post(s_) - T
                if new:
                    s2 = list(new)[0]
                    V.append(s2)
                    T.add(s2)
                else:
                    V.pop()
        return cycle_found


    def reachable_cycle(s):#找路径。外部深搜。
        nonlocal cycle_found
        #U用来记录路径，R用来记录经过的状态。
        U.append(s)
        R.add(s)
        while True:
            if not U or cycle_found:
                break
            s_ = U[-1]
            new = TS.post(s_) - R
            if new:
                s2 = list(new)[0]
                U.append(s2)
                R.add(s2)
            else:
                U.pop()
                if s_[1] in F:
                    cycle_found = cycle_check(s_)

    while set(TS.I()) - R and not cycle_found:#从还没有被考虑过的初始状态开始找带圈路径。（详细定义见报告）
        new = set(TS.I()) - R
        s = list(new)[0]
        reachable_cycle(s)

    if not cycle_found:
        return True
    else:
        U.extend(V)#边和圈连起来。
        U = [i[0] for i in U]#投射到原来TS上。
        return (False, U)