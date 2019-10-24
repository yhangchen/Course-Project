from LTL_Formula import LTL
from Subsets import subsets
class TransitionSystem:         
    def __init__(self, S, Act, Arrow, I, AP, L):  
        self._S = tuple(S)
        self._Act = tuple(Act)
        self._Arrow = tuple(Arrow)
        self._I = tuple(I)
        self._AP = tuple(AP)
        self._L = L

    def S(self):
        return self._S

    def Act(self):
        return self._Act

    def Arrow(self):
        return self._Arrow

    def I(self):
        return self._I

    def AP(self):
        return self._AP

    def L(self):
        return self._L

    def items(self):
        return (self._S, self._Act, self._Arrow, self._I, self._AP, self._L)

    @staticmethod
    def input(states_list, transitions_list, initial_states):
        S = tuple(states_list)
        Act = [t[1] for t in transitions_list]
        Arrow = transitions_list
        I = tuple(initial_states)
        AP = tuple(LTL(i) for i in states_list) #这里假设AP就为S。即LTL是直接关于S的状态描述。
        L = lambda x: LTL(x)#形式化定义，此后并无用处。
        return TransitionSystem(S, Act, Arrow, I, AP, L)


    def post(self, s):#状态s的下一个状态。
        ans = []
        for i in self.Arrow():
            if i[0] == s:
                ans.append(i[2])
        return set(ans)

    def complete(self):#若有终状态，转换为无终状态的TS。
        ter = []
        S = list(self.S())
        Act = list(self.Act())
        Arrow = list(self.Arrow())
        for i in S:
            if not self.post(i):
                ter.append(i)
        if not ter:
            return self
        else:
            S.append(None)
            Arrow.append((None, ' ', None))#因为边不可能为空格，状态不会是None，所以这里不会造成重复。
            for i in ter:
                Arrow.append((i, ' ', None))
            return TransitionSystem(S, Act, Arrow, self.I(), self.AP(), self.L())


if __name__ == '__main__':#试验。
    from input import *
    filename = 'LTL.txt'
    input_data = read_TS(filename)
    TS = TransitionSystem.input(input_data[0], input_data[1], input_data[2])
    for i in TS.AP():
        i.prt()
    print(TS.items())