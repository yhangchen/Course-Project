# 定义全体变量operators表示所有的运算符（包括括号）
operators = {'(', ')','true', 'or', 'and', 'not', 'next', 'until', 'always', 'weakuntil', 'eventually', 'forall', 'exist'}

#定义栈
class StackUnderflow(ValueError):
    pass
class SStack():
    def __init__(self):
        self._elems = []
    def is_empty(self):
        return self._elems == []
    def top(self):
        if self._elems == []:
            raise StackUnderflow("in SStack.top()")
        return self._elems[-1]
    def push(self, elem):
        self._elems.append(elem)
    def pop(self):
        if self._elems == []:
            raise StackUnderflow("in SStack.pop()")
        return self._elems.pop()

# 读取文件
def read_TS(filename):     
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

#定义TS类
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

    # 定义构造一个TS的静态方法（输入所有状态，迁移列表和初始状态）
    @staticmethod
    def input(states_list, transitions_list, initial_states):
        S = tuple(states_list)
        Act = [t[1] for t in transitions_list]
        Arrow = transitions_list
        I = tuple(initial_states)
        AP = tuple(i for i in states_list) #这里假设AP就为S。即是直接关于S的状态描述
        L = lambda x: x 
        return TransitionSystem(S, Act, Arrow, I, AP, L)

    # 状态s的后继状态
    def post(self, s):
        ans = []
        for i in self.Arrow():
            if i[0] == s:
                ans.append(i[2])
        return set(ans)
    # 状态s的前驱状态
    def pre(self, s):
        ans = []
        for i in self.Arrow():
            if i[2] == s:
                ans.append(i[0])
        return set(ans)

# 构造CTL类（将公式存储在一棵二叉树中）
class CTL:
    def __init__(self, expression, TS):
        self._formula_tree = self.construct_tree(self.infix_order(TS.S(), operators, expression), self.trans_infix_suffix(TS.S(), operators, expression))
        self.sign = self._formula_tree.sign
        self.left = self._formula_tree.left
        self.right = self._formula_tree.right

    # 生成器函数，读取line中所有的运算符和运算元并逐一输出
    @staticmethod   
    def tokens(Propositions, Operators, line):
        i, llen = 0, len(line)
        while i < llen:
            while line[i].isspace():
                i += 1
                if i >= llen:
                    break
            judge = 0 # 判断符judge，判断是否在内层循环中yield一个值
            for j in range(i+1,llen+1):
                if line[i: j] in Operators or line[i: j] in Propositions:
                    yield line[i: j]
                    i = j
                    judge = 1
                    break
            if judge == 0:
                i += 1
                
    # line中运算符和运算元入树后的中根序列               
    @staticmethod
    def infix_order(Propositions, Operators, line):
        exp = []
        for x in CTL.tokens(Propositions, Operators, line):
            if x != '(' and x != ')':
                exp.append(x)
        return exp

    # 从中根序列到后根序列的转化    
    @staticmethod     
    def trans_infix_suffix(Propositions, Operators, line):
        st = SStack()
        exp = []

        for x in CTL.tokens(Propositions, Operators, line):
            if x in Propositions or x == 'true':
                exp.append(x)
            elif x == ')':
                while not st.is_empty() and st.top() != '(':
                    exp.append(st.pop())
                if st.is_empty():
                    raise SyntaxError("Missing '('.")
                st.pop()
            else:
                st.push(x)

        while not st.is_empty():
            if st.top() == '(':
                raise SyntaxError("Extra '('.")
            exp.append(st.pop())

        return exp

    # 由中根序列和后根序列生成二叉树    
    @staticmethod              
    def construct_tree(infix_order, suffix_order):
        if len(infix_order) == 0:
            return None
        root_data = suffix_order[-1]
        for i in range(len(infix_order)):
            if infix_order[i] == root_data:
                break
        left = CTL.construct_tree(infix_order[:i], suffix_order[:i])
        right = CTL.construct_tree(infix_order[i+1:], suffix_order[i: -1])
        return Node(root_data, left, right)

# 结点类（用于二叉树的构造与使用）               
class Node:                          
    def __init__(self, d, l, r):
        self.sign = d
        self.left = l
        self.right = r

# 后根序遍历二叉树           
def suf_order_tran(tree, proc):
    if tree is None:
        return
    suf_order_tran(tree.left, proc)
    suf_order_tran(tree.right,proc)
    proc(tree)

# 转换操作函数   
def proc(tree):
    if tree.sign == 'forall':
        if tree.right.sign == 'next':
            tree = Node('not', None, Node('exist', None, Node('next', None, Node('not', None, tree.right.right))))
        if tree.right.sign == 'and':
            tree = Node('and', Node('not', None, Node('exist', None, Node('until', Node('not', None, tree.right.right), Node('and', Node('not', None, tree.right.left), Node('not', None, tree.right.right))))), Node('not', None, Node('exist', None, Node('always', None, Node('not', None, tree.right.right)))))
        if tree.right.sign == 'eventually':
            tree = Node('not', None, Node('exist', None, Node('always', None, Node('not', None, tree.right.right))))
        if tree.right.sign == 'always':
            tree = Node('not', None, Node('exist', None, Node('until', true, Node('not', None, tree.right.right))))
    if tree.sign == 'weakuntil':
        tree = Node('not', None, Node('until', Node('and', tree.left, Node('not', None, tree.right)), Node('and', Node('not', None, tree.left), Node('not', None, tree.right))))
    if tree.sign == 'or':
        tree = Node('not', None, Node('and', Node('not', None, tree.left), Node('not', None, tree.right)))
    if tree.sign == 'eventually':
        tree = Node('until', Node('true', None, None), tree.right)

# 求交集函数
def intersection(A, B):
    C = {}
    for x in A:
        if x in B:
            C.add(x)
    return C

# 求补集函数
def complementary_set(A, S):
    S1 = {}
    for x in S:
        if x not in A:
            S1.add(x)
    return S1

# 应答函数
def answer(TS, a):
    ans = set()
    for s in TS.S():
        if a in TS.L()(s):
            ans.add(s)
    return ans

# 应答函数1   
def answer1(TS, tree):
    ans = set()
    for s in TS.S():
        if not intersection(TS.post(s), Sat(TS, tree)):
            ans.add(s)
    return ans

# 应答函数2
def answer2(TS, tree1, tree2):
    E = Sat(TS, tree2)
    T = E.copy()
    while E != set():
        s1 = E.pop()
        for s in TS.pre(s1):
            if s in Sat(TS, tree1) and s not in T:
                E.add(s)
                T.add(s)
    return T

# 应答函数3
def answer3(TS, tree):
    E = complementary_set(Sat(TS, tree), TS.S())
    T = Sat(TS, tree)
    count = {}
    for s in Sat(TS, tree):
        count[s] = len(TS.post(s))
    while E != set():
        s1 = E.pop()
        for s in TS.pre(s1):
            if s in T:
                count[s] -= 1
                if count[s] == 0:
                    T.remove(s)
                    E.add(s)
    return T

# Sat()计算函数   
def Sat(TS, tree):
    if tree.sign == 'true':
        return set(TS.S())
    if tree.left == None and tree.right == None:
        return answer(TS, tree.sign) 
    if tree.sign == 'and':
        return intersection(Sat(TS, tree.left),Sat(TS, tree.right))
    if tree.sign == 'not':
        return complementary_set(Sat(TS, tree.right), TS.S())
    if tree.sign == 'exist' and tree.right.sign == 'next':
        return answer1(TS, tree.right) 
    if tree.sign == 'exist' and tree.right.sign == 'until':
        return answer2(TS, tree.right.left, tree.right.right)
    if tree.sign == 'exist' and tree.right.sign == 'always':
        return answer3(TS, tree.right.right)

# 输出最终结果的函数
def FinalAnswer(TS, tree):
    if not Sat(TS, tree):
        return False
    for x in TS.I():
        if x not in Sat(TS, tree):
            return False
    return True

# 检验部分              
filename = 'CTL.txt'
input_data = read_TS(filename)
TS = TransitionSystem.input(input_data[0], input_data[1], input_data[2])
ctls = []
for ctl in input_data[3]:
    ctls.append(ctl)
check_answer = []
for ctl in ctls:
    temp = CTL(ctl, TS)
    suf_order_tran(temp, proc)
    check_answer.append(FinalAnswer(TS, temp))
# check_answer的每一位对应相应位置公式的的检验结果（True or False）
print(check_answer)  
