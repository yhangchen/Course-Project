from Subsets import subsets
class LTL:
    def __init__(self, data, left=None, right=None):
        self.data = data
        self.left = left
        self.right = right

    @staticmethod
    def make_tree(exp):#输入，分析详见报告。
        test = exp.strip()[1:-1]
        testl = test.split()
        if len(testl) == 1:# 检查LTL公式是否由一个单独的基本命题构成
            return LTL(testl[0], None, None)
        st = [] # 存储算符的栈
        exl = []# 存储表达式的表
        exp = exp.strip()
        if not exp:
        	return None
        if exp[-1] != ')':# 给不标准的公式末尾加上括号
        	exp = exp + ')'
        exp = exp.replace('weakuntil', 'weak') # 将weakuntil换为weak方便处理
        while exp:
            if exp[0] == '(':# 打开括号
                exp = exp[1:].strip()
            elif exp[0] == ')':# 处理表达式
                data = st.pop()
                r, l = exl.pop(), exl.pop()
                exl.append(LTL(data, l, r))
                exp = exp[1:].strip()
            elif exp[:3] == 'not':# 处理not
                exl.append(None)
                st.append('not')
                exp = exp[3:].strip()
                if exp[0] == '(':
                    continue
                else:
                    c = exp.find(')')
                    exl.append(LTL(exp[:c].strip(), None, None))
                    exp = exp[c:]
                    continue
            elif exp[:4] == 'next':# 处理next
                exl.append(None)
                st.append('next')
                exp = exp[4:].strip()
                if exp[0] == '(':
                    continue
                else:
                    c = exp.find(')')
                    exl.append(LTL(exp[:c].strip(), None, None))
                    exp = exp[c:]
                    continue
            elif exp[:6] == 'always':#处理always
                exl.append(None)
                st.append('always')
                exp = exp[6:].strip()
                if exp[0] == '(':
                    continue
                else:
                    c = exp.find(')')
                    exl.append(LTL(exp[:c].strip(), None, None))
                    exp = exp[c:]
                    continue
            elif exp[:10] == 'eventually':# 处理eventually
                exl.append(None)
                st.append('eventually')
                exp = exp[10:].strip()
                if exp[0] == '(':
                    continue
                else:
                    c = exp.find(')')
                    exl.append(LTL(exp[:c].strip(), None, None))
                    exp = exp[c:]
                    continue
            elif 'and' in exp or 'until' in exp or 'or' in exp or 'weak' in exp:
                a, b, d, w = exp.find('and'), exp.find('until'), exp.find('or'), exp.find('weak')
                if b != -1 and (a == -1 or b < a) and (d == -1 or b < d) and (w == -1 or b < w):    # until
                    if b == 0:
                        st.append('until')
                        exp = exp[b+5:].strip()
                        if exp[0] == '(':
                            continue
                        else:
                            c = exp.find(')')
                            exl.append(LTL(exp[:c].strip(), None, None))
                            exp = exp[c:]
                            continue
                    else:
                        exl.append(LTL(exp[: b].strip(), None, None))
                        st.append('until')
                        exp = exp[b+5:].strip()
                        if exp[0] == '(':
                            continue
                        else:
                            c = exp.find(')')
                            exl.append(LTL(exp[:c].strip(), None, None))
                            exp = exp[c:]
                            continue
                elif a != -1 and (b == -1 or a < b) and (d == -1 or a < d) and (w == -1 or a < w):     # and
                    if a == 0:
                        st.append('and')
                        exp = exp[a+3:].strip()
                        if exp[0] == '(':
                            continue
                        else:
                            c = exp.find(')')
                            exl.append(LTL(exp[:c].strip(), None, None))
                            exp = exp[c:]
                            continue
                    else:
                        exl.append(LTL(exp[: a].strip(), None, None))
                        st.append('and')
                        exp = exp[a+3:].strip()
                        if exp[0] == '(':
                            continue
                        else:
                            c = exp.find(')')
                            exl.append(LTL(exp[:c].strip(), None, None))
                            exp = exp[c:]
                            continue
                elif d != -1 and (a == -1 or d < a) and (b == -1 or d < b) and (w == -1 or d < w):   # or
                    if d == 0:
                        st.append('or')
                        exp = exp[d+2:].strip()
                        if exp[0] == '(':
                            continue
                        else:
                            c = exp.find(')')
                            exl.append(LTL(exp[:c].strip(), None, None))
                            exp = exp[c:]
                            continue
                    else:
                        exl.append(LTL(exp[: d].strip(), None, None))
                        st.append('or')
                        exp = exp[d+2:].strip()
                        if exp[0] == '(':
                            continue
                        else:
                            c = exp.find(')')
                            exl.append(LTL(exp[:c].strip(), None, None))
                            exp = exp[c:]
                            continue
                elif w != -1 and (a == -1 or w < a) and (b == -1 or w < b) and (d == -1 or w < d): # weakuntil
                    if w == 0:
                        st.append('weakuntil')
                        exp = exp[w+4:].strip()
                        if exp[0] == '(':
                            continue
                        else:
                            c = exp.find(')')
                            exl.append(LTL(exp[:c].strip(), None, None))
                            exp = exp[c:]
                            continue
                    else:
                        exl.append(LTL(exp[: w].strip(), None, None))
                        st.append('weakuntil')
                        exp = exp[w+4:].strip()
                        if exp[0] == '(':
                            continue
                        else:
                            c = exp.find(')')
                            exl.append(LTL(exp[:c].strip(), None, None))
                            exp = exp[c:]
                            continue
            else:
                raise ValueError('Wrong_input!')
        #以下为若干等价替换方法。
        def trans_weakuntil(ft):
            queue = [ft]
            while queue:
                t = queue.pop()
                if t.data == 'weakuntil':
                    l = t.left
                    r = t.right
                    lct = LTL('until', l, r)
                    rct = LTL('always', None, r)
                    t.data = 'or'
                    t.left = lct
                    t.right = rct
                if t.left is not None:
                    queue = [t.left] + queue
                if t.right is not None:
                    queue = [t.right] + queue
            return ft

        def trans_or(ft):
            queue = [ft]
            while queue:
                t = queue.pop()
                if t.data == 'or':
                    l = t.left
                    r = t.right
                    rclct = LTL('not', None, l)
                    rcrct = LTL('not', None, r)
                    rct = LTL('and', rclct, rcrct)
                    t.data = 'not'
                    t.left = None
                    t.right = rct
                if t.left is not None:
                    queue = [t.left] + queue
                if t.right is not None:
                    queue = [t.right] + queue
            return ft

        def trans_always(ft):
            queue = [ft]
            while queue:
                t = queue.pop()
                if t.data == 'always':
                    r = t.right
                    rcrct = LTL('not', None, r)
                    rct = LTL('eventually', None, rcrct)
                    t.data = 'not'
                    t.right = rct
                if t.left is not None:
                    queue = [t.left] + queue
                if t.right is not None:
                    queue = [t.right] + queue
            return ft

        def trans_eventually(ft):
            queue = [ft]
            while queue:
                t = queue.pop()
                if t.data == 'eventually':
                    t.data = 'until'
                    t.left = LTL('true', None, None)
                if t.left is not None:
                    queue = [t.left] + queue
                if t.right is not None:
                    queue = [t.right] + queue
            return ft

        formula_tree = exl[0]
        formula_tree = trans_weakuntil(formula_tree)
        formula_tree = trans_or(formula_tree)
        formula_tree = trans_always(formula_tree)
        formula_tree = trans_eventually(formula_tree)

        return formula_tree

    #以下是检验基本子集的三个静态方法。
    @staticmethod
    def consistent(subset, closure):
        if LTL('true') in closure and LTL('true') not in subset:
            return False
        for i in subset:
            if i.negative() in subset:
                return False
        for i in subset:
            if i.data == 'and':
                l = i.left
                r = i.right
                if not (l in subset and r in subset):
                    return False
        for i in subset:
            for j in subset:
                if i.together(j) in closure and i.together(j) not in subset:
                    return False
        return True 

    @staticmethod
    def locally_consistent(subset, closure):
        for i in closure:
            if i.data == 'until':
                i1 = i.left
                i2 = i.right
                if i2 in subset and i not in subset:
                    return False
                if i in subset and i2 not in subset and i1 not in subset:
                    return False
        return True

    @staticmethod
    def maximal(subset, closure):
        for i in closure:
            if i not in subset and i.negative() not in subset:
                return False
        return True

    def negative(self):#给出一个论断的否定。且仅仅把not not a与a等同。
        if self.data == 'not':
            return self.right
        else:
            return LTL('not', None, self)

    def together(self, another):#将两个论断以and连接。
        return LTL('and', self, another)

    def until(self, another):#将论断前加until。
        return LTL('until', self, another)

    def sub(self):#给出LTL公式的子公式。
        ans = []
        s = []
        t = self
        while t is not None or s:
            while t is not None:
                ans.append(t)
                s.append(t.right)
                t = t.left
            t = s.pop()
        return ans

    def closure(self):#给出LTL公式的闭包。
        subformula = set(self.sub())
        s = []
        for i in subformula:
            s.append(i.negative())
        s = set(s)
        return s.union(subformula)

    def elem_sub(self):#给出LTL公式的基本子集。注意调用了三个方法。
        closure = self.closure()
        all_subsets = subsets(closure)
        elem_sub = []
        for subset in all_subsets:
            if LTL.consistent(subset,closure) and LTL.locally_consistent(subset,closure) and LTL.maximal(subset,closure):
                elem_sub.append(subset)
        return elem_sub

    def __eq__(self,another):
        return hash(self) == hash(another)

    def pre(self):#LTL公式的前缀表示。
        ans = []
        s = []
        t = self
        while t is not None or s:
            while t is not None:
                ans.append(t.data)
                s.append(t.right)
                t = t.left
            t = s.pop()
        return ans

    def middle(self):#LTL公式的中缀表示
        ans = []
        s = []
        t = self
        while t is not None or s:
            while t is not None:
                s.append(t)
                t = t.left
            t = s.pop()
            ans.append(t.data)
            t = t.right
        return ans

    def __hash__(self): #利用先根和中根可以完全确定二叉树。
        pre = self.pre()
        middle = self.middle()
        pre.extend(middle)
        return hash(tuple(pre))

    def store(self, exp): #为记录exp的原始形式。便于最后输出。
        self.content = exp

    def prt(self):            # debug用，打印树,先根序
        print(self.pre())




if __name__ == '__main__': 
    ltl = LTL.make_tree('(true)')
    for i in ltl.elem_sub():
        for j in i:
            j.prt()
