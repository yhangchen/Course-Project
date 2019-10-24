from LTL_Formula import LTL
from Subsets import subsets
from TransitionSystem import TransitionSystem
def LTL2GNBA(ltl, AP):#LTL转换为GNBA。意义见报告。
    Q = ltl.elem_sub()
    cl = ltl.closure()
    Q0 = []
    AP = set(AP)
    for i in Q:
        if ltl in i and i not in Q0:
            Q0.append(i)
    F = []
    for i in ltl.closure():
        if i.data == 'until':
            f = []
            i1 = i.left
            i2 = i.right
            for j in Q:
                if i not in j or i2 in j:
                    f.append(j)
            F.append(f)
    sigma = subsets(AP)
    def delta(B, A):
        if set(A) & cl != set(B) & AP:
            return []
        else:
            ans = []
            for B_ in Q:
                for i in ltl.closure():
                    if i.data == 'next' and B_ is not None:
                        try:
                            assert (i in B) == (i.right in B_)
                        except:
                            B_ = None
                for j in ltl.closure():
                    if j.data == 'until' and B_ is not None:
                        j1 = j.left
                        j2 = j.right
                        try:
                            assert (j in B) == (j2 in B or (j1 in B and j in B_))
                        except:
                            B_ = None
                if B_ is not None and B_ not in ans:
                    ans.append(B_)
            return ans
    if not F:
        F = [Q]
    return (Q, sigma, delta, Q0, F)


def GNBA2NBA(GNBA):#GNBA转换为NBA。意义见报告。
    (Q, sigma, delta, Q0, F) = GNBA
    k = len(F)
    Q_ = [(q,i) for q in Q for i in range(k)]
    Q0_ = [(q,0) for q in Q0]
    F_ = [(f,0) for f in F[0]]
    def delta_(x, A):
        (q, i) = x
        if q not in F[i]:
            return [(q_, i) for q_ in delta(q, A)]
        else:
            return [(q_, i+1 % k) for q_ in delta(q, A)]
    return (Q_, sigma, delta_, Q0_, F_)

def LTL2NBA(ltl, AP):
    return GNBA2NBA(LTL2GNBA(ltl, AP))

def product(TS, NBA):#构造乘积TS。意义见报告。
    (S, Act, Arrow, I, AP, L) = TS.items()
    (Q, sigma, delta, Q0, F) = NBA
    S_ = [(s,q) for s in S for q in Q]
    Arrow_ = []
    for (s, alpha, t) in Arrow:
        for q in Q:
            for p in delta(q, [L(t)]):
                Arrow_.append(((s,q),alpha,(t,p)))
    I_ = []
    for s0 in I:
        for q0 in Q0:
            for q in delta(q0, [L(s0)]):
                I_.append((s0, q))
                
    AP_ = [LTL(i) for i in Q]#形式化定义。以后无用。
    def L_(x):#形式化定义。以后无用。
        (s, q) = x
        return set(q)

    new_TS = TransitionSystem(S_, Act, Arrow_, I_, AP_, L_)
    return (new_TS, F)

