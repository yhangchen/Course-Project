def subsets(s):#求一个集合的子集，返回以元组的列表形式给出。
    n = len(s)
    s = list(s)
    subsets0 = []
    for i in range(2**n):
        subset = []
        for j in range(n):
            if (i>>j)%2 == 1:
                subset.append(s[j])
        subsets0.append(tuple(subset))
    return subsets0
