from sympy import Symbol, S, symbols
from sympy.printing import pretty
import numpy as np

S1  = symbols('phi:3',   commutative=False)
S1T = symbols('P:3',     commutative=False)
P0  = symbols('theta:2', commutative=False)
P0T = symbols('T:2',     commutative=False)

X = S1T; Y = S1
reps = [(X[i]*Y[j], int(i==j)) for i in range(3) for j in range(3)]
X = P0T; Y = P0
reps += [(X[i]*Y[j], int(i==j)) for i in range(2) for j in range(2)]
def eval(expr): return  expr.expand().subs(reps)

for v in S1T:
    for u in S1:
        print('<%s, %s> = %s'%(pretty(v), pretty(u), eval(v*u)))
print()
for v in P0T:
    for u in P0:
        print('<%s, %s> = %s'%(pretty(v), pretty(u), eval(v*u)))
print()

def transpose(expr):
    expr = expr.expand()
    expr = expr.subs(zip(S1, symbols('x:3')))
    expr = expr.subs(zip(S1T, S1))
    expr = expr.subs(zip(symbols('x:3'), S1T))

    expr = expr.subs(zip(P0, symbols('x:2')))
    expr = expr.subs(zip(P0T, P0))
    expr = expr.subs(zip(symbols('x:2'), P0T))
    return expr

print(pretty(transpose(S1[0]))+'.T == '+pretty(transpose(transpose(S1[0]))))
print(pretty(transpose(P0[0]))+'.T == '+pretty(transpose(transpose(P0[0]))))
print()

def transposeOperator(d):
    T = transpose
    ants = d.__annotations__.copy()
    codom = ants.pop('return')
    dom = ants.popitem()[1]
    return lambda s: lambda p: d(p)(s)

def g(u : S1) -> P0:
    im = [P0[0], P0[1], -P0[0]-P0[1]] # [grad S1[0], grad S1[1], grad S1[2]]
    return eval(sum(eval(S1T[i]*u)*im[i] for i in range(3)))
for u in S1:
    print("g(%s) = %s"%(pretty(u), pretty(g(u))))
print()

def d(u : P0) -> S1T:
    T = transpose
    return lambda v: eval(T(u)*g(v))

for v in S1:
    for u in S1:
        print('<d(g(%s), %s> = %s'%(pretty(v), pretty(u), eval(d(g(v))(u))))
print()

T = transpose
dT = transposeOperator(d)

for v in S1T:
    for u in S1:
        print('<dT(%s), g(%s)> = %s'%(pretty(v), pretty(u), dT(T(v))(g(u))))
print()

def localMatrix(dT,g):
    ret = np.empty((3,3))
    for i in range(3):
        for j in range(3):
            ret[i,j] = eval(dT(T(S1T[i]))(g(S1[j])))
    return ret

print(localMatrix(dT, g))
