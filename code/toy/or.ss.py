# here we analytically solve an expression for 2x2 OR at steady state
# i.e. among "very old" people / ignoring exposure/outcome status among model entrants

from sympy import symbols,Eq,solve,simplify

# notation: (exp = exposure, out = outcome)
# {variables}: y00: exp=0,out=0 | y10: exp=1,out=0 | y01: exp=0,out=1 | y11: exp=1,out=1
# {base rates}: a: exp.onset | b: exp.recov | c: out.onset | b: out.recov
# {IRRs}: K: IRR out.onset while exp | R: IRR out.recov while exp
y00,y10,y01,y11 = symbols('y00 y10 y01 y11',positive=True)
a,b,c,d,K,R,OR = symbols('a b c d K R OR',positive=True)

eqs = [ # system of equations
  Eq(1, y00+y10+y01+y11),                                 # sum = 1
  Eq(0, -y00*(a+c) +y10* b      +y01*   d              ), # d/dt(y00) = 0
  Eq(0, +y00* a    -y10*(b+c*K)            +y11*   d*R ), # d/dt(y10) = 0
  Eq(0, +y00*   c               -y01*(a+d) +y11* b     ), # d/dt(y01) = 0
  Eq(0,            +y10*   c*K  +y01* a    -y11*(b+d*R)), # d/dt(y11) = 0
  # Eq(OR, (y00*y11)/(y01+y10)), # slows solving for some reason
]
X = solve(eqs,[y00,y10,y01,y11]) # get analytic expressions for y00,y10,y01,y11
# for y in X: print(y,':',simplify(X[y])) # yikes!
print('OR :',simplify( (X[y00]*X[y11]) / (X[y10]*X[y01]) )) # luckily, many terms cancel
