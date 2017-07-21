import sympy
#import numpy
import pylab as py
import random

E,z = sympy.symbols('E z')
A = 4*10**(6) #g/ucm^2keV^n
p = 3.21*10**18 #g/ucm^3 can assume pure SiC because radioactive source is only a small amount
n = 1.6
m = 1.28
median_depth = A*(E**n)*(10**6)/p

def rejection_sampling(probability_function,domain_sup,codomain_sup,substitute):
    x = random.uniform(0,domain_sup)
    y = random.uniform(0,codomain_sup)
    while probability_function.subs(substitute,x).evalf()<y:
        x = random.uniform(0,domain_sup)
        y = random.uniform(0,codomain_sup)
    return x
    
p_energy = 0.0002*E #0-100 for probability of 1
zprime = 2.58*median_depth
p_depth = lambda: (sympy.diff((-(sympy.exp(-((z*10**(-6)/zprime)*(1+z*10**(-6)/zprime)**2)**m))), z)).subs(E, rejection_sampling(p_energy,100,0.02,E))

list_of_depths = [float(rejection_sampling(p_depth(),50,0.2,z)) for _ in xrange(10)]

#py.hist(list_of_depths,bins=100)
#py.xlabel("energy")
#py.ylabel("number of occurances")
#py.show()

