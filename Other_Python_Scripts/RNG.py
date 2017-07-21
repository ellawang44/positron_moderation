import pylab as py
import sympy
import random

E,z = sympy.symbols('E z')
A = 4*10**(6) #g/ucm^2keV^n
p = 3.21*10**18 #g/ucm^3 can assume pure SiC because radioactive source is only a small amount
n = 1.6
m = 1.28
median_depth = A*(E**n)*(10**6)/p

#energy function
p_energy = 0.0002*E
integrate_energy = sympy.integrate(p_energy,E)

#RNG
#r_energy = lambda: sympy.nsolve(integrate_energy - random.betavariate(0.1,1), E, 0)
#(a,b) a = lowest value, b = highest value, screw this up and your probability isn't 1 anymore
r_energy = lambda: sympy.nsolve(integrate_energy - random.uniform(0,1),E,(0,100), solver='bisect')

zprime = 2.58*median_depth
p_depth = (sympy.diff((-(sympy.exp(-((z*10**(-6)/zprime)*(1+z*10**(-6)/zprime)**2)**m))), z)).subs(E, r_energy())

sum_energy = sympy.sum(p_depth, (z,1,150))

#loops random energy generator
#list_of_energies = [r_energy for _ in range(100)]
#list_of_energies = [float(r_energy()) for _ in range(1000)]

#py.hist(list_of_energies,bins=100)
#py.xlabel("energy")
#py.ylabel("number of occurances")
#py.show()

print sum_energy