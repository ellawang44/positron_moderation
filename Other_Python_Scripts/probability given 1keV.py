from __future__ import division
import pylab as py
import sympy
import numpy

numpy.seterr(over='raise')
depths = numpy.arange(0,700,0.01)

A = 4*10**(6) #g/ucm^2keV^n
p = 3.21*10**18 #g/ucm^3 can assume pure SiC because radioactive source is only a small amount
n = 1.6
m = 1.28
z, E = sympy.symbols('z E')
median_depth = A*(E**n)*(10**6)/p

zprime = 2.58*median_depth
probability = (sympy.diff((-(sympy.exp(-((z*10**(-6)/zprime)*(1+z*10**(-6)/zprime)**2)**m))), z))

py.plot(depths,map(sympy.lambdify(z,probability.subs(E,10),"numpy"),depths),'b')
py.plot(depths,map(sympy.lambdify(z,probability.subs(E,20),"numpy"),depths),'g')
py.plot(depths,map(sympy.lambdify(z,probability.subs(E,30),"numpy"),depths),'red')
blue_patch = py.matplotlib.patches.Patch(color='blue', label='E=10keV')
green_patch = py.matplotlib.patches.Patch(color='green', label='E=20keV')
red_patch = py.matplotlib.patches.Patch(color='red', label='E=30keV')
py.matplotlib.pyplot.legend(handles=[blue_patch,green_patch,red_patch])
py.ylabel("P(z,E)")
py.xlabel("depth (ucm)")
py.title("probability of thermalising at a given depth \n for positrons of different energies")
py.show()