from __future__ import division
import sympy
import pylab
import numpy
import math
from matplotlib import pyplot

c = 3*10**8 #m/s
#m = 9.109e-37 #MeV/(c^2)
m = 0.511
Q = 1.4743 #MeV
unity = 1/0.748611
E = sympy.symbols ('E')
p_energy = unity*((E**2+2*E*m)**(1/2))*((Q-E)**2)*(E+m) #C/c^5 is unity
x = numpy.arange(0,1.473,0.01) 
x_log = numpy.arange(0.01,1.4743,0.01)

probability = [float(p_energy.subs(E,i)) for i in x]
#probability_log = [float(p_energy.subs(E,i)) for i in x_log]
#log_probability = numpy.log(numpy.array(probability_log))
#log_energy = numpy.log(x_log)
'''log_energy = []
for i in x:
    loge = [math.log(i)]
    log_energy += loge
'''    
#pylab.plot(x, probability)
pylab.axes(xscale='log',yscale='log')
pylab.plot(x, probability)
pylab.xlabel("energy(MeV)")
pylab.ylabel("N(E)")
pylab.title("energy probability function")
pylab.show()
