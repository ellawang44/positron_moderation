import pylab as py
import sympy
import numpy

numpy.seterr(over='raise')
data = py.loadtxt("ASC_python.csv", delimiter=",")

A = 4*10**(6) #g/ucm^2keV^n
p = 3.21*10**18 #g/ucm^3 can assume pure SiC because radioactive source is only a small amount
n = 1.6
m = 1.28
z, E = sympy.symbols('z E')
median_depth = A*(E**n)*(10**6)/p

zprime = 2.58*median_depth
probability = (sympy.diff((-(sympy.exp(-((z*10**(-6)/zprime)*(1+z*10**(-6)/zprime)**2)**m))), z))
depth_range = data[:,1]

probability_function = sympy.lambdify((E, z),probability,"numpy")

z_range = range(1,41)
e_range = range(1,31)

import csv
with open('lookup_table.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['','implantation depth (ucm)'] + z_range)
    writer.writerow(['energy (keV)'])
    for E in e_range:
        writer.writerow([E,''] + [probability_function(E,z) for z in z_range])


#py.plot(depth_range,probability_function(depth_range))
#py.ylabel("P(z,E)")
#py.xlabel("depth (ucm)")
#py.title("probability of implanting a certain depth for 1keV")
#py.show()