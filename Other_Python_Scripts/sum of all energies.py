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
sum_probabilities = sympy.summation(probability,(E, 1, 100)) #difference between 90 and 100keV is quite small, negligible.
depth_range = data[:,1]


probability_function = sympy.lambdify(z,sum_probabilities,"numpy")
py.plot(depth_range,probability_function(depth_range))
py.xlabel("depth (ucm)")
py.ylabel("P(z,E)")
py.title("P(z,E) as implantation depth varies \n(A=4*10^6 g ucm^-2 keV^-n, n=1.6, p=3.21*10^18 g ucm^-3 m=1.28", y=1.03)
py.show()