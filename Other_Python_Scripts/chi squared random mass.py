from __future__ import division
import pylab as py
import sympy
import numpy
import scipy.stats as stats
import math

data = py.loadtxt('random_masses.csv',delimiter=',')
if len(data.shape) > 1:
    depths = data[:,0]
else:
    depths = data

m,x,r = sympy.symbols('m x r')
make_my_code_run_fast = 10**9
rho = 3210
min_mass = 4.6637066*10**(-26)+1.9944235*10**(-26)
max_mass = 10**(-13)
mass_to_radius = (3*x/(4*rho*math.pi))**(1/3)
min_radius = make_my_code_run_fast*mass_to_radius.subs(x,min_mass)
max_radius = make_my_code_run_fast*mass_to_radius.subs(x,max_mass)
radius_range = numpy.arange(0,max_radius,10**(0))
F = sympy.Piecewise((0,r<min_radius),(0,r>max_radius),((4*r**3*math.pi*rho/3)**(-3/2),True))
howdoyoutell = numpy.arange(1.7*10**(-10),3.35*10**(-9),0.5*10**(-11))

scaling_range = xrange(13000,15000)
scaling_factor = 10000 #lowest chi squared value

def scaled_model(x):
    """Scale a model by a given factor"""
    return scaling_factor*F

hist_values,_,_,_ = stats.histogram(depths,numbins=318)
    
def chisquared_calc(scaling_value):
    
    """Calculate the chi squared value when the model is multiplied
    by the scaling factor
    """
    
    total = dof = 0
    for (hn,fn) in zip(hist_values,(F.subs(m,i) for i in howdoyoutell)):
        if hn == 0: break
        dof +=1
        total += (hn - scaling_value * fn)**2 / hn
    return total/dof

if __name__ == '__othermain__':    
    py.plot(howdoyoutell, [scaling_factor*F.subs(m,i) for i in howdoyoutell])
    py.hist(depths,bins=318)
    py.xlabel('mass (kg)')
    py.ylabel('number of occurances')
    py.title("scaled model on top of randomly generated masses \n scaling factor =1332, chi squared = 1.9758")
    py.show()

if __name__ == '__main__':
    py.plot(scaling_range,[chisquared_calc(x) for x in scaling_range])
    py.xlabel('scaling factor')
    py.ylabel('chi squared')
    py.title("scaling factor and corresponding chi squared values")
    py.show()

    print min([chisquared_calc(x) for x in scaling_range])
    print min(scaling_range,key=chisquared_calc)