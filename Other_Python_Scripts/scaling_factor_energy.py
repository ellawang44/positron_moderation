from __future__ import division
import pylab as py
import sympy
import numpy
import scipy.stats as stats

data = py.loadtxt('random_energies.csv',delimiter=',')
if len(data.shape) > 1:
    depths = data[:,0]
else:
    depths = data
    
E = sympy.symbols('E')
c = 3*10**8 #m/s
m = 0.511 #MeV/(c^2)
Q = 1.4743 #MeV
p_energy = ((E**2+2*E*m)**(1/2))*((Q-E)**2)*(E+m) #C/c^5 is scaling factor roughly 1

scaling_range = xrange(1000,1500)
scaling_factor = 1332 #lowest chi squared value

def scaled_model(x):
    """Scale a model by a given factor"""
    return scaling_factor*p_energy

hist_values,_,_,_ = stats.histogram(depths,numbins=145)
    
def chisquared_calc(scaling_value):
    
    """Calculate the chi squared value when the model is multiplied
    by the scaling factor
    """
    
    total = dof = 0
    for (hn,fn) in zip(hist_values,(p_energy.subs(E,i) for i in numpy.arange(0.005,1.475,0.01))):
        if hn == 0: break
        dof +=1
        total += (hn - scaling_value * fn)**2 / hn
    return total/dof

if __name__ == '__main__':    
    py.plot(numpy.arange(0.005,1.475,0.01), [scaling_factor*p_energy.subs(E,i) for i in numpy.arange(0.005,1.475,0.01)])
    py.hist(depths,bins=145)
    py.xlabel('energy (MeV)')
    py.ylabel('number of occurances')
    py.title("scaled model on top of randomly generated energy values \n scaling factor =1332, chi squared = 1.9758")
    py.show()

if __name__ == '__othermain__':
    py.plot(scaling_range,[chisquared_calc(x) for x in scaling_range])
    py.xlabel('scaling factor')
    py.ylabel('chi squared')
    py.title("scaling factor and corresponding chi squared values")
    py.show()

    print min([chisquared_calc(x) for x in scaling_range])
    print min(scaling_range,key=chisquared_calc)