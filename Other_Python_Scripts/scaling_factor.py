from __future__ import division
import pylab as py
import sympy
import numpy
import scipy.stats as stats
from slice_sampling import (z, E, probability)

data = py.loadtxt('random_depths.csv',delimiter=',')
if len(data.shape) > 1:
    depths = data[:,0]
else:
    depths = data
    
sum_probabilities = sympy.lambdify(z,sympy.summation(probability,(E, 1, 100)),'numpy') 

scaling_range = xrange(8000,10500)
scaling_factor = 9527 #lowest chi squared value

def scaled_model(x):
    """Scale a model by a given factor"""
    return scaling_factor*sum_probabilities(x)

hist_values,_,_,_ = stats.histogram(depths,numbins=502)
    
def chisquared_calc(scaling_value):
    
    """Calculate the chi squared value when the model is multiplied
    by the scaling factor
    """
    
    total = dof = 0
    for (hn,fn) in zip(hist_values,(sum_probabilities(numpy.arange(5,5020,10)))):
        if hn == 0: break
        dof +=1
        total += (hn - scaling_value * fn)**2 / hn
    return total/dof

if __name__ == '__main__':    
    py.plot(numpy.arange(1,5020), scaled_model(numpy.arange(1,5020)))
    py.hist(depths,bins=502)
    py.xlabel('depth (ucm)')
    py.ylabel('number of occurances')
    py.title("scaled model on top of randomly generated depth values")
    py.show()

if __name__ == '__othermain__':
    py.plot(scaling_range,[chisquared_calc(x) for x in scaling_range])
    py.xlabel('scaling factor')
    py.ylabel('chi squared')
    py.title("scaling factor and corresponding chi squared values")
    py.show()

    print min([chisquared_calc(x) for x in scaling_range])
    print min(scaling_range,key=chisquared_calc)