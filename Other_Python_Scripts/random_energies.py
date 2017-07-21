""" Slice sampling for approximating the penetration depth of
positrons released from Sc-44
"""

from __future__ import division
import pylab as py
import sympy
import random
import csv
import math
import scipy.optimize as opt

def slice_sampling(probability_function,substitute,initial_value=None,
                   width=100,burn=300,skip=0):
    
    """Returns a generator which generates values from a given probability
    distribution (which is defined using sympy). Like most Monte-Carlo
    Markov Chain methods, this requires a burn-in period to approach the correct
    distribution. This will be achieved faster if the initial_value supplied is
    close the the maximum of the function (if none is given, this will attempt
    to find it).
    
    Also due to being a Monte-Carlo Markov Chain method, this results 
    given are not statistically independent. If this is required increase
    the skip value as required. 
    """

    # Generate real python function
    p = sympy.lambdify(substitute,probability_function,"mpmath")
    
    # Calculate a good initial_value if none is given
    if initial_value is None:
        min_val = opt.minimize_scalar(lambda a:(-p(a)))
        initial_value = float(min_val.x)
    
    skip_count = 0
    while True:
        
        v_init = p(initial_value)
        v_slice = random.uniform(0,v_init)
        
        # Stepping out
        step_left = initial_value - width
        step_right = initial_value + width
        while p(step_left) > v_slice:
            step_left = step_left-width
        while p(step_right) > v_slice:
            step_right = step_right+width
        
        while True:
            assert(step_left<step_right)
            r_step = random.uniform(step_left,step_right)
            vr_step = p(r_step)
            if vr_step > v_slice: break
            
            # Shrinking
            if r_step > initial_value:
                step_right = r_step
            else:
                step_left = r_step
        initial_value = r_step
        
        # If we've finished the burn-in period
        if burn == 0: 
            # And we don't need to skip this
            if skip_count==skip:
                #print initial_value
                yield initial_value
                skip_count=0
            else: 
                skip_count = skip_count+1
        else:
            burn = burn-1
            
def take(n,g):
    """Take at most the first n elements of an iterable"""
    iter = g.__iter__()
    for _ in xrange(n):
            yield (iter.next())

E = sympy.symbols('E')
c = 3*10**8 #m/s
m = 0.511 #MeV/(c^2)
Q = 1.4743 #MeV

p_energy = sympy.Piecewise((0,E<0),(0,E>1.4743),(((E**2+2*E*m)**(1/2))*((Q-E)**2)*(E+m),True)) #C/c^5 is scaling factor roughly 1
number_of_samples = 10**5
output_csv_file = 'random_energies.csv' 

list_of_energies = take(number_of_samples,slice_sampling(p_energy,E))

with open(output_csv_file, 'w') as f:
    writer = csv.writer(f,dialect='excel',lineterminator='\n')
    for row in list_of_energies:  
        writer.writerow([row])