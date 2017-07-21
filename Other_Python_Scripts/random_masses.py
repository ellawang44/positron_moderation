""" Slice sampling for approximating the penetration depth of
positrons released from Sc-44
"""

from __future__ import division
import pylab as py
import sympy
import random
import csv
import numpy
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
                print initial_value
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

m_universal,x,r = sympy.symbols('m_universal x r')
make_my_code_run_fast = 10**7 #don't even ask... just trust me on this one
rho = 3210 #density of SiC in kg/m^3
#min_radius = make_my_code_run_fast*mass_to_radius.subs(x,min_mass)
min_radius = make_my_code_run_fast*10**(-8)
#1.5*10**(-8) Guessom 2005 alpha = -3/2
#0.2*10**(-9) arbitary, may need to be increased to 1-2 because lattice spacing of SiC
max_radius = make_my_code_run_fast*0.25*10**(-6)
radius_range = numpy.arange(0,max_radius,10**(2)) #this was just used to plot universal size distribution
F = sympy.Piecewise((0,r<min_radius),(0,r>max_radius),(r**(-3/2),True))
number_of_samples = 10**7
output_csv_file = 'random_masses.csv' 

list_of_masses = (((float(i)/make_my_code_run_fast) for i in list(take(number_of_samples,slice_sampling(F,r,width=10**(-1))))))

with open(output_csv_file, 'w') as f:
    writer = csv.writer(f,dialect='excel',lineterminator='\n')
    for row in list_of_masses:  
        writer.writerow([row])