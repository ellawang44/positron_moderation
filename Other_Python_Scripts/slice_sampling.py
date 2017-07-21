""" Slice sampling for approximating the penetration depth of
positrons released from Sc-44
"""

from __future__ import division
import pylab as py
import sympy
import random
import csv
import math
import numpy
import scipy.optimize as opt

def slice_sampling(probability_function,substitute,initial_value=None,
                   width=100,burn=300,skip=0, magnitude=1.0):
    
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

    # Generate python function
    p = sympy.lambdify(substitute,probability_function,"mpmath")
    
    # Calculate a good initial_value if none is given
    if initial_value is None:
        min_val = opt.minimize_scalar(lambda a:((1/magnitude) * -p(a)))
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
    giter = g.__iter__()
    for _ in xrange(n):
            yield (giter.next())

def genzip(g1,g2):
    """Zip a pair of iterables returning a generator"""
    giter1 = g1.__iter__()
    giter2 = g2.__iter__()
    while True:
        yield (giter1.next(),giter2.next())

E,z = sympy.symbols('E z')
A = 4*10**(6) #g/ucm^2keV^n
p = 3.21*10**18 #g/ucm^3 can assume pure SiC because radioactive source is only a small amount
n = 1.6
m = 1.28
c = 3*10**8 #m/s
m_e = 0.511 #MeV/(c^2)
Q = 1.4743 #MeV
keV=1000
p_energy = sympy.Piecewise((0,E<0),(0,E>1.4743),(((E**2+2*E*m)**(1/2))*((Q-E)**2)*(E+m),True)) #C/c^5 is scaling factor roughly 1
median_depth = A*(E**n)*(10**6)/p
zprime = 2.58*median_depth
p_depth_noE = sympy.diff(
    (-(sympy.exp(-((z*10**(-6)/zprime)*(1+z*10**(-6)/zprime)**2)**m))), z)
number_of_samples = 10**5
output_csv_file = 'Monte_Carlo_simulation rounding to 0.1keV.csv' 

m_universal,x,r = sympy.symbols('m_universal x r')
make_my_code_run_fast = 10**9
rho = 3210
min_mass = 4.6637066*10**(-26)+1.9944235*10**(-26)
max_mass = 10**(-13)
mass_to_radius = (3*x/(4*rho*math.pi))**(1/3)
#min_radius = make_my_code_run_fast*mass_to_radius.subs(x,min_mass)
min_radius = make_my_code_run_fast*0.2*10**(-9) #may need to be increased to 1-2 because lattice spacing of SiC
max_radius = make_my_code_run_fast*10**(-4)
radius_range = numpy.arange(0,max_radius,10**(2))
F = sympy.Piecewise((0,r<min_radius),(0,r>max_radius),((4*r**3*math.pi*rho/3)**(-3/2),True))

def p_depth(energy_value):    

    """Given an energy value, return a probability
    distrubution function for penetration depth
    """  
        
    return sympy.Piecewise((0,z<=0),(p_depth_noE,z>0)).subs(E,energy_value)

if __name__ == "__main__":
    mass_generator = slice_sampling(F,r,width=10)
    energy_generator = slice_sampling(p_energy,E)
    depth_functions = [slice_sampling((p_depth(generator)),z,width=10000) for generator in numpy.arange(0.1,1474.3,0.1)]
    list_of_depths = [(lambda d: (e,(d*10**(-8)),(m*10**(-9)),0 if (d*10**(-8)) < (m*10**(-9)) else 1))(depth_functions[math.ceil(e*keV*10)/10].next()) 
                        for m,e in take(number_of_samples,genzip(mass_generator,energy_generator))]
                    #list indices must be integers, not float - 
                    #trying to increase accuracy of ceiling energy to 0.1keV as opposed to 1keV, not sure what I broke

    with open(output_csv_file, 'w') as f:
        writer = csv.writer(f,dialect='excel',lineterminator='\n')
        for row in list_of_depths:
            writer.writerow(row)