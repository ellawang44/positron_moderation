"""A module containing slice sampling
"""

import positron_constants
import random
import sympy
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
                if positron_constants.config.debug: print initial_value
                yield initial_value
                skip_count=0
            else: 
                skip_count = skip_count+1
        else:
            burn = burn-1