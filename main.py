"""The main entry point for this program"""

from __future__ import division
import csv
from depth_functions import DepthFunctions
from itertools import izip, islice
from positron_constants import * # sorry python :'(
from slice_sampling import slice_sampling
import struct
import time

s_exp = 5 #exponent for the number of samples
number_of_samples = 10**s_exp

def generate_file_name(n):
    """Generates an informative file name that doesn't clash with python or file naming syntax 
    while making this function unreadable in the process"""
    return ''.join(['MCsim ',str(config.accuracy),'keV accuracy ',str(config.r_mul),'x10^',str(config.r_exp),
        ' RMax ', str(config.alpha_value), ' alpha 10^',str(s_exp),' samples trial ',str(n),'.',output_type])

output_type = 'bin'

def run_test(trial_number):
    """runs the Monte Carlo simulation once"""
    start_time = time.time()
    radius_generator = slice_sampling(config.size_distribution,rad,width=10**(1))
    energy_generator = slice_sampling(p_energy,E)
    depth_functions = DepthFunctions(config.accuracy) #accuracy of rounding
    depths = ((e,(d*10**(-8)),(r/scaling_factor),0 if (d*10**(-8)) < (r/scaling_factor) else 1) 
                        for r,e in islice(izip(radius_generator,energy_generator),number_of_samples)
                        for d in (depth_functions[e*keV],))
    if output_type == 'csv':
        with open(generate_file_name(trial_number+1), "w") as f:
            writer = csv.writer(f,dialect='excel',lineterminator='\n')
            for row in depths:
                writer.writerow(row)
    else:
        with open(generate_file_name(trial_number+1), 'wb') as f:
            for e,d,r,p in depths:
                f.write(struct.pack('fffb', e,d,r,bool(p)))
    
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == "__main__":
    config.accuracy = 0.1
    number_of_trials = 1 #runs this number of trials
    config.debug = False #prints the values that it's generating to terminal, use "False" if you want it to stop, 
                     #helps debug (makes sure the code is actually running)
    config.r_mul = 1
    config.r_exp = -6
    #is equal to radius*10**radius_exponent
    #file naming issues
    config.alpha_value = -3/2
    trials = xrange(number_of_trials)
    for trial in trials:
        run_test(trial)