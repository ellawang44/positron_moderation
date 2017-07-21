import struct
import pylab
import math
from helper_functions.folds import combine, select, when, before, collect_list, histogram

binary_file = 'MCsim 0.1keV accuracy 1x10^-6 RMax -1.5 alpha 10^7 samples trial 1.bin'

def read_bin_file(filename):
    with open(filename,'rb')as f:
        while True:
            data = f.read(struct.calcsize('fffb'))
            if len(data) != struct.calcsize('fffb'): return 
            yield struct.unpack_from('fffb',data)
            
#positrons are ejected with a positive energy because of properties of SiC
positive_work_function = 0.0022*1000

#annihilation data is in 4th column, so select that column of data to get information on annihilation
get_annihilated = before(collect_list(),select(3))
#gets the energies of the escaped positrons
get_energies = when(lambda x: x[3] ==1, before(collect_list(), select(0)))
#sets the energy of the annihilated positrons to the positive work function
get_annihilated_data = when(lambda x: x[3]==0, before(collect_list(), lambda _: positive_work_function))
#builds a function to process the data
processing_function = combine(get_annihilated, get_energies, get_annihilated_data)
#runs the function
annihilated,energies,annihilated_data = processing_function(read_bin_file(binary_file))
no_of_bins = 100

bins = pylab.linspace(0,1.5,no_of_bins)
annihilated_scaling = 100
annihilated_bins = pylab.linspace(0,1.5,no_of_bins*annihilated_scaling) #make thinner

pylab.hist(energies,bins) #histograms the escaped positrons based on their energy
pylab.hist(annihilated_data*annihilated_scaling #make taller
    ,annihilated_bins,facecolor='g') #histograms the annihilated positrons 
pylab.axes(xscale='log',yscale='log') #comment out this line for non-log scale
pylab.xlabel("Energy (keV)")
pylab.ylabel("no. of occurances")
pylab.title("energy of emitted positrons with annihilated")
pylab.show()

#prints the number of annihilated positrons
print len(annihilated_data)