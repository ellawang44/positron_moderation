import pylab
import scipy.stats
from helper_functions import genzip


data = pylab.loadtxt('MCsim 10^6.csv',delimiter=",")

energy = data[:,0]
depth = data[:,1]
radius = data[:,2]
escape = data[:,3]

most_common_radius = (scipy.stats.mode(map(lambda x : round(x,12),radius)).mode)[0]
energy_values = [en for en,rad in genzip(energy,radius) if round(rad,12) == most_common_radius]

pylab.hist(energy_values,bins=100)
pylab.xlabel("energy (MeV)")
pylab.ylabel("no. of occurances")
pylab.show()

print most_common_radius
print len(energy_values)