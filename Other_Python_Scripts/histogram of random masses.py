import pylab

data = pylab.loadtxt('random_masses.csv',delimiter=',')

if len(data.shape) > 1:
    masses = data[:,0]
else:
    masses = data

pylab.hist(masses,bins=1000)
pylab.xlabel("radius (m)")
pylab.ylabel("number of occurances")
pylab.title("random masses generated with alpha value (-3/2)")
pylab.show()