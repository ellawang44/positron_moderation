import pylab

data = pylab.loadtxt('random_depths.csv',delimiter=',')
energies = data[:,0]
depths = data[:,1]

#if len(data.shape) > 1:
#    depths = data[:,0]
#else:
#    depths = data

pylab.hist(depths*10**(-3),bins=100)
pylab.xlabel("depth (mcm)")
pylab.ylabel("number of occurances")
pylab.title("random depths generated with proper energy distribution function for Sc-44")
pylab.show()