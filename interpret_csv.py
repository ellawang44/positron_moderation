import pylab

data = pylab.loadtxt('MCsim 0.1keV accuracy 10^-6 Rmax -1.5 alpha 10^7 samples.csv',delimiter=",")

energy = data[:,0]
depth = data[:,1]
radius = data[:,2]
escape = data[:,3]

annihilated = 0
for a in escape:
    if a == 0:
        annihilated += 1
    
print annihilated