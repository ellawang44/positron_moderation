import pylab as py

data = py.loadtxt("ASC_python.csv", delimiter=",")
E = data[:,0]

A = 4*10**6 #g/ucm^2keV^n
p = 3.21*10**18 #g/ucm^3 can assume pure SiC because radioactive source is only a small amount
n = 1.6
median_depth = A*(E**n)/p

py.plot (E,median_depth)
py.xlabel("energy (keV)")
py.ylabel("median implantation depth (ucm)")
py.title("median implantation depth as energy varies \n(A=4*10^6 g ucm^-2 keV^-n, n=1.6, p=3.21*10^18 g ucm^-3", y=1.03)
py.show()