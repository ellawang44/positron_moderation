from __future__ import division
import pylab
import numpy
import sympy
import math

m_universal,x,r = sympy.symbols('m_universal x r')
scaling_factor = 10**7 #To deal with issues with scipy.opt.minimize_scalar's handling of small numbers
rho = 3210 #density of SiC in kg/m^3
min_mass = 4.6637066*10**(-26)+1.9944235*10**(-26) #kg, mass of C and Si added together
max_mass = 10**(-13) #kg not sure why this is here, it's not even being used
mass_to_radius = (3*x/(4*rho*math.pi))**(1/3)
#min_radius = make_my_code_run_fast*mass_to_radius.subs(x,min_mass)
min_radius = scaling_factor*1.5*10**(-8) #Taken from Guessom 2005
max_radius = scaling_factor*3.0*10**(-7) #varies
radius_range = numpy.arange(0,max_radius+0.1,10**(-2)) #this was just used to plot universal size distribution
F = sympy.Piecewise((0,r<min_radius),(0,r>max_radius),((-3/2)*r**(-5.0/2.0),True))

pylab.plot(list(radius_range),list(F.subs(r,i) for i in radius_range))
pylab.title('dN/dR as radius varies \n alpha = -3/2')
pylab.xlabel('radius (nm)')
pylab.ylabel('dN/dR')
pylab.show()