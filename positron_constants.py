"""Useful constants (and functions) for this project
"""

from __future__ import division
import sympy

E,z,rad = sympy.symbols('E z rad')
A = 4*10**(6) #g(ucm^-2)(keV^-n)
p = 3.21*10**18 #g(ucm^-3) can assume pure SiC density because radioactive source is only a small amount
n = 1.6 #taken from http://www.positronannihilation.net/index_files/Positron%20Beam.pdf section 4.1
m = 1.28 #http://iopscience.iop.org/article/10.1088/0953-8984/3/22/021/pdf
c = 3*10**8 #m/s
m_e = 0.511 #MeV/(c^2)
Q = 1.4743 #MeV
keV = 1000 #scaling factor needed for outputs to be in same units
#probability distribution function for positron emission energies
p_energy = sympy.Piecewise((0,E<0),(0,E>1.4743),(((E**2+2*E*m)**(1/2))*((Q-E)**2)*(E+m),True))
median_depth = A*(E**n)*(10**6)/p
zprime = 2.58*median_depth
#Makhovian equation for positron implantation depths in one dimension
#taken from http://iopscience.iop.org/article/10.1088/0953-8984/3/22/021/pdf
p_depth_energy = sympy.diff(
    (-(sympy.exp(-((z*10**(-6)/zprime)*(1+z*10**(-6)/zprime)**2)**m))), z)

m_universal,x = sympy.symbols('m_universal x')
scaling_factor = 10**7 #To deal with issues with scipy.opt.minimize_scalar's handling of small numbers
min_radius = scaling_factor*1.5*10**(-8) #Taken from Guessom 2005

#fragmentation size distribution equation
#size_distribution = sympy.Piecewise((0,rad<min_radius),(0,rad>max_radius),(rad**(alpha_value),True))

def p_depth(energy_value):    

    """Given an energy value, return a probability
    distrubution function for penetration depth
    """  
        
    return sympy.Piecewise((0,z<=0),(p_depth_energy,z>0)).subs(E,energy_value)
    
class Simulation_Config():
    def __init__(self):
        self.r_mul = 10 #radius multiplier
        self.r_exp = -7 #radius exponent
        self.alpha_value = -3/2 #fragmentation distribution, can be changed with set_alpha
        self.accuracy = 0.1
        self.debug = False
    # Treat max_radius and size_distribution as fields even though
    # they aren't really
    def __getattr__(self, name):
        if name == 'max_radius':
            return scaling_factor * self.r_mul * (10**self.r_exp)        
        elif name == 'size_distribution':
            return sympy.Piecewise((0,rad<min_radius),(0,rad>self.max_radius),(rad**(self.alpha_value),True))
        raise AttributeError
   
config = Simulation_Config()
