# -*- coding: utf-8 -*-
"""
Hydrogenic.py
A module for representing hydrogenic electron orbitals.
Angular and radial functions are accessible separately, as well as radius
which encloses 90 percent of the electron density.
Note that to be compatable with scipy.special.sph_harm spherical coordinates
are designated according to the physics/chemistry convention rather than the
mathematics convention. i.e. theta is the angle of the projection on the xy
plane and takes values from 0 to 2*pi. phi is the angle with the z axis and
takes values from 0 to pi.

@author: Matthew B Rowley
"""
from scipy.special import sph_harm as sh
import numpy as np
from numpy import cos, sin, sqrt, pi, exp

sqrt2 = sqrt(2)
sqrt3 = sqrt(3)
sqrt5 = sqrt(5)
sqrt6 = sqrt(6)
sqrt30 = sqrt(30)
sqrtpi = sqrt(pi)

d_dict = {'z^2': lambda theta, phi: sqrt(5/pi)/4*(3*cos(phi)*cos(phi)-1),
          'xz': lambda theta, phi: sqrt(15/pi)/2*cos(theta)*sin(phi)*cos(phi),
          'yz': lambda theta, phi: sqrt(15/pi)/2*sin(theta)*sin(phi)*cos(phi),
          'x^2-y^2': lambda theta, phi: sqrt(15/pi)/4*cos(2*theta)*sin(phi)*sin(phi),
          'xy': lambda theta, phi: sqrt(15/pi)/4*sin(2*theta)*sin(phi)*sin(phi),
          2: lambda theta, phi: sqrt(15/2/pi)/4*exp(2j*theta)*sin(phi)*sin(phi),
          1: lambda theta, phi: -sqrt(15/2/pi)/2*exp(1j*theta)*sin(phi)*cos(phi),
          0: lambda theta, phi: sqrt(5/pi)/4*(3*cos(phi)*cos(phi)-1),
          -1: lambda theta, phi: sqrt(15/2/pi)/2*exp(-1j*theta)*sin(phi)*cos(phi),
          -2: lambda theta, phi: sqrt(15/2/pi)/4*exp(-2j*theta)*sin(phi)*sin(phi)
          }

f_dict = {'z^3': lambda theta, phi: sqrt(7/pi)/4*(5*cos(phi)*cos(phi)*cos(phi)-3*cos(phi)),
          'f1': lambda theta, phi: 1./sqrt2*(sh(1, 3, theta, phi)+sh(-1, 3, theta, phi)),
          'f-1': lambda theta, phi: 1./1j/sqrt2*(sh(1, 3, theta, phi)-sh(-1, 3, theta, phi)),
          'f2': lambda theta, phi: 1./sqrt2*(sh(2, 3, theta, phi)+sh(-2, 3, theta, phi)),
          'f-2': lambda theta, phi: 1./1j/sqrt2*(sh(2, 3, theta, phi)-sh(-2, 3, theta, phi)),
          'f3': lambda theta, phi: 1./sqrt2*(sh(3, 3, theta, phi)+sh(-3, 3, theta, phi)),
          'f-3': lambda theta, phi: 1./1j/sqrt2*(sh(3, 3, theta, phi)-sh(-3, 3, theta, phi)),
          3: lambda theta, phi: -sqrt(35/pi)/8*exp(3j*theta)*sin(phi)*sin(phi)*sin(phi),
          2: lambda theta, phi: sqrt(105/2/pi)/4*exp(2j*theta)*sin(phi)*sin(phi)*cos(phi),
          1: lambda theta, phi: -sqrt(21/pi)/8*exp(1j*theta)*sin(phi)*(5*cos(phi)*cos(phi)-1),
          0: lambda theta, phi: sqrt(7/pi)/4*(5*cos(phi)*cos(phi)*cos(phi)-3*cos(phi)),
          -1: lambda theta, phi: sqrt(21/pi)/8*exp(-1j*theta)*sin(phi)*(5*cos(phi)*cos(phi)-1),
          -2: lambda theta, phi: sqrt(105/2/pi)/4*exp(-2j*theta)*sin(phi)*sin(phi)*cos(phi),
          -3: lambda theta, phi: sqrt(35/pi)/8*exp(-3j*theta)*sin(phi)*sin(phi)*sin(phi)
          }

g_dict = {4: lambda theta, phi: sqrt(35/2/pi)*3/16*exp(4j*theta)*sin(phi)*sin(phi)*sin(phi)*sin(phi),
          3: lambda theta, phi: -sqrt(35/pi)*3/8*exp(3j*theta)*sin(phi)*sin(phi)*sin(phi)*cos(phi),
          2: lambda theta, phi: sqrt(5/2/pi)*3/8*exp(2j*theta)*sin(phi)*sin(phi)*(7*cos(phi)*cos(phi)-1),
          1: lambda theta, phi: -sqrt(5/pi)*3/8*exp(1j*theta)*sin(phi)*(7*cos(phi)*cos(phi)*cos(phi)-3*cos(phi)),
          0: lambda theta, phi: 3/16/sqrtpi*(35*cos(phi)*cos(phi)*cos(phi)*cos(phi)-30*cos(phi)*cos(phi)+3),
          -1: lambda theta, phi: sqrt(5/pi)*3/8*exp(-1j*theta)*sin(phi)*(7*cos(phi)*cos(phi)*cos(phi)-3*cos(phi)),
          -2: lambda theta, phi: sqrt(5/2/pi)*3/8*exp(-2j*theta)*sin(phi)*sin(phi)*(7*cos(phi)*cos(phi)-1),
          -3: lambda theta, phi: sqrt(35/pi)*3/8*exp(-3j*theta)*sin(phi)*sin(phi)*sin(phi)*cos(phi),
          -4: lambda theta, phi: sqrt(35/2/pi)*3/16*exp(-4j*theta)*sin(phi)*sin(phi)*sin(phi)*sin(phi)
          }

p_dict = {'x': lambda theta, phi: sqrt(3/pi)/2*sin(phi)*cos(theta),
          'y': lambda theta, phi: sqrt(3/pi)/2*sin(phi)*sin(theta),
          'z': lambda theta, phi: sqrt(3/pi)/2*cos(phi),
          1: lambda theta, phi: -sqrt(3/2/pi)/2*exp(1j*theta)*sin(phi),
          0: lambda theta, phi: sqrt(3/pi)/2*cos(phi),
          -1: lambda theta, phi: sqrt(3/2/pi)/2*exp(-1j*theta)*sin(phi)
          }

radials = {'10': lambda r, z: 2*z**(3./2)*np.exp(-z*r),
           '20': lambda r, z: 2*(z/2.)**(3./2)*(1-z*r/2)*np.exp(-z*r/2),
           '21': lambda r, z: 1/sqrt3*(z/2.)**(3./2)*z*r*np.exp(-z*r/2),
           '30': lambda r, z: 2*(z/3.)**(3./2)*(1-2*z*r/3+2*z*z*r*r/27)*np.exp(-z*r/3),
           '31': lambda r, z: 8./27/sqrt6*z**(3./2)*(z*r-z*z*r*r/6)*np.exp(-z*r/3),
           '32': lambda r, z: 4./81/sqrt30*z**(7./2)*r*r*np.exp(-z*r/3),
           '40': lambda r, z: 2,
           '41': lambda r, z: 2,
           '42': lambda r, z: 2,
           '43': lambda r, z: 2
           }


class Orbital:
    '''
    A class to represent a hydrogenic orbital
    '''
    def __init__(self, n, l, m, s=1, z=1):
        '''
        Initialize an instance of Orbital with quantum numbers for complex
        wavefunctions. "2px", for example, will have to be initialized
        '''
        self.n = n
        self.l = l
        self.s = s
        self.z = z
        self.m = m
        self.angular = get_angular(self.l, self.m, self.s)
        self.radial = get_radial(self.n, self.l, self.z)
        self.r_90p = get_90p(self.radial)
        self.psi = lambda r, theta, phi: (self.radial(r) *
                                          self.angular(theta, phi))


def get_angular(l, m, s):
    '''
    Return the angular wavefunction for the given l, m values.
    Performance is 60x faster for equations explicitly written in the
    dictionaries, but as a last resort just returning scipy.special.sph_harm
    will get the job done for arbitrarily high quantum numbers.
    '''
    if l == 0 or l == 's':
        angular = lambda theta, phi: np.ones_like(theta)*s*0.5/sqrtpi
    elif l == 1 or l == 'p':
        angular = lambda theta, phi: s*p_dict[m](theta, phi)
    elif l == 2 or l == 'd':
        angular = lambda theta, phi: s*d_dict[m](theta, phi)
    elif l == 3 or l == 'f':
        angular = lambda theta, phi: s*f_dict[m](theta, phi)
    elif l == 4 or l == 'g':
        angular = lambda theta, phi: s*g_dict[m](theta, phi)
    else:
        angular = lambda theta, phi: s * sh(m, l, theta, phi)
    return angular


def get_radial(n, l, z):
    '''
    Return the radial wavefunction for the given n, z, z values
    '''
    if n <= 4:
        # Typing in the functions explicitly allows speed optimization
        radial = lambda r: radials["{}{}".format(n, l)](r, z)
    else:
        # General form for arbitrarily large quantum numbers
        radial = lambda r: r
        # TODO http://quantummechanics.ucsd.edu/ph130a/130_notes/node237.html#derive:Hradial
    return radial


def get_90p(radial):
    '''
    Numerically integrate the radial function. Return radius
    which encloses 90%.
    '''
    r_range = np.arange(0, 30, 0.01)
    total = 0
    index = 0
    r_val = 0
    while(total < 0.9):
        r_val = r_range[index]
        y_val = radial(r_val)
        total = total + y_val*y_val*r_val*r_val/100
        index += 1
    return r_val

orbitals = {'1s': Orbital(1, 0, 0),
            '2s': Orbital(2, 0, 0),
            '2pz': Orbital(2, 1, 0),
            '2px': Orbital(2, 1, 'x'),
            '2py': Orbital(2, 1, 'y'),
            '2p1': Orbital(2, 1, 1),
            '2p-1': Orbital(2, 1, -1),
            '3s': Orbital(3, 0, 0),
            '3pz': Orbital(3, 1, 'z'),
            '3px': Orbital(3, 1, 'x'),
            '3py': Orbital(3, 1, 'y'),
            '3p1': Orbital(3, 1, 1),
            '3p-1': Orbital(3, 1, -1)}
