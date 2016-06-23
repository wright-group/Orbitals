# -*- coding: utf-8 -*-
# This program is licenced under an MIT license. Full licence is at the end of
# this file.
"""
Hydrogenic.py
A module for representing hydrogenic electron orbitals.
Angular and radial functions are accessible separately, as well as radius
which encloses 90 percent of the electron density integrating along radial
lines.

Note that to be compatable with scipy.special.sph_harm spherical coordinates
are designated according to the physics/chemistry convention rather than the
mathematics convention. i.e. theta is the angle of the projection on the xy
plane and takes values from 0 to 2*pi. phi is the angle with the z axis and
takes values from 0 to pi.

@author: Matthew B Rowley
"""
from __future__ import division
from scipy.special import sph_harm as sh
import numpy as np
from numpy import cos, sin, sqrt, pi, exp

sqrt2 = sqrt(2)
sqrt3 = sqrt(3)
sqrt5 = sqrt(5)
sqrt6 = sqrt(6)
sqrt7 = sqrt(7)
sqrt15 = sqrt(15)
sqrt21 = sqrt(21)
sqrt24 = sqrt(24)
sqrt27 = sqrt(27)
sqrt30 = sqrt(30)
sqrt35 = sqrt(35)
sqrt105 = sqrt(105)
sqrtpi = sqrt(pi)

d_dict = {'z^2': lambda theta, phi: sqrt5/sqrtpi/4*(3*cos(phi)*cos(phi)-1),
          'xz': lambda theta, phi: sqrt15/sqrtpi/2*cos(theta)*sin(phi)*cos(phi),
          'yz': lambda theta, phi: sqrt15/sqrtpi/2*sin(theta)*sin(phi)*cos(phi),
          'x^2-y^2': lambda theta, phi: sqrt15/sqrtpi/4*cos(2*theta)*sin(phi)*sin(phi),
          'xy': lambda theta, phi: sqrt15/sqrtpi/4*sin(2*theta)*sin(phi)*sin(phi),
          2: lambda theta, phi: sqrt15/sqrt2/sqrtpi/4*exp(2j*theta)*sin(phi)*sin(phi),
          1: lambda theta, phi: -sqrt15/sqrt2/sqrtpi/2*exp(1j*theta)*sin(phi)*cos(phi),
          0: lambda theta, phi: sqrt5/sqrtpi/4*(3*cos(phi)*cos(phi)-1),
          -1: lambda theta, phi: sqrt15/sqrt2/sqrtpi/2*exp(-1j*theta)*sin(phi)*cos(phi),
          -2: lambda theta, phi: sqrt15/sqrt2/sqrtpi/4*exp(-2j*theta)*sin(phi)*sin(phi)
          }

# alternate F orbital names:
# y^3-3yx^2
# x^3-3xy^2
# 5yz^2-yr^2
# 5xz^2-3xr^2
# zx^2-zy^2
# xyz
# 5z^3-3zr^2
# TODO: Simplify real f orbitals.
f_dict = {'z^3': lambda theta, phi: sqrt7/sqrtpi/4*(5*cos(phi)*cos(phi)*cos(phi)-3*cos(phi)),
          'xz^2': lambda theta, phi: 1/1j/sqrt2*(-sqrt21/sqrtpi/8*exp(1j*theta)*sin(phi)*(5*cos(phi)*cos(phi)-1)+sqrt21/sqrtpi/8*exp(-1j*theta)*sin(phi)*(5*cos(phi)*cos(phi)-1)),
          'yz^2': lambda theta, phi: 1/sqrt2*(-sqrt21/sqrtpi/8*exp(1j*theta)*sin(phi)*(5*cos(phi)*cos(phi)-1)-sqrt21/sqrtpi/8*exp(-1j*theta)*sin(phi)*(5*cos(phi)*cos(phi)-1)),
          'xyz': lambda theta, phi: 1/sqrt2*(sqrt105/sqrt2/sqrtpi/4*exp(2j*theta)*sin(phi)*sin(phi)*cos(phi)+sqrt105/sqrt2/sqrtpi/4*exp(-2j*theta)*sin(phi)*sin(phi)*cos(phi)),
          'z(x^2-y^2)': lambda theta, phi: 1/1j/sqrt2*(sqrt105/sqrt2/sqrtpi/4*exp(2j*theta)*sin(phi)*sin(phi)*cos(phi)-sqrt105/sqrt2/sqrtpi/4*exp(-2j*theta)*sin(phi)*sin(phi)*cos(phi)),
          'x(x^2-3y^2)': lambda theta, phi: 1/1j/sqrt2*(-sqrt35/sqrtpi/8*exp(3j*theta)*sin(phi)*sin(phi)*sin(phi)+sqrt35/sqrtpi/8*exp(-3j*theta)*sin(phi)*sin(phi)*sin(phi)),
          'y(3x^3-y^2)': lambda theta, phi: 1/sqrt2*(-sqrt35/sqrtpi/8*exp(3j*theta)*sin(phi)*sin(phi)*sin(phi)-sqrt35/sqrtpi/8*exp(-3j*theta)*sin(phi)*sin(phi)*sin(phi)),
          3: lambda theta, phi: -sqrt35/sqrtpi/8*exp(3j*theta)*sin(phi)*sin(phi)*sin(phi),
          2: lambda theta, phi: sqrt105/sqrt2/sqrtpi/4*exp(2j*theta)*sin(phi)*sin(phi)*cos(phi),
          1: lambda theta, phi: -sqrt21/sqrtpi/8*exp(1j*theta)*sin(phi)*(5*cos(phi)*cos(phi)-1),
          0: lambda theta, phi: sqrt7/sqrtpi/4*(5*cos(phi)*cos(phi)*cos(phi)-3*cos(phi)),
          -1: lambda theta, phi: sqrt21/sqrtpi/8*exp(-1j*theta)*sin(phi)*(5*cos(phi)*cos(phi)-1),
          -2: lambda theta, phi: sqrt105/sqrt2/sqrtpi/4*exp(-2j*theta)*sin(phi)*sin(phi)*cos(phi),
          -3: lambda theta, phi: sqrt35/sqrtpi/8*exp(-3j*theta)*sin(phi)*sin(phi)*sin(phi)
          }

g_dict = {4: lambda theta, phi: sqrt35/sqrt2/sqrtpi*3/16*exp(4j*theta)*sin(phi)*sin(phi)*sin(phi)*sin(phi),
          3: lambda theta, phi: -sqrt35/sqrtpi*3/8*exp(3j*theta)*sin(phi)*sin(phi)*sin(phi)*cos(phi),
          2: lambda theta, phi: sqrt5/sqrt2/sqrtpi*3/8*exp(2j*theta)*sin(phi)*sin(phi)*(7*cos(phi)*cos(phi)-1),
          1: lambda theta, phi: -sqrt5/sqrtpi*3/8*exp(1j*theta)*sin(phi)*(7*cos(phi)*cos(phi)*cos(phi)-3*cos(phi)),
          0: lambda theta, phi: 3/16/sqrtpi*(35*cos(phi)*cos(phi)*cos(phi)*cos(phi)-30*cos(phi)*cos(phi)+3),
          -1: lambda theta, phi: sqrt5/sqrtpi*3/8*exp(-1j*theta)*sin(phi)*(7*cos(phi)*cos(phi)*cos(phi)-3*cos(phi)),
          -2: lambda theta, phi: sqrt5/sqrt2/sqrtpi*3/8*exp(-2j*theta)*sin(phi)*sin(phi)*(7*cos(phi)*cos(phi)-1),
          -3: lambda theta, phi: sqrt35/sqrtpi*3/8*exp(-3j*theta)*sin(phi)*sin(phi)*sin(phi)*cos(phi),
          -4: lambda theta, phi: sqrt35/sqrt2/sqrtpi*3/16*exp(-4j*theta)*sin(phi)*sin(phi)*sin(phi)*sin(phi)
          }

p_dict = {'x': lambda theta, phi: sqrt3/sqrtpi/2*sin(phi)*cos(theta),
          'y': lambda theta, phi: sqrt3/sqrtpi/2*sin(phi)*sin(theta),
          'z': lambda theta, phi: sqrt3/sqrtpi/2*cos(phi),
          1: lambda theta, phi: -sqrt3/sqrt2/sqrtpi/2*exp(1j*theta)*sin(phi),
          0: lambda theta, phi: sqrt3/sqrtpi/2*cos(phi),
          -1: lambda theta, phi: sqrt3/sqrt2/sqrtpi/2*exp(-1j*theta)*sin(phi)
          }

# TODO: 4x dont have z dependence and 5x are totally wrong
radials = {'10': lambda r, z: 2*z**(3/2)*np.exp(-z*r),
           '20': lambda r, z: 2*(z/2)**(3/2)*(1-z*r/2)*np.exp(-z*r/2),
           '21': lambda r, z: 1/sqrt3*(z/2)**(3/2)*z*r*np.exp(-z*r/2),
           '30': lambda r, z: 2*(z/3)**(3/2)*(1-2*z*r/3+2*z*z*r*r/27)*np.exp(-z*r/3),
           '31': lambda r, z: 8/27/sqrt6*z**(3/2)*(z*r-z*z*r*r/6)*np.exp(-z*r/3),
           '32': lambda r, z: 4/81/sqrt30*z**(7/2)*r*r*np.exp(-z*r/3),
           '40': lambda r, z: 1/4*(1-3*r/4+r*r/8-r*r*r/192)*np.exp(-r/4),
           '41': lambda r, z: sqrt5/16/sqrt3*(1-r/4+r*r/80)*r*np.exp(-r/4),
           '42': lambda r, z: r*r/64/sqrt5*(1-r/12)*np.exp(-r/4),
           '43': lambda r, z: r*r*r/768/sqrt35*np.exp(-r/4),
           '50': lambda r, z: 2,
           '51': lambda r, z: 2,
           '52': lambda r, z: 2,
           '53': lambda r, z: 2,
           '54': lambda r, z: 2
           }


class Orbital:
    '''
    A class to represent a hydrogenic orbital
    '''
    def __init__(self, n, l, m, s=1, z=1, bohr=1):
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
        self.bohr = bohr

    def setBohr(self, bohr):
        '''Define the bohr oscillation'''
        self.bohr = bohr


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
    Return the radial wavefunction for the given n, z, z values.
    Performance is improved for the equations which are explicitly written
    in the dictionary, but as a last resort the Lagrange polynomial is
    generated programmatically.
    '''
    if n <= 5:
        # Typing in the functions explicitly allows speed optimization
        radial = lambda r: radials["{}{}".format(n, l)](r, z)
    elif l == n - 1:
        print("L is N-1")
        radial = lambda r: r**(n-1)*np.exp(-z*r/n)
        rvals = np.arange(0,200,0.01)
        integral = 0
        for r in rvals:
            radial_val = radial(r)
            integral = integral + r*r*radial_val*radial_val/100
        radial = lambda r: r**(n-1)*np.exp(-z*r/n)/integral
    else:
        # General form for arbitrarily large quantum numbers
        radial = lambda r: 2
        # TODO http://quantummechanics.ucsd.edu/ph130a/130_notes/node237.html#derive:Hradial
    return radial


def get_90p(radial):
    '''
    Numerically integrate the radial function. Return radius
    which encloses 90%.
    '''
    r_range = np.arange(0, 50, 0.01)
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
            '3p-1': Orbital(3, 1, -1),
            '3dz^2': Orbital(3, 2, 0),
            '3dxz': Orbital(3, 2, 'xz'),
            '3dyz': Orbital(3, 2, 'yz'),
            '3dx^2-y^2': Orbital(3, 2, 'x^2-y^2'),
            '3dxy': Orbital(3, 2, 'xy'),
            '4fz^3': Orbital(4,3, 'z^3'),
            '4fxz^2': Orbital(4, 3, 'xz^2'),
            '4fyz^2': Orbital(4, 3, 'yz^2'),
            '4fxyz': Orbital(4, 3, 'xyz'),
            '4fz(x^2-y^2)': Orbital(4, 3, 'z(x^2-y^2)'),
            '4fx(x^2-3y^2)': Orbital(4, 3, 'x(x^2-3y^2)'),
            '4fy(3x^2-y^2)': Orbital(4, 3, 'y(3x^3-y^2)')}

# The MIT License (MIT)
#
# Copyright (c) 2015 Matthew B. Rowley
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
