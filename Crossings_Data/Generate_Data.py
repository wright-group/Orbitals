# -*- coding: utf-8 -*-
# This program is licences under an MIT license. Full licence is at the end of
# this file.
"""
Orbitals_UI.py
A user interface for generating and interacting with various animations of
hydrogenic wavefunctions.

@author: Matthew B Rowley
"""
from __future__ import division
import pickle
import numpy as np

phi, theta = np.mgrid[0:np.pi:20j, 0:2*np.pi:40j]
r = 0.75
x = r * np.sin(phi) * np.cos(theta)
y = r * np.sin(phi) * np.sin(theta)
z = r * np.cos(phi)

bond_lengths = np.linspace(1.5,5,100)

Low_A_Curve = []
High_A_Curve = []
Low_D_Curve = []
High_D_Curve = []
Low_A_MO = []
High_A_MO = []
Low_D_MO = []
High_D_MO = []

for length in bond_lengths:
    # Some dummy curves so I can get the GUI ready
    Low_A_Curve.append(20 + length)
    High_A_Curve.append(30 + length)
    Low_D_Curve.append(21 - length)
    High_D_Curve.append(29 - length)
    # Some dummy MOs, too
    Low_A_MO.append((x, y + 1, z))
    High_A_MO.append((x, y - 1, z + length))
    Low_D_MO.append((x, y - 1, z))
    High_D_MO.append((x, y + 1, z + length))

Low_A_Curve = np.array(Low_A_Curve)
High_A_Curve = np.array(High_A_Curve)
Low_D_Curve = np.array(Low_D_Curve)
High_D_Curve = np.array(High_D_Curve)
Low_A_MO = np.array(Low_A_MO)
High_A_MO = np.array(High_A_MO)
Low_D_MO = np.array(Low_D_MO)
High_D_MO = np.array(High_D_MO)

pickled_data = {'bond_lengths': bond_lengths,
                'Low_A_Curve': Low_A_Curve,
                'High_A_Curve': High_A_Curve,
                'Low_D_Curve': Low_D_Curve,
                'High_D_Curve': High_D_Curve,
                'Low_A_MO': Low_A_MO,
                'High_A_MO': High_A_MO,
                'Low_D_MO': Low_D_MO,
                'High_D_MO': High_D_MO}

with open("pickled_data.p", "wb") as file:
    pickle.dump(pickled_data, file)

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
