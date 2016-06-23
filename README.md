Orbitals
========

A program for students to explore time-dependent wavefuntions.

#####Author - Matthew Rowley

Released under an MIT license. Full license is contained in a file
named 'LICENSE'.

Running Orbitals
----------------
To run, first ensure that you have all the dependencies installed.
Then run:

``python Start_Orbitals.py``

###Dependencies

  * PyQt4

  * scipy

  * numpy

  * mayavi

  * pyface

  * traits

  * traitsui
  
  * pyqtgraph

To install dependencies, you may want to use a python package manager, such as pip. For example, try running the following (with the necessary privileges to make changes to your python installation directory):

``pip install PyQt4 scipy numpy mayavi pyface traits traitsui pyqtgraph``

Alternatively, you can navigate to the project home page of each dependency and install separately.

Rendering Orbitals in Three Dimensions
--------------------------------------
There are several ways to render atomic orbitals in three dimensions. The method used within the Orbitals application makes some reasonable sacrifices in exchange for dramatically faster computation time. Because of this, animations can be calculated and rendered in real-time. For pre-generating images or videos, however, one of the several other methods might be preferred. "Plot_Orbitals.py" includes code snippets for rendering atomic orbitals in many different ways, and can be a starting point for generating many types of images and videos. This code is based on the fine tutorials on the Mayavi website by Gael Varoquaux <gael.varoquaux@normalesup.org>.

The "scipy.special.sph_harm" function is useful as a general function (for arbitrarily large values of l and m), but it seems to be quite slow. The "Hydrogenic.py" module hard codes the most common spherical harmonics, and also includes the most common radial wavefunctions. It is used by both the Orbitals application and the "Plot_Orbitals.py" scripts for rapid and convenient access to atomic wavefunctions.
