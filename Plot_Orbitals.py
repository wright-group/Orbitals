# Much of the core functionality in these scripts comes from
# the fine tutorial demonstrations on the Mayavi website.
# I've included the original boilerplate for attribution
# below.

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# Copyright (c) 2008, Enthought, Inc.
# License: BSD Style.

# Create the data ############################################################
#from __future__ import division
import numpy as np
import Hydrogenic as hyd
from mayavi import mlab

x, y, z = np.mgrid[- 20:20:150j, - 20:20:150j, - 20:20:150j]
r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
phi = np.arccos(z/r)
theta = np.arctan(y/x)

orbital = hyd.Orbital(3, 1, 0, 1)
L = orbital.radial
A = orbital.angular
Psi = L(r) * A(theta, phi)
Phi = (np.conj(Psi)*Psi).real
phase = np.angle(Psi)

x2, y2, z2 = np.mgrid[- 20:20:150j, - 0.1:20:75j, - 20:20:150j]
r2 = np.sqrt(x2 ** 2 + y2 ** 2 + z2 ** 2)
phi2 = np.arccos(z2/r2)
theta2 = np.arctan(y2/x2)

Psi2 = L(r2) * A(theta2, phi2)
Phi2 = (np.conj(Psi2)*Psi2).real
phase2 = np.angle(Psi2)

if False: # Iso phase coloring
    # Plot it ####################################################################
    
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    # We create a scalar field with the module of Phi as the scalar
    src = mlab.pipeline.scalar_field(x, y, z, Phi)

    # And we add the phase of Phi as an additional array
    # This is a tricky part: the layout of the new array needs to be the same
    # as the existing dataset, and no checks are performed. The shape needs
    # to be the same, and so should the data. Failure to do so can result in
    # segfaults.
    src.image_data.point_data.add_array(phase.T.ravel())
    # We need to give a name to our new dataset.
    src.image_data.point_data.get_array(1).name = 'angle'
    # Make sure that the dataset is up to date with the different arrays:
    src.update()

    # We select the 'scalar' attribute, ie the norm of Phi
    src2 = mlab.pipeline.set_active_attribute(src)#,
                                        #point_scalars='scalar')

    # Cut isosurfaces of the norm
    contour = mlab.pipeline.contour(src2)
    contour.filter.contours= [0.00005]

    # Now we select the 'angle' attribute, ie the phase of Phi
    contour2 = mlab.pipeline.set_active_attribute(contour,
                                        point_scalars='angle')

    # And we display the surface. The colormap is the current attribute: the phase.
    mlab.pipeline.surface(contour2, colormap='hsv', vmax=np.pi, vmin=-np.pi)
    mlab.colorbar(title='Phase', orientation='vertical', nb_labels=5)
    mlab.view(-10, 90)
    mlab.show()


if False: # iso density colors
    # Plot it ####################################################################
    
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    contour = mlab.contour3d(x, y, z, Phi, colormap='hsv')
    source = contour.mlab_source
    source2 = mlab.pipeline.scalar_field(x, y, z, phase)
    #contour2 = mlab.pipeline.set_active_attribute(contour, source2)                                                 
    #mlab.pipeline.surface(contour2, colormap='hsv', vmax=np.pi, vmin=-np.pi)
    mlab.colorbar(title='Density', orientation='vertical', nb_labels=5)
    mlab.view(-10, 90)
    mlab.show()    
    
if False: # volume
    # Plot it ####################################################################
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    # We create a scalar field with the module of Phi as the scalar
    src = mlab.pipeline.scalar_field(x, y, z, Phi)

    # Cut isosurfaces of the norm
    volume1 = mlab.pipeline.volume(src, vmin=0)
    #volume2 = mlab.pipeline.volume(src, vmax=-0.008)
    src3 = mlab.pipeline.scalar_field(x, y, z, Phi)
    #mlab.pipeline.image_plane_widget(src3, plane_orientation='x_axes', slice_index=100, colormap='hsv', vmax=np.pi, vmin=-np.pi)
    mlab.colorbar(title='Density', orientation='vertical', nb_labels=5, label_fmt='%.0e')
    mlab.view(-10, 90)
    mlab.show()
    
if False: # Iso and volume
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    # iso
    src = mlab.pipeline.scalar_field(x2, y2, z2, np.abs(Phi2))
    src.image_data.point_data.add_array(phase2.T.ravel())
    src.image_data.point_data.get_array(1).name = 'angle'
    src.update()
    src2 = mlab.pipeline.set_active_attribute(src)
    contour = mlab.pipeline.contour(src2)
    contour.filter.contours= [0.00005]
    contour2 = mlab.pipeline.set_active_attribute(contour, point_scalars='angle')
    mlab.pipeline.surface(contour2, colormap='hsv', vmax=np.pi, vmin=-np.pi)
    src3 = mlab.pipeline.scalar_field(x, y, z, Phi)
    volume1 = mlab.pipeline.volume(src3, vmin=0)
    mlab.view(-10, 90)
    mlab.show()
    
if False: # iso + Scaled Sphere
    #-------------------------------------------------------------------------------------
    #iso
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    src = mlab.pipeline.scalar_field(x2, y2, z2, np.abs(Phi2))
    src.image_data.point_data.add_array(phase2.T.ravel())
    src.image_data.point_data.get_array(1).name = 'angle'
    src.update()
    src2 = mlab.pipeline.set_active_attribute(src)
    contour = mlab.pipeline.contour(src2)
    contour.filter.contours= [0.00005]
    contour2 = mlab.pipeline.set_active_attribute(contour, point_scalars='angle')
    mlab.pipeline.surface(contour2, colormap='hsv', vmax=np.pi, vmin=-np.pi)
    # Scaled Sphere
    phis, thetas = np.mgrid[0:np.pi:80j, np.pi:2*np.pi:160j]
    rs = orbital.r_90p
    xs = rs * np.sin(phis) * np.cos(thetas)
    ys = rs * np.sin(phis) * np.sin(thetas)
    zs = rs * np.cos(phis)
    psis = orbital.angular(thetas, phis)
    rs = (np.conj(psis) * psis).real
    rs = rs / np.max(rs)
    xs, ys, zs = rs*xs, rs*ys, rs*zs
    angles = np.abs(np.angle(psis))-np.pi
    mesh = mlab.mesh(xs, ys, zs, scalars=angles, colormap='hsv', vmax=np.pi, vmin=-np.pi)
    # This low res colormap is simpler and less distracting
    my_map = [[255, 0, 0, 255], [140, 0, 140, 255], [140, 0, 140, 255],
              [0, 0, 255, 255], [0, 0, 255, 255], [255, 153, 18, 255],
              [255, 153, 18, 255], [255, 0, 0, 255]]
    # Comment out the next line to use a smoother high res colormap
    #mesh.module_manager.scalar_lut_manager.lut.table = my_map
    #mlab.colorbar(title='Phase', orientation='vertical', nb_labels=5)
    mlab.view(-10, 90)
    mlab.show()
if False: # Scaled Sphere
    #-------------------------------------------------------------------------------------
    phis, thetas = np.mgrid[0:np.pi:80j, 0:2*np.pi:160j]
    rs = orbital.r_90p
    xs = rs * np.sin(phis) * np.cos(thetas)
    ys = rs * np.sin(phis) * np.sin(thetas)
    zs = rs * np.cos(phis)
    psis = orbital.angular(thetas, phis)
    rs = (np.conj(psis) * psis).real
    rs = rs / np.max(rs)
    xs, ys, zs = rs*xs, rs*ys, rs*zs
    angles = np.abs(np.angle(psis))-np.pi
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    mesh = mlab.mesh(xs, ys, zs, scalars=angles, colormap='hsv', vmax=np.pi, vmin=-np.pi)
    # This low res colormap is simpler and less distracting
    my_map = [[255, 0, 0, 255], [140, 0, 140, 255], [140, 0, 140, 255],
              [0, 0, 255, 255], [0, 0, 255, 255], [255, 153, 18, 255],
              [255, 153, 18, 255], [255, 0, 0, 255]]
    # Comment out the next line to use a smoother high res colormap
    #mesh.module_manager.scalar_lut_manager.lut.table = my_map
    mlab.colorbar(title='Phase', orientation='vertical', nb_labels=5)
    mlab.view(0, 90)
    mlab.show()
#-------------------------------------------------------------------------------------
if False: # points
    num = 400j
    xp, yp, zp = np.ogrid[- 35:35:num, - 35:35:num, - 35:35:num]
    rp = np.sqrt(xp ** 2 + yp ** 2 + zp ** 2)
    phip = np.arccos(zp/rp)
    thetap = np.arctan(yp/xp)
    psip = orbital.angular(thetap, phip) * orbital.radial(rp)
    anglep = np.angle(psip)
    densityp = (np.conj(psip)*psip).real   
    xpoints = []
    ypoints= []
    zpoints = []
    psipoints=[]
    max_density = np.max(densityp) * 20
    xp = np.ndarray.flatten(xp)
    yp = np.ndarray.flatten(yp)
    zp = np.ndarray.flatten(zp)
    for i, this_x in enumerate(xp):
        for j, this_y in enumerate(yp):
            for k, this_z in enumerate(zp):
                this_density = densityp[i][j][k]
                this_angle = - psip[i][j][k].real
                if np.random.rand()*max_density < this_density:
                    xpoints.append(this_x)
                    ypoints.append(this_y)
                    zpoints.append(this_z)
                    psipoints.append(this_angle)
    xpoints = np.array(xpoints)
    ypoints = np.array(ypoints)
    zpoints = np.array(zpoints)
    psipoints = - np.array(psipoints)  # take the negative to match with the phase colormap colors 
    print(len(xpoints), len(ypoints), len(zpoints), len(psipoints))
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    mlab.points3d(xpoints, ypoints, zpoints, -psipoints, colormap="jet", scale_mode='none', scale_factor=0.2)
    mlab.colorbar(title='Psi', orientation='vertical', nb_labels=5, label_fmt='%.2f')
    mlab.view(-10, 90)
    mlab.show()
if False: # iso + points
    #-------------------------------------------------------------------------------------
    #iso
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    src = mlab.pipeline.scalar_field(x2, y2, z2, np.abs(Phi2))
    src.image_data.point_data.add_array(phase2.T.ravel())
    src.image_data.point_data.get_array(1).name = 'angle'
    src.update()
    src2 = mlab.pipeline.set_active_attribute(src)
    contour = mlab.pipeline.contour(src2)
    contour.filter.contours= [0.00005]
    contour2 = mlab.pipeline.set_active_attribute(contour, point_scalars='angle')    
    mlab.pipeline.surface(contour2, colormap='hsv', vmax=np.pi, vmin=-np.pi)
    mlab.points3d(xpoints, ypoints, zpoints, -psipoints, colormap="jet", scale_mode='none', scale_factor=0.2)
    mlab.view(-10, 90)
    mlab.show()
