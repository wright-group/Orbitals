# -*- coding: utf-8 -*-
# This program is licences under an MIT license. Full licence is at the end of
# this file.
"""
Orbitals_UI.py
A user interface for generating and interacting with various animations of
hydrogenic wavefunctions.

@author: Matthew B Rowley
"""

# A few global variables

import os
os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui, QtCore
from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor
import Hydrogenic as hyd
import time
import numpy as np
#os.environ['ETS_TOOLKIT'] = 'qt4'

sqrt2 = np.sqrt(2)

class MainWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.setWindowTitle("Orbitals!")
        self.frame = QtGui.QFrame(self)
        self.layout = QtGui.QHBoxLayout(self.frame)
        self.left_layout = QtGui.QVBoxLayout()
        self.tabs = QtGui.QTabWidget(self, )
        self.tabs.addTab(StationaryPanel(self), 'Stationary States')
        self.tabs.addTab(CoherencePanel(self), 'Coherences')
        self.tabs.addTab(CrossingPanel(self), 'Avoided Crossing')
        self.tabs.setMaximumWidth(550)
        self.tabs.currentChanged.connect(self.changeTab)
        self.left_layout.addWidget(self.tabs)
        self.author = AuthorPanel(self)
        self.left_layout.addWidget(self.author)
        self.layout.addLayout(self.left_layout)
        self.mayavi_widget = MayaviQWidget(self.frame)
        self.layout.addWidget(self.mayavi_widget)
        self.setCentralWidget(self.frame)

    def closeEvent(self, evt):
        global calculator
        QtGui.QMainWindow.closeEvent(self, evt)
        calculator.animation_timer.stop()

    def changeTab(self):
        global calculator
        calculator.writeMode(self.tabs.tabText(self.tabs.currentIndex()))


class OrbitalsPanel(QtGui.QWidget):
    def __init__(self, my_orbital, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.my_orbital = my_orbital
        self.n = 0
        self.l = 0
        self.m = 0
        self.s = 1
        self.layout = QtGui.QVBoxLayout(self)
        self.invert = QtGui.QCheckBox(self, text='Invert Spin')
        self.layout.addWidget(self.invert)
        self.orbitals_layout = QtGui.QHBoxLayout()
        self.real = QtGui.QVBoxLayout()
        self.real_label = QtGui.QLabel(self, text='Real')
        self.real.addWidget(self.real_label)
        self.s1 = OrbitalButton(self.my_orbital, self, text='1s')
        self.real.addWidget(self.s1)
        self.s2 = OrbitalButton(self.my_orbital, self, text='2s')
        self.real.addWidget(self.s2)
        self.p2_layout = QtGui.QHBoxLayout()
        self.p2x = OrbitalButton(self.my_orbital, self, text='2px')
        self.p2_layout.addWidget(self.p2x)
        self.p2y = OrbitalButton(self.my_orbital, self, text='2py')
        self.p2_layout.addWidget(self.p2y)
        self.p2z = OrbitalButton(self.my_orbital, self, text='2pz')
        self.p2_layout.addWidget(self.p2z)
        self.real.addLayout(self.p2_layout)
        self.real.addItem(VerticalSpacer())
        self.layout.addLayout(self.real)
        self.layout.addWidget(HorizontalLine(self))
        self.numbers_layout = QtGui.QHBoxLayout()
        self.num_label = QtGui.QLabel(self,
                                      text='Set Quantum Numbers Manually')
        self.layout.addWidget(self.num_label)
        self.n_label = QtGui.QLabel(self, text='n')
        self.numbers_layout.addWidget(self.n_label)
        self.num_n = QtGui.QSpinBox(self)
        self.num_n.setMinimum(1)
        self.num_n.valueChanged.connect(self.testNumbers)
        self.numbers_layout.addWidget(self.num_n)
        self.numbers_layout.addItem(HorizontalSpacer())
        self.numbers_layout.addWidget(VerticalLine(self))
        self.numbers_layout.addItem(HorizontalSpacer())
        self.l_label = QtGui.QLabel(self, text='l')
        self.numbers_layout.addWidget(self.l_label)
        self.num_l = QtGui.QSpinBox(self)
        self.num_l.valueChanged.connect(self.testNumbers)
        self.numbers_layout.addWidget(self.num_l)
        self.numbers_layout.addItem(HorizontalSpacer())
        self.numbers_layout.addWidget(VerticalLine(self))
        self.numbers_layout.addItem(HorizontalSpacer())
        self.m_label = QtGui.QLabel(self, text='m')
        self.numbers_layout.addWidget(self.m_label)
        self.num_m = QtGui.QSpinBox(self)
        self.num_m.valueChanged.connect(self.testNumbers)
        self.numbers_layout.addWidget(self.num_m)
        self.numbers_layout.addItem(HorizontalSpacer())
        self.numbers_layout.addWidget(VerticalLine(self))
        self.numbers_layout.addItem(HorizontalSpacer())
        self.s_label = QtGui.QLabel(self, text='s')
        self.numbers_layout.addWidget(self.s_label)
        self.num_s = QtGui.QSpinBox(self)
        self.num_s.valueChanged.connect(self.testNumbers)
        self.num_s.setRange(-1, 1)
        self.num_s.setValue(1)
        self.numbers_layout.addWidget(self.num_s)
        self.layout.addLayout(self.numbers_layout)
        self.layout.addItem(VerticalSpacer())

    def chooseNumbers(self):
        global signals
        self.my_orbital.writeNumbers([self.n, self.l, self.m, self.s])
        signals.orbital_change.emit()

    def testNumbers(self):
        self.valid_numbers = True
        self.n = self.num_n.value()
        self.num_l.setMaximum(self.n - 1)
        self.l = self.num_l.value()
        self.num_m.setRange(-self.l, self.l)
        self.m = self.num_m.value()
        if self.num_s.value() == 0:
            self.s = -self.s
            self.num_s.setValue(self.s)
        self.chooseNumbers()


class StationaryPanel(QtGui.QWidget):
    def __init__(self, parent=None):
        global stationary_orbital
        QtGui.QWidget.__init__(self, parent)
        self.layout = QtGui.QVBoxLayout(self)
        self.inst_text = '''Stationary State -
        A stationary state is a population in a single quantum state. Select a
        common named orbital (either real, or complex), or input arbitrary
        quantum numbers.'''
        self.instructions = QtGui.QLabel(self, text=self.inst_text)
        self.layout.addWidget(self.instructions)
        self.layout.addWidget(HorizontalLine(self))
        self.orbitals = OrbitalsPanel(stationary_orbital, self)
        self.layout.addWidget(self.orbitals)


class CoherencePanel(QtGui.QWidget):
    def __init__(self, parent=None):
        global ket_orbital, bra_orbital, signals
        QtGui.QWidget.__init__(self, parent)
        self.layout = QtGui.QVBoxLayout(self)
        self.inst_text = '''Coherence -
        A coherence is an oscillating wavefunction which is composed of two
        states. Usually an external electric field (light) induces a coherence
        between the ground and excited state for the transition which is
        resonant with the light.'''
        self.instructions = QtGui.QLabel(self, text=self.inst_text)
        self.layout.addWidget(self.instructions)
        self.layout.addWidget(HorizontalLine(self))
        self.options = QtGui.QHBoxLayout()
        self.rabi = QtGui.QCheckBox(self, text='Show Rabi Cycle')
        self.options.addWidget(self.rabi)
        self.layout.addLayout(self.options)
        self.functions = QtGui.QHBoxLayout()
        self.ket_layout = QtGui.QVBoxLayout()
        self.ket_label = QtGui.QLabel(self, text='|ket>')
        self.ket_layout.addWidget(self.ket_label)
        self.ket_layout.addWidget(HorizontalLine())
        self.ket_orbitals = OrbitalsPanel(ket_orbital, self)
        self.ket_layout.addWidget(self.ket_orbitals)
        self.bra_layout = QtGui.QVBoxLayout()
        self.bra_label = QtGui.QLabel(self, text='<bra|')
        self.bra_layout.addWidget(self.bra_label)
        self.bra_layout.addWidget(HorizontalLine())
        self.bra_orbitals = OrbitalsPanel(bra_orbital, self)
        self.bra_layout.addWidget(self.bra_orbitals)
        self.functions.addLayout(self.ket_layout)
        self.functions.addWidget(VerticalLine())
        self.functions.addLayout(self.bra_layout)
        self.layout.addLayout(self.functions)
        self.animate_button = QtGui.QPushButton(self,
                                                text='Start/Stop Animation')
        self.layout.addWidget(self.animate_button)
        self.animate_button.clicked.connect(signals.animate_clicked.emit)


class CrossingPanel(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.layout = QtGui.QVBoxLayout(self)
        self.inst_text = '''Avoided Crossing -
        An avoided crossing can illustrate the consequences of the breakdown of
        the Born-Oppenheimer approximation.'''
        self.instructions = QtGui.QLabel(self, text=self.inst_text)
        self.layout.addWidget(self.instructions)
        self.layout.addWidget(HorizontalLine(self))
        self.layout.addItem(VerticalSpacer())


class AuthorPanel(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.layout = QtGui.QVBoxLayout(self)
        self.author = QtGui.QLabel(self,
                                   text='Written by: Matthew Rowley - 2015')
        self.layout.addWidget(self.author)
        self.email_label = QtGui.QLabel(self,
                                        text='email: Matt.B.Rowley@gmail.com')
        self.layout.addWidget(self.email_label)


class VerticalSpacer(QtGui.QSpacerItem):
    def __init__(self):
        QtGui.QSpacerItem.__init__(self, 20, 20, QtGui.QSizePolicy.Minimum,
                                   QtGui.QSizePolicy.Expanding)


class HorizontalSpacer(QtGui.QSpacerItem):
    def __init__(self):
        QtGui.QSpacerItem.__init__(self, 20, 20, QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Minimum)


class VerticalLine(QtGui.QFrame):
    def __init__(self, parent=None):
        QtGui.QFrame.__init__(self, parent)
        self.setFrameShape(QtGui.QFrame.VLine)
        self.setFrameShadow(QtGui.QFrame.Sunken)


class HorizontalLine(QtGui.QFrame):
    def __init__(self, parent=None):
        QtGui.QFrame.__init__(self, parent)
        self.setFrameShape(QtGui.QFrame.HLine)
        self.setFrameShadow(QtGui.QFrame.Sunken)


class OrbitalButton(QtGui.QPushButton):
    def __init__(self, my_orbital, parent=None, text=''):
        QtGui.QPushButton.__init__(self, parent=None, text=text)
        self.clicked.connect(self.writeOrbital)
        self.target_text = text
        self.my_orbital = my_orbital

    def writeOrbital(self):
        global signals
        self.my_orbital.writeName(self.target_text)
        signals.orbital_change.emit()


class Visualization(HasTraits):
    '''
    Mayavi visualization
    '''
    scene = Instance(MlabSceneModel, ())

    @on_trait_change('scene.activated')
    def createPlot(self):
        # This function is called when the view is opened. We don't
        # populate the scene when the view is not yet open, as some
        # VTK features require a GLContext.

        # We can do normal mlab calls on the embedded scene.
        global points
        self.scene.mlab.clf()
        x, y, z, psi = points.readPoints()
        self.mesh = self.scene.mlab.mesh(x, y, z, scalars=psi, colormap='jet', vmax=np.pi, vmin=-np.pi)
        #jet_map = self.mesh.module_manager.scalar_lut_manager.lut.table.to_array()
        #jet_map2 = jet_map[::-1]
        #jet_map2[:,1] = jet_map[:,2]
        #jet_map2[:,2] = jet_map[:,1]
        #jet_map2 = jet_map2[::-1]
        #my_map = np.concatenate((jet_map, jet_map2))
        my_map = [[255,0,0,255],[140,0,140,255],[140,0,140,255],[0,0,255,255],[0,0,255,255],[255,153,18,255],[255,153,18,255],[255,0,0,255]]
        self.mesh.module_manager.scalar_lut_manager.lut.table = my_map
        self.axes = self.scene.mlab.axes()
        self.fig = self.scene.mlab.gcf()
        self.source = self.mesh.mlab_source
        self.scene.mlab.view(0, 90, np.max((x, y, z))*5, (0, 0, 0))
        self.scene.reset_zoom()

    # the layout of the dialog screated
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True  # We need this to resize with the parent widget
                )

    def updatePoints(self):
        global points, zoom
        x, y, z, psi = points.readPoints()
        self.source.set(x=x, y=y, z=z, scalars=psi)
        if zoom.read():
            self.scene.reset_zoom()
        self.axes.remove()
        self.axes = self.scene.mlab.axes()
        #zoom.write(False)


class MayaviQWidget(QtGui.QWidget):
    '''
    Widget containing the visualization
    '''
    def __init__(self, parent=None):
        global signals
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.visualization = Visualization()
        self.ui = self.visualization.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)
        signals.update_visualization.connect(self.visualization.updatePoints)


class OrbitalCalculator(QtCore.QMutex):
    '''
    Calculates the appropriate surface to graph for stationary states
    and coherences. Avoided crossings are handled separately.
    '''
    def __init__(self, stationary_orbital, ket_orbital, bra_orbital):
        global signals
        QtCore.QMutex.__init__(self)
        self.stationary_orbital = stationary_orbital
        self.ket_orbital = ket_orbital
        self.bra_orbital = bra_orbital
        self.mode = 'Stationary States'
        self.rabi = False
        self.animating = False
        '''
        phi = []
        theta = []
        num = 15
        phi_vals = np.linspace(0, 1 * np.pi, num / 2)
        for phival in phi_vals:
            theta_vals = np.linspace(0, 2 * np.pi, np.sin(phival) * num + 1)
            theta_offset = np.random.rand()
            for thetaval in theta_vals:
                phi.append(phival)
                theta.append(thetaval + theta_offset)
        new_phi=[]
        new_theta=[]
        length = len(phi)
        for i, phi_val in enumerate(phi):
            new_theta.append(theta)
            new_phi.append(np.ones_like(phi)*phi_val)
        self.phi = np.array(new_phi)
        self.theta = np.array(new_theta)
        '''
        self.phi, self.theta = np.mgrid[0:np.pi:50j, 0:2*np.pi:100j]
        signals.orbital_change.connect(self.orbitalChange)
        signals.animate_clicked.connect(self.setupAnimation)
        self.orbitalChange()
        self.animation_timer = QtCore.QTimer()
        self.animation_timer.timeout.connect(self.runAnimation)

    def writeMode(self, new_mode):
        self.lock()
        self.mode = new_mode
        self.unlock()

    def orbitalChange(self):
        global signals, points
        if self.mode == 'Stationary States':
            orbital = self.stationary_orbital.read()
            r = orbital.r_90p
            x = r * np.sin(self.phi) * np.cos(self.theta)
            y = r * np.sin(self.phi) * np.sin(self.theta)
            z = r * np.cos(self.phi)
            rs = (np.conj(orbital.angular(self.theta, self.phi)) *
                  orbital.angular(self.theta, self.phi)).real
            psi = np.angle(orbital.angular(self.theta, self.phi))
            data = (rs*x, rs*y, rs*z, psi)
            points.writePoints(data)
            zoom.write(True)
            signals.update_visualization.emit()

    def setupAnimation(self):
        global zoom
        if not self.animating:
            self.ket = self.ket_orbital.read()
            self.bra = self.bra_orbital.read()
            if self.rabi is False:
                self.radial = self.radialNoCycle
                self.angular = self.angularNoCycle
                self.angularConjugate = self.angularNoCycleConjugate
                self.times = np.linspace(0, 2*np.pi, 100)
            else:
                self.radial = self.radialRabiCycle
                self.angular = self.angularRabiCycle
                self.angular = self.angularRabiCycleConjugate
                self.times = np.linspace(0, 2*np.pi, 400)
            self.i = 0
            self.animation_timer.start(25)
            self.animating = True
            zoom.write(False)
        else:
            self.animation_timer.stop()
            self.animating = False

    def runAnimation(self):
        global zoom
        t = self.times[self.i % len(self.times)]
        r = self.radial(t)
        x = r * np.sin(self.phi) * np.cos(self.theta)
        y = r * np.sin(self.phi) * np.sin(self.theta)
        z = r * np.cos(self.phi)
        rs = (self.angularConjugate(self.theta, self.phi, t) *
              self.angular(self.theta, self.phi, t)).real
        psi = np.angle(self.angular(self.theta, self.phi, t))
        data = (rs*x, rs*y, rs*z, psi)
        points.writePoints(data)
        signals.update_visualization.emit()
        self.i = self.i + 1
        if not self.mode == 'Coherences':
            self.animating = False
            self.animation_timer.stop()

    def radialNoCycle(self, t):
        return 0.5/sqrt2*(self.ket.r_90p + self.bra.r_90p)

    def angularNoCycle(self, theta, phi, t):
        return (np.exp(1j*3*t)*self.ket.angular(theta, phi) +
                np.exp(1j*1*t)*self.bra.angular(theta, phi))

    def angularNoCycleConjugate(self, theta, phi, t):
        return (np.exp(-1*1j*3*t)*self.ket.angular(theta, phi) +
                np.exp(-1*1j*1*t)*self.bra.angular(theta, phi))

    def radialRabiCycle(self, t):
        return 0.5/sqrt2*(np.sin(t)*np.sin(t)*self.ket.r_90p +
                    np.cos(t)*np.cos(t)*self.bra.r_90p)

    def angularRabiCycle(self, theta, phi, t):
        return (np.sin(t)*np.exp(1j*5*t)*self.ket.angular(theta, phi) +
                np.cos(t)*np.exp(1j*1*t)*self.bra.angular(theta, phi))

    def angularRabiCycleConjugate(self, theta, phi, t):
        return (np.sin(t)*np.exp(1j*-5*t)*self.ket.angular(theta, phi) +
                np.cos(t)*np.exp(1j*-1*t)*self.bra.angular(theta, phi))


class OrbitalMutex(QtCore.QMutex):
    '''
    Stores a Hydrogenic.Orbital object
    '''
    def __init__(self):
        QtCore.QMutex.__init__(self)
        self.orbital = hyd.orbitals['2pz']

    def read(self):
        return self.orbital

    def writeName(self, new_name):
        self.lock()
        self.orbital = hyd.orbitals[new_name]
        self.unlock()

    def writeNumbers(self, new_numbers):
        self.lock()
        n = new_numbers[0]
        l = new_numbers[1]
        m = new_numbers[2]
        s = new_numbers[3]
        self.orbital = hyd.Orbital(n, l, m, s)
        self.unlock()


class PointsMutex(QtCore.QMutex):
    '''
    Stores the points for visualization
    '''
    def __init__(self):
        QtCore.QMutex.__init__(self)
        self.points = (0, 0, 0, 0)

    def readPoints(self):
        return self.points

    def writePoints(self, new_points):
        self.lock()
        self.points = new_points
        self.unlock()

class BoolMutex(QtCore.QMutex):
    '''
    Stores a boolean variable
    '''
    def __init__(self, init_value=False):
        QtCore.QMutex.__init__(self)
        self.value = init_value

    def read(self):
        return self.value

    def write(self, new_value):
        self.lock()
        self.value = new_value
        self.unlock()

class Signals(QtCore.QObject):
    '''
    A QObject that holds all the pyqtSignals
    '''
    orbital_change = QtCore.pyqtSignal()
    animate_clicked = QtCore.pyqtSignal()
    update_visualization = QtCore.pyqtSignal()
    frame_rendered = QtCore.pyqtSignal()

signals = Signals()
stationary_orbital = OrbitalMutex()
bra_orbital = OrbitalMutex()
ket_orbital = OrbitalMutex()
points = PointsMutex()
zoom = BoolMutex(init_value=True)
calculator = OrbitalCalculator(stationary_orbital, ket_orbital, bra_orbital)

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
