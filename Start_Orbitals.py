# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 13:38:01 2015

# Create the GUI and start the application

@author: Matt B Rowley
"""
from pyface.qt import QtGui, QtCore
import Orbitals_UI
import os
os.environ['ETS_TOOLKIT'] = 'qt4'

def main():
    window = Orbitals_UI.MainWindow()
    window.setAttribute(QtCore.Qt.WA_DeleteOnClose, True)
    window.showMaximized()
    return window

app = QtGui.QApplication.instance()
window = main()
app.exec_()
