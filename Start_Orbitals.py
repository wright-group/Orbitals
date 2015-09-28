# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 13:38:01 2015
# Instantiate the application


# Create the GUI and start the application

@author: Matt B Rowley
"""
import sys
import PyQt4.QtCore as QtCore
import PyQt4.QtGui as QtGui
import Orbitals_UI


def main():
    window = Orbitals_UI.MainWindow()
    window.setAttribute(QtCore.Qt.WA_DeleteOnClose, True)
    window.showMaximized()
    return window

#app = QtGui.QApplication(sys.argv)
window = main()
#app.exec_()
