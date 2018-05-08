import numpy as np
import matplotlib.pyplot as plt

import dvrClass 


calculator=dvrClass.dvr(leftX=3.8, rightX=6.0,particle_in_motion='Ether-Hydronium')

T=calculator.calculateT()

VScan=calculator.calculateV('scanData',fileName='OOScanData.dat',lengthUnit='angstrom',energyUnit='Eh',functionalFit='morse',plot=True)

E,Psi_n=calculator.calculateWavefunctions(T+VScan)

print 'E in cm-1', E[0:20]*219474.6314

for o in range(1,6):
    print 'deltaE(', o, '-', o-1, ') =', (E[o]-E[o-1])*219474.6314, '(1/cm),'

calculator.plotV(VScan,Psi_n,E)


