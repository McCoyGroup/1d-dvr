import numpy as np
import matplotlib.pyplot as plt

import dvrClass 


calculator=dvrClass.dvr()
T=calculator.calculateT()

print '###### Harmonic Scan Data ########'
VScan=calculator.calculateV('scanData',fileName='HOScanData.dat')
EScan,Psi_n=calculator.calculateWavefunctions(T+VScan)
print 'E in cm-1', EScan[0:20]*219474.6314

for o in range(1,6):
        print 'deltaE(', o, '-', o-1, ') =', (EScan[o]-EScan[o-1])*219474.6314, '(1/cm),'

calculator.plotV(VScan,Psi_n,EScan)



print '###### Harmonic Function ########'
V=calculator.calculateV('harmonic potential')

E,Psi_n=calculator.calculateWavefunctions(T+V)

print 'E in cm-1', E[0:20]*219474.6314

for o in range(1,6):
        print 'deltaE(', o, '-', o-1, ') =', (E[o]-E[o-1])*219474.6314, '(1/cm),'

calculator.plotV(V,Psi_n,E)

np.savetxt('HOScanData.dat', zip(calculator.x,np.diag(V)))

print '#################'
V2=calculator.calculateV('morse potential')
E2,Psi2_n=calculator.calculateWavefunctions(T+V2)
print 'E in cm-1', E2[0:20]*219474.6314

for o in range(1,6):
        print 'deltaE(', o, '-', o-1, ') =', (E2[o]-E2[o-1])*219474.6314, '(1/cm),'

calculator.plotV(V2,Psi2_n,E2)


