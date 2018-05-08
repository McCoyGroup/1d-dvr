import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def calculateT(dx,N,interval):

    T=np.eye(N)
    if interval=='infinite plane':
        for ii,row in enumerate(T):
            for jj,col in enumerate(row):
                i=float(ii)
                j=float(jj)
                T[ii,jj]=hbar*hbar/(2.0*m*dx*dx)*(-1.0)**(i-j)
                if i!=j:
                    T[ii,jj]= T[ii,jj]*2.0/(i-j)**2
                else:
                    T[ii,jj]=T[ii,jj]* np.pi**2.0/3.0
    return T

def Morse(x, De, beta,x0):
    return De*(1-np.exp(-beta*x-x0))**2

def Quartic(x, A, B, x0,x1,G):
    return A*(x-x0)**4 +B*(x-x1)**2+G

def fitPotential(fileName):
    rawData=np.loadtxt(fileName)
    xData=rawData[:,0]
    yData=rawData[:,1]
    plt.scatter(xData,yData)
    func=Quartic
    popt, pcov = scipy.optimize.curve_fit(func, xData, yData)
    print 'opt values', popt
    xvalues=np.linspace(np.min(xData), np.max(xData),num=50)
    plt.plot(xvalues, func(xvalues, *popt), 'r-', label='fit')
    plt.show()
    #    print rawData.shape
    return func,popt

def V(x,potential,center=0.0,omega=1.0,fileName=None):

    v=np.zeros(x.size)
    if potential=='morse':
        D_e=3
        alpha=2
        v=D_e*(1-np.exp(alpha*x))**2

    if potential=='harmonic':
        k=1.0
        v=0.5*k*(x*x)
    elif potential=='double-well':
        omega=1.0
        k=omega**2
        b=4.0
        a=.50
        g=b**2/(4.0*a)
        v=a*x**4-b*x**2+g
    elif potential=='box':
        v=np.zeros(x.size)
    elif potential=='half-harmonic':
        k=1.0
        for i,guy in enumerate(x):
            if guy<0:
                v[i]=1000
            else:
                v[i]=0.5*k*(guy*guy)
    elif potential=='scanData':
        print 'doing ScanData, ', fileName
        pot=fitPotential(fileName)
        v,popt=pot(x)

    elif potential=='triple-well':
        print 'triple well!'
        a=np.shape(x)[0]
        b=10.0
        v=x*0
        v[16*a/50:17*a/50]=b
        v[33*a/50:34*a/50]=b
    
    v=np.diag(v)
    return v


print 'start'

hbar=1.0
m=1.0
a=-10.
b=10.
N=503
#x=np.linspace(a,b,N)

x=[a+(b-a)/N*float(i) for i in range(1,int(N-1))]
x=np.array(x[1:])
dx=(b-a)/(N)
print 'N', N
print 'dx', dx

T=calculateT(dx,N-3,'infinite plane')

#print 'T'
#plt.imshow(T)
#plt.colorbar()
#plt.show()

k=1.0

V=V(x,'box')
#V=V(x,'harmonic')
#V=V(x,'morse')
#V=V(x,'triple-well')
#V=V(x,'half-harmonic')
#V=V(x,'double-well')
#V=V(x,'scanData',fileName='RvsE.data')

#print 'shape of V', V.shape
#print 'plot the V'
#plt.imshow(V)
#plt.show()
plt.plot(x,np.diag(V))
plt.show()

H=T+V #add them
#print 'plot of H'
#plt.imshow(H)
#plt.colorbar()
#plt.show()

eigenVals,eigenVects=np.linalg.eigh(H) #diagonalize it

#print 'H', H[0:3,0:3]
print 'E', eigenVals[0:20]
print 'plot of E'
#plt.plot(range(eigenVals.size)[0:25],eigenVals[0:25])
#plt.show()
#plt.plot(range(eigenVals.size),eigenVals)
#plt.show()
#print 'Higher Es', eigenVals[240:260]
#print 'some of the vects',eigenVects[0,240:260]
#print eigenVects[0].size
#print x.size
#ExactE=[0.5+1.0*float(i) for i in range(0,int(N-3))]
#print 'ExactE',ExactE
#plt.plot(ExactE[0:20]-eigenVals[0:20])
#plt.show()
#plt.plot(x,np.diag(V))

fig, ax1 = plt.subplots()
ax2=ax1.twinx()
ax2.set_ylabel('Potential Energy')
ax1.set_ylabel('Eigenvalue')
plt.ylim([0,10])
ax2.set_autoscaley_on(False)
ax2.plot(x,np.diag(V))

for i in range(3):
    ax1.plot(x,(10.0*(eigenVects[:,i])+eigenVals[i]))
    
plt.show()

print 'All done!'
