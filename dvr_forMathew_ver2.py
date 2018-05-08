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

#This is the method for calculating the rotational energy levels, without adding V(n)
def R(i,j):
	r = np.zeros(x.size)
	rotE = (((1.05457180013904E-034)**2)*(j*(j+1))/(2*8.36998689267015E-028*((i)**2)))/(4.35974465E-18)
	r = r+rotE
	return r

#This method is used to plot V+T+R
def plotVR(eigenVects, eigenVals, Rsize,i):
	#setting up plot	
	fig, ax1 = plt.subplots()
	ax2=ax1.twinx()
	ax2.set_ylabel('Potential Energy')
	ax1.set_ylabel('Eigenvalue')
	ax1.set_ylim([-.3,.4])
	ax2.set_ylim([-.3,.4])
	ax2.set_autoscaley_on(False)
	ax2.plot(x,np.diag(V))	
	plt.title('T+V+R')
	for n in range(0,i):
		print 'n = ',n
		for j in range(0,Rsize):
			#i in R(i,j) is in m
			#r = R(0.756005*0.0000000001, j) #harmonic calc2
			r = R(0.741400*0.0000000001, j) #harmonic theoretical
			#r = R(0.741440*0.0000000001, j) #harmonic experimental
			#r = R(0.754331*0.0000000001, j) #morse calculated
			#r = R(0.741400*0.0000000001, j) #morse theoretical
			#r = R(0.754331*0.0000000001, j) #morse D2
			print 'j = ', j
			print 'r', r[0:20]*219474.631371725
			print 'V', eigenVals[0:20]*219474.631371725
			print 'V+r', (eigenVals[0:20]+r[0])*219474.631371725
			print ''
			ax1.plot(x, 5*(eigenVals[n]+r))	#multiplying by five to see j spacing
	plt.show()

#This method is used to plot V+T
def plotV(eigenVects, eigenVals):
	#setting up plot	
	fig, ax1 = plt.subplots()
	ax2=ax1.twinx()
	ax2.set_ylabel('Potential Energy')
	ax1.set_ylabel('Eigenvalue')
	ax1.set_ylim([-.3,.4])
	ax2.set_ylim([-.3,.4])
	ax2.set_autoscaley_on(False)
	ax2.plot(x,np.diag(V))
	plt.title('V+T')
	for i in range(3):
		ax1.plot(x,((0.05*eigenVects[:,i])+eigenVals[i]))
	plt.show()

def calculateV(x,potential,center=0.0,omega=1.0,fileName=None):
    v=np.zeros(x.size)
    if potential=='harmonic':
        k=1.0
        v=0.5*k*(x*x)
    elif potential == 'harmonic potential':
	#calc k	
		#k = 536.4273415 #(kg/s^2)
		#B = -1.1534912077 #(Eh)
		#i = 0.7462*0.0000000001 #(m)
	#calc2 k, this is the more accurate calculation using a harmonic oscillator
		#k = 544.5495458 #(kg/s^2)
		#B = -1.1534912077 #(Eh)
		#i = 0.756005*0.0000000001 #(m)
	#theoretical k
		k = 575.8823348 #(kg/s^2)
		i = 0.7414*0.0000000001 #(m)
	#experimental k
		#k = 575.267849 #(kg/s^2)
		#i = 0.74144*0.0000000001 #(m)
		v = (0.5*k*(x*0.000000000052917721067-i)**2)/(4.35974465E-18) #(Eh)
	#more accurate model than H, especially at big (x-i) and big n.    
    elif potential == 'morse potential':
	# calc: 
		#be = 1.98841 #(1/angstrom)
		#D = 0.160904 #(Eh)
		#i = 0.754331 #(angstrom)
		#B = -1.1534912077 #(Eh)
	# theoretical:
		#be = 1.93 #(1/angstrom)
		#D = 0.1745515073 #( Eh)
		#i = 0.7414 #(angstrom) ##resave theoretical
	# D2 calc:
		be = 1.98841 #(1/angstrom)
		D = 0.160904 #(Eh)
		i = 0.754331 #(angstrom)
		B = -1.1534912077 #(Eh)
		#while vrot > 0
		v = D*(1 - np.exp(-be*(x*0.5291772107-i)))**2 #(Eh)
    elif potential =='double-well':
        omega=1.0
        k=omega**2
        b=4.0
        a=.50
        g=b**2/(4.0*a)
        v=a*x**4-b*x**2+g
    elif potential =='box':
        v=np.zeros(x.size)
    elif potential =='half-harmonic':
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
    
    v=np.diag(v)
    return v

print 'start'

hbar=1.0
m= 918.8306976 #me not amu
a=-6.
b=6.
N=300
#x=np.linspace(a,b,N)

x=[a+(b-a)/N*float(i) for i in range(1,int(N-1))]
x=np.array(x[1:])

print x[0:10]
print x[x.size-10:x.size]
print x[-1]
print x.size

dx=(b-a)/(N)
print 'N', N
print 'dx', dx
print 'midpoints?', (N/2)-4, x[(N/2)-4:(N/2)+4]
print 'x[3]-x[2], dx',x[3]-x[2], dx
#i=3
T=calculateT(dx,N-3,'infinite plane')
#V=calculateV()
#print 'T'
#plt.imshow(T)
#plt.colorbar()
#plt.show()

k=1.0

V = calculateV(x, 'harmonic potential')
#V = V(x, 'morse potential')
#V=V(x,'box')
#V=V(x,'harmonic')
#V=V(x,'half-harmonic')
#V=V(x,'double-well')
#V=V(x,'scanData',fileName='RvsE.data')

print 'shape of V', V.shape
print 'plot the V'
#plt.imshow(V)
#plt.show()#deltaE(1-0) 0.830601287168

#plt.plot(x,np.diag(V))
#plt.show()

#This is the rotational energy (j transitions)
Rset = np.zeros(x.size)
Rsize = 5

H=T+V
eigenVals,eigenVects=np.linalg.eigh(H)

print 'H', H[0:3,0:3]
print 'plot of H'
#plt.imshow(H)
#plt.colorbar()
#plt.show()

print 'E in cm-1', eigenVals[0:20]*219474.6314

for o in range(1,6):
	print 'deltaE(', o, '-', o-1, ') =', (eigenVals[o]-eigenVals[o-1])*219474.6314, '(1/cm),'

print 'plot of E'
#plt.plot(range(eigenVals.size)[0:25],eigenVals[0:25])
#plt.show()
#plt.plot(range(eigenVals.size),eigenVals)
#plt.show()
print 'Higher Es', eigenVals[240:260]

print 'some of the vects',eigenVects[0,240:260]
print eigenVects[0].size
print x.size
#ExactE=[0.5+1.0*float(i) for i in range(0,int(N-3))]

#print 'ExactE',ExactE
#plt.plot(ExactE[0:20]-eigenVals[0:20])
#plt.show()
#plt.plot(x,np.diag(V))

#plot V+T
plotV(eigenVects,eigenVals)

#This method plots V+T+R for n = 0 to i
plotVR(eigenVects,eigenVals, Rsize,3)

print 'All done!'
