import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
massH=1.00782503223
massD=2.0141017778
massO=15.99491561957
massC=12.00000000000
massF=18.9984032
massConversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28#1822.88839                                           
ang2bohr=1.88973
bohr2ang=1.000/ang2bohr
rad2deg=180.000/np.pi
au2wn=219474.631371725


class dvr (object):

    ###
    ###  This DVR calculator can be used on 1 dimensional systems involving H motion or Hydronium/Ether motion
    ###  It requires either:
    ###          1) an analytical potential energy function.  
    ###       or 2) a file of scan data with the first column being a coordinate in Angstrom or bohr and the 
    ###             second column being in Hartree or wavenumbers (cm^-1).
    ###  General usage:
    ###       (default is Harmonic oscillator approximation for H2)
    ###       In a new python file: 
    ###           Import the dvrClass:                      import dvrClass
    ###           Import Numpy                              import numpy as np
    ###           Make a DVR calculator:                    calculator=dvrClass.dvr()
    ###           Calculate the Kinetic Energy              T=calculator.calculateT()
    ###           Calculate Potential Energy                V=calculator.calculateV('harmonic potential')
    ###           Calculate the Hamiltonian                 H=V+T
    ###           Calculate the energies and wavefunctions  E,Psi_n=calculator.calculateWavefunctions(T+V)
    ###           Save the data                             np.savetxt('Energies.data',E)
    ###                                                     np.savetxt('GroundStateVibrationalWavefunction.data',Psi_n[0]
    ###                                                     np.savetxt('FirstExcitedStateVibrationalWavefunction.data',Psi_n[1])
    ###    
    ###         
    ###   Other notes:      
    ###       E is a list of energies E[0] is the ZPE.  Psi_n is a list of the wavefunctions Psi_n[0] is the first wavefunction
    ###       I highly recommend reading about each function and what the optional inputs are.  You will likely have to specify    
    ###       additional parameter for systems that are not H2.

    def __init__(self,leftX=-6.0,rightX=6.0,numPoints=300,particle_in_motion='H2'):
        #leftX and rightX make up the lowest(left) and highest (right) x values on the DVR grid.  The DVR grid is in bohr!
        #numPoints is related to the number of points on the DVR grid.  The actual number is numPoints-3 (need odd number of points)
        #particle_in_motion can be either 'H2', 'H', or 'Ether-Hydronium'

        self.hbar=1.0 #because atomic units.

        self.m = self.set_mass(particle_in_motion) #in atomic units

        self.a=leftX 
        self.b=rightX
        self.N=numPoints

        
        self.x=[self.a+(self.b-self.a)/self.N*float(i) for i in range(1,int(self.N-1))]
        self.x=np.array(self.x[1:])
        print 'sizes', self.x.size, self.N
        self.dx=(self.b-self.a)/(self.N)
        


######################################################
#### BEGINNING OF ESSENTIAL CALCULATOR FUNCTIONS #####
######################################################

    def calculateWavefunctions(self,H):
        #H= T+V

        eigenVals,eigenVects=np.linalg.eigh(H)
        
        self.wavefunctions=eigenVects.T
        self.EnergyEigenvalues=eigenVals

        return eigenVals,eigenVects.T


    def calculateT(self,interval='infinite plane'):
        #interval can only be 'infinite plane'. Maybe someday I will implement other intervals for the DVR grid.
        T=np.eye(self.x.size)
        if interval=='infinite plane':
            for ii,row in enumerate(T):
                for jj,col in enumerate(row):
                    i=float(ii)
                    j=float(jj)
                    T[ii,jj]=self.hbar*self.hbar/(2.0*self.m*self.dx*self.dx)*(-1.0)**(i-j)
                    if i!=j:
                        T[ii,jj]= T[ii,jj]*2.0/(i-j)**2
                    else:
                        T[ii,jj]=T[ii,jj]* np.pi**2.0/3.0
        self.T=T
        return T

##########################################################
######## BEGINNING OF POTENTIAL SECTION ##################
##########################################################
    def calculateV(self,potential='harmonic potential',center=0.0,omega=1.0,fileName=None,lengthUnit='bohr',energyUnit='Eh',functionalFit='harmonic',plot='False'):
        #potential can either be 'harmonic', 'harmonic potential','morse potential','double-well','box','half-harmonic', or 'scanData'
        #if potential is 'scanData' then it will open the file under "fileName", intrepret it with the "lengthUnit" or "energyUnit" you specify", and try to fit it to the function you specify with "functionalFit". 


        v=np.zeros(self.x.size)
        if potential=='harmonic':
            k=1.0
            v=0.5*k*(self.x*self.x)
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
		v = (0.5*k*(self.x*0.000000000052917721067-i)**2)/(4.35974465E-18) #(Eh)
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
		v = D*(1 - np.exp(-be*(self.x*0.5291772107-i)))**2 #(Eh)
        elif potential =='double-well':
            omega=1.0
            k=omega**2
            b=4.0
            a=.50
            g=b**2/(4.0*a)
            v=a*self.x**4-b*self.x**2+g
        elif potential =='box':
            v=np.zeros(self.x.size)
        elif potential =='half-harmonic':
            k=1.0
            for i,xcoord in enumerate(self.x):
                if xcoord<0:
                    v[i]=1000
                else:
                    v[i]=0.5*k*(xcoord*xcoord)
        elif potential=='scanData':
            print 'fitting a function to the Scan Data from ', fileName
            pot,popt=self.fitPotential(fileName,lengthUnit=lengthUnit,energyUnit=energyUnit,fitTo=functionalFit,plot=plot)
            print 'results from fitPotential', pot, popt
            v=pot(self.x, *popt)

    
        self.V=np.diag(v)
        return self.V
    

    def fitPotential(self,fileName,lengthUnit='bohr',energyUnit='Eh',fitTo='harmonic',plot=False):
        #This is for fitting a predicted potential energy function to some data.  You will fit for the parameter values.
        #Current options for "fitTo" include: harmonic, quartic, morse, lennardjones.  
        #'lengthUnit', and 'energyUnit' refer to the units in 'fileName' that has the first column being the coordinate data,
        #and the second column being the energy data.
        #plot=True will make an xy plot showing you how well the fit matches the data

        rawData=np.loadtxt(fileName)
        xData=rawData[:,0]
        print 'fitting!'
        if lengthUnit=='angstrom':
            xData=xData*ang2bohr
        yData=rawData[:,1]
        if energyUnit=='wavenumbers' or energyUnit=='cm-1' or energyUnit=='cm^-1':
            yData=yData/au2wn

        if fitTo=='harmonic':
            print 'fitting to a Harmonic Function'
            yData=yData-np.min(yData) #Shifted so that the min is a 0
            func=self.Harmonic
            min=np.argmin(yData)
            initialparams=[(yData[min-1]-2.0*yData[min]+yData[min+1])/((xData[1]-xData[0])**2),xData[min],0]

        elif fitTo=='quartic':
            print 'fitting to a Quartic Function'
            yData=yData-np.min(yData)
            initialparams=[1.0,1.0,1.0,3.0,1.0] 
            func=self.Quartic

        elif fitTo=='morse':
            print 'fitting to a Morse Function'
            yData=yData-np.min(yData)
            func=self.Morse
            initialparams=[yData[-1],.7,xData[np.argmin(yData)]]

        elif fitTo=='lennardjones' or fitTo=='LennardJones' or fitTo=='LJ' or fitTo=='6-12':
            print 'fitting to a Lennard-Jones or 6-12 potential'
            func=self.LennardJones
            yData=yData-yData[-1]
            initiaparams=[18.0,xData[np.argmin(yData)]]

        print 'initialParams', initialparams
        if plot : plt.scatter(xData,yData) 
        if plot : plt.plot(xData,func(xData,*initialparams))

        popt, pcov = scipy.optimize.curve_fit(func, xData, yData,p0=initialparams)
        print 'opt values', popt
        xvalues=np.linspace(np.min(xData), np.max(xData),num=50)
        if plot : plt.plot(xvalues, func(xvalues, *popt), 'r-', label='fit') 
        if plot : plt.show() 

        return func,popt

    ###Below here are the functions that can currently be fit with the fitPotential() function

    def Morse(self,x, De, beta,x0):
        return De*(1-np.exp(-beta*(x-x0)))**2

    def Quartic(self, x,A, B, x0,x1,G):
        return A*(x-x0)**4 +B*(x-x1)**2+G
    
    def Harmonic(self,x,A,C,x0):
        return A*(x-x0)**2+C

    def LennardJones(self,x,epsilon, rm):
        return epsilon*(rm/x)**12-2.0*(rm/x)**6

##########################################################
#################END of POTENTIAL SECTION#################
##########################################################

    def set_mass(self,particle):
        #Sets the mass for the particle that is in motion (not fixed in the rigid scan)
        if particle=='H2':
            return 0.5*massH*massConversionFactor# 918.8306976 #in units of me not amu
        elif particle=='H'or particle=='hydrogen' or particle=='proton':
            return 1.0*massH*massConversionFactor# 
        elif particle=='Ether-Hydronium':
            reducedMass= (1.00000000/(4.0*massC+6.0*massH+massO))+(1.00000000/(3.0*massH+massO))
            print 'mass!', 1.0/reducedMass ,'amu', massConversionFactor/reducedMass , 'me'
            return massConversionFactor/reducedMass
        else:
            print '******ERROR!!!******** \n That "particle_in_motion" is not implemented yet! Resorting to default, mass of H',
            print 'but this is mostly likely a bad calculation!!! Please address the issue in the set_mass function'
            print '******ERROR!!!********'
            return 1.0*massH*massConversionFactor
            

##########################################################
##############BEGINNING OF ROTATIONAL SECTION#############
##########################################################

    # This is the method for calculating the rotational energy levels, without adding V(n)
    def R(self,i,j):
	r = np.zeros(self.x.size)
	rotE = (((1.05457180013904E-034)**2)*(j*(j+1))/(2*8.36998689267015E-028*((i)**2)))/(4.35974465E-18)
	r = r+rotE
	return r

##########################################################
#################END of ROTATIONAL SECTION################
##########################################################

##########################################################
######### END OF ESSENTIAL CALCULATOR FUNCTIONS ##########
##########################################################



##########################################################
##############BEGINNING OF PLOTTING SECTION###############
##########################################################

    # This method is used to plot V+T+R
    def plotVR(self,eigenVects, eigenVals, Rsize,i):
	#setting up plot	
	fig, ax1 = plt.subplots()
	ax2=ax1.twinx()
	ax2.set_ylabel('Potential Energy')
	ax1.set_ylabel('Eigenvalue')
	ax1.set_ylim([-.3,.4])
	ax2.set_ylim([-.3,.4])
	ax2.set_autoscaley_on(False)
	ax2.plot(self.x,np.diag(self.V))	
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
			print 'r', r[0:20]*au2wn
			print 'V', eigenVals[0:20]*au2wn
			print 'V+r', (eigenVals[0:20]+r[0])*au2wn
			print ''
			ax1.plot(self.x, 5*(eigenVals[n]+r))	#multiplying by five to see j spacing
	plt.show()
        return


    # This method is used to plot V+T
    def plotV(self,V,eigenVects, eigenVals):
	#setting up plot	
	fig, ax1 = plt.subplots()
	ax2=ax1.twinx()
	ax2.set_ylabel('Potential Energy')
	ax1.set_ylabel('Eigenvalue')
        energies=eigenVals*au2wn

        minE=-energies[0]
        maxE=5.0*(energies[1]-energies[0])

	ax1.set_ylim([minE,maxE])
	ax2.set_ylim([minE,maxE])
	ax2.set_autoscaley_on(False)
	ax2.plot(self.x,np.diag(V)*au2wn)
	plt.title('V+T')

        print 'E0=',energies[0],energies[1]-energies[0], -.0005*au2wn, .006*au2wn
        for i in range(3):
            ax1.plot(self.x,(eigenVects[i]**2*(energies[1]-energies[0])*10.0+energies[i]))
#	for i in range(3):
#            ax1.plot(self.x,((1.0*(eigenVects[i]+eigenVals[i]*au2wn)))
	plt.show()
        return

##########################################################
#################END of PLOTTING SECTION##################
##########################################################

    





