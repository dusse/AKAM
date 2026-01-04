import numpy as np

class Initializer:

    
    def __init__(self):
        #   0: run from scratch
        #   1: restart, set input file
        self.runType = 0
        self.inputFile = "./output/_restart.h5"
        # BOX
        self.dim=3
        
        dx=1.0
        dy=3.0
        dz=3.0
        self.boxSize    = [70.0,90.,90.]# in d0 - ion inertia length
        self.boxSizePxl = [int(self.boxSize[0]/dx), int(self.boxSize[1]/dy), int(self.boxSize[2]/dz) ]

        self.pxSize  = self.boxSize[0]/self.boxSizePxl[0]
        self.pxSizeY = self.boxSize[1]/self.boxSizePxl[1]
        
        self.bcType = [1,1,1] #1 - periodic 0 - ideal
        
        self.partclBcType =[1,1,1] #2 - reflect 1 - periodic 0 - outflow
        
        self.dampingBoundaryWidth = [[0,0], [0,0], [0,0]]
         
        self.mpiCores  = [1,2,1]
        
        # time
        self.ts = 0.05
        self.maxtsnum = 401
        self.outputStride = 40

        # output. need to create it before
        self.outputDir = "./output/"
        self.fileTemplate = "test_"
    
        # Particles
        self.ppc4load = 100
        self.ppc = [1, 100]
        self.ppcMinDens = 1.0

        self.numOfSpecies = 2
        self.masses  = [1, 1]
        self.charges = [1, 1]
        self.dens = [0.001, 1.]
        self.vel1 = 0.0
        self.vfl1 = [0.0,0.0,0.0]
        self.vel2 = 0.0
        self.vfl2 = [0.0,0.0,0.0]
       
        self.Pele0 = 0.0
        # Laser spots
        self.spotsNum = 1

        self.criticalDensity = 100.0
        self.potentialThreshold4Interaction = 1e-2
        self.eiCollFreq = 1.

        self.pulseDuration = 150
        self.type2Load = 2
        self.thermVion = 0.01
	
        self.vfl2Load = [0.0,0.0,0.0]	

        self.dens2sustain = self.dens[1]
        self.temp2sustain = 1e-10
    
        #magnetic field magnitude
        self.Bfield = [0.0, 0.0, 0.0]
        
        #ohm's law
        self.resistivity = 0.0

        # pressure tensor
        self.emass = 0.1
        self.tau   = 0.01 # ~ time
        self.relaxFactor = 0.0
        self.smoothStride = int(4)

        #misc
        self.ZERO = 0.0
        self.PI = 3.14159265358979323846
        self.FROZEN = 1
    
    #   0: run from scratch
    #   1: restart, set input file
    def getRunType(self):
        return self.runType
    
    def getInputFile(self):
        return self.inputFile
    
    #   spatial: left lower corner is (0,0,0)
    #   spatial: total box length in normalized units
    def getXright(self):
        return self.boxSize[0]
    
    def getYright(self):
        return self.boxSize[1]
    
    def getZright(self):
        return self.boxSize[2]
    
    
    # if number of pixels is less than 2, there is no direvative
    #   total box length in pixels Lx
    def getXresolution(self):
        return self.boxSizePxl[0]
    #   total box length in pixels Ly
    def getYresolution(self):
        return self.boxSizePxl[1]
    #   total box length in pixels Lz
    def getZresolution(self):
        return self.boxSizePxl[2]
    
    def getXmpiDomainNum(self):
        return self.mpiCores[0]
    
    def getYmpiDomainNum(self):
        return self.mpiCores[1]
    
    def getZmpiDomainNum(self):
        return self.mpiCores[2]
    
    def getElectronPressureXX(self, x, y, z):
        return 1e-8

    def getElectronPressureYY(self, x, y, z):
        return self.getElectronPressureXX(x, y, z)
    
    def getElectronPressureZZ(self, x, y, z):
        return self.getElectronPressureXX(x, y, z)

    
    #   BC
    def getFieldBCTypeX(self):
        return self.bcType[0]
    
    def getFieldBCTypeY(self):
        return self.bcType[1]
    
    def getFieldBCTypeZ(self):
        return self.bcType[2]

    def getDampingBoundaryWidthXleft(self):
        return self.dampingBoundaryWidth[0][0]

    def getDampingBoundaryWidthXright(self):
        return self.dampingBoundaryWidth[0][1]

    def getDampingBoundaryWidthYleft(self):
        return self.dampingBoundaryWidth[1][0]

    def getDampingBoundaryWidthYright(self):
        return self.dampingBoundaryWidth[1][1]

    def getDampingBoundaryWidthZleft(self):
        return self.dampingBoundaryWidth[2][0]

    def getDampingBoundaryWidthZright(self):
        return self.dampingBoundaryWidth[2][1]
    
    def getParticleBCTypeX(self):
        return self.partclBcType[0]
    
    def getParticleBCTypeY(self):
        return self.partclBcType[1]
    
    def getParticleBCTypeZ(self):
        return self.partclBcType[2]
    
    
    #   time
    def getTimestep(self):
        return self.ts
    
    def getMaxTimestepsNum(self):
        return self.maxtsnum
    
    #   output
    def getOutputDir(self):
        return self.outputDir
    
    def getOutputFilenameTemplate(self):
        return self.fileTemplate
    
    def getOutputTimestep(self):
        return self.outputStride

    def getElectronPressure(self, x, y, z):
        return self.Pele0
    

    #   physics: particles
#   first set number of used species
    def getNumOfSpecies(self):
        return self.numOfSpecies
    
    def getMinimumDens2ResolvePPC(self):
        return self.ppcMinDens
    
    #           species 1
    def getPPC4species1(self):
        return self.ppc[0]


    def getMass4species1(self):
        return self.masses[0]
    
    def getCharge4species1(self):
        return self.charges[0]

    def getDensity4species1(self, x, y, z):
        return self.dens[0]

    def getIfParticleTypeIsFrozen4species1(self):
        return 1


    # species 1: modulus of velocity for Maxwell distribution
    def getVelocityX4species1(self, x, y, z):
        return self.vel1
    
    def getVelocityY4species1(self, x, y, z):
        return self.vel1
    
    def getVelocityZ4species1(self, x, y, z):
        return self.vel1


    # species 1: modulus of fluid velocity
    def getFluidVelocityX4species1(self, x, y, z):
        return self.vfl1[0]
    
    def getFluidVelocityY4species1(self, x, y, z):
        return self.vfl1[1]
    
    def getFluidVelocityZ4species1(self, x, y, z):
        return self.vfl1[2]


    #           species 2
    def getPPC4species2(self):
        return self.ppc[1]
    
    def getMass4species2(self):
        return self.masses[1]
    
    def getCharge4species2(self):
        return self.charges[1]
        
    def getDensity4species2(self, x, y, z):
        if (x > 0.3*self.boxSize[0] and x < 0.6*self.boxSize[0]):    
            return self.dens[1]
        else:
            return self.ZERO

    def getIfParticleTypeIsFrozen4species2(self):
        return 0


    # species 2: modulus of velocity for Maxwell distribution
    def getVelocityX4species2(self, x, y, z):
        return self.vel2
    
    def getVelocityY4species2(self, x, y, z):
        return self.vel2
    
    def getVelocityZ4species2(self, x, y, z):
        return self.vel2
    
    # species 2: modulus of fluid velocity
    def getFluidVelocityX4species2(self, x, y, z):
        return self.ZERO
    
    def getFluidVelocityY4species2(self, x, y, z):
        return self.vfl2[1]
    
    def getFluidVelocityZ4species2(self, x, y, z):
        return self.vfl2[2]

    #   physics: laser
    def getNumberOfLaserSpots(self):
        return self.spotsNum
    
    #   physics: laser pulse duration in timesteps
    def getLaserPulseDuration(self):
        return self.pulseDuration
    
    # particles type (1,2,...) is needed to copy mass and charge for injected particles
    def getParticleType2Load(self):
        return self.type2Load

    # can specify ppc for the injected fraction
    def getPPC4loadedParticles(self):
        return self.ppc4load
    
    def getDFtype4InjectedParticles(self):
        return 0

    # injected particles: modulus of fluid velocity
    def getFluidVelocityX4InjectedParticles(self, x, y, z):
        return self.vfl2Load[0]
    
    def getFluidVelocityY4InjectedParticles(self, x, y, z):
        return self.vfl2Load[1]
    
    def getFluidVelocityZ4InjectedParticles(self, x, y, z):
        return self.vfl2Load[2];

    # injected particles: modulus of thermal velocity
    def getVelocityX4InjectedParticles(self, x, y, z):
        return self.thermVion
    
    def getVelocityY4InjectedParticles(self, x, y, z):
        return self.getVelocityX4InjectedParticles(x, y, z)
    
    def getVelocityZ4InjectedParticles(self, x, y, z):
        return self.getVelocityX4InjectedParticles(x, y, z)

    # set density profile to sustain by ablation operator
    def getTargetIonDensity2sustain(self, x, y, z):
        return self.ZERO

    # smooth polynom by roch smets 2014 PoP
    def polynomByRochSmets(self, x): # x = |x|
        return -6.0*x**5+15.0*x**4-10.0*x**3+1.0

    def getIfWeUseIonization(self):
        return 1

    def getIfWeUseElectronIonCollisions(self):
        return 1

    def getElectronIonCollisionFrequency(self):
        return self.eiCollFreq

    def getLaserBeamEnvelopeProfile(self, x, y, z, time):
        omega = 10.
        waist = 12.
        laser_fwhm = 11.0
        sigma = (0.5*laser_fwhm)**2/0.3
        a0 = 1.0
        center = 0.25*self.boxSize[0]
        focus = [ 0.45*self.boxSize[0], 0.5*self.boxSize[1], 0.5*self.boxSize[2] ]
        Zr = omega * waist**2/2.0
        w  = (1./(1.+   ( (x-focus[0])/Zr  )**2 ) )**0.5
        coeff = omega * ( x-focus[0] ) * w**2 / (2.*Zr**2)
        phase = coeff * ( (y-focus[1])**2 + (z-focus[2])**2 )
        exponential_with_total_phase = np.exp(1j*(phase-np.arctan((x-focus[0])/Zr )))
        invWaist2 = (w/waist)**2
        spatial_amplitude = a0* w * np.exp( -invWaist2*(  (y-focus[1])**2 + (z-focus[2])**2 )  )
        time_envelope = np.exp( -( time-center )**2 /sigma )
        space_time_envelope = spatial_amplitude * time_envelope
        return space_time_envelope *  exponential_with_total_phase


    def getPotentialThreshold4Interaction(self):
        return self.potentialThreshold4Interaction

    def getCriticalDensity(self):
        return self.criticalDensity

     # set electron pressure profile to sustain by ablation operator
    def getElectronPressure2sustain(self, x, y, z):
        return self.ZERO

    #   physics: magnetic field
    def getBfieldX(self, x, y, z):
        return self.Bfield[0]
    
    def getBfieldY(self, x, y, z):
        return self.Bfield[1]
    
    def getBfieldZ(self, x, y, z):
        return self.Bfield[2]
    
    #   physics: ohm's law resistivity
    def getResistivity(self):
        return self.resistivity
    
    #   physics: pressure evolution
    def getElectronMass(self):
        return self.emass

    #   physics: pressure evolution : isotropization
    def getRelaxFactor(self):
        return self.relaxFactor

    #   physics: pressure evolution : smoothing
    def getElectronPressureSmoothingStride(self):
        return self.smoothStride

    #   physics: pressure evolution : 1 - isothermal (optional), 0 - evolution equation
    def getIfWeUseIsothermalClosure(self):
        return 0
    #   physics: isothemal closure : electron temperature
    def getElectronTemperature4IsothermalClosure(self):
        return 0.0

    def getDefaultColoumbLogarithm(self):
        return 10

    def getIonIonCollisionFrequencyFactor(self):
        return 0.0

    def getCellBreakdownEfieldFactor(self):
        return 1e10

    def getCriticalPressure(self):
        return 1e10