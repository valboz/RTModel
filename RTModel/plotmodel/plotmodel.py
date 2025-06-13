import VBMicrolensing
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.animation as animation
from matplotlib.patches import Circle
import os
import numpy as np
import shutil
from PIL import Image
from tqdm import tqdm
import sys
import inspect
import glob
from tabulate import tabulate


class plotmodel:
    def __init__(self, eventname,model = '', tmin = '', tmax = '', magmin = '', magmax = '', tlabel='t', maglabel='mag', reslabel = 'Res',
                 referencephot = 0, timesteps = 300, 
                 modelfile = None, parameters = [], line = 0,printpars = True, animate = False,interval = 1000, 
                 satellitedir = '.', accuracy = 0.01, colors = None, satellitecolors = None):
        self.satellitedir = satellitedir
        self.parameters = parameters
        filin=inspect.getfile(VBMicrolensing)
        self.eventname = eventname
        self.model = model
        self.tmin = tmin
        self.tmax = tmax
        self.inputmagmin = magmin
        self.inputmagmax = magmax
        self.tlabel = tlabel
        self.maglabel = maglabel
        self.reslabel = reslabel
        self.timesteps = timesteps
        self.referencephot = referencephot
        self.modelfile = modelfile
        self.line = line
        self.printpars = printpars
        self.animate = animate
        if(self.model == ''):
            if(modelfile == None):
                print('Please specify a model or a modelfile')
                return
            self.model = os.path.basename(self.modelfile)
        if(modelfile == None and parameters == []):
            print('Please specify model parameters')
            return
        self.vbm = VBMicrolensing.VBMicrolensing()
        self.vbm.Tol = accuracy
        self.vbm.SetMethod(VBMicrolensing.VBMicrolensing.Multipoly)
        # General information on models
        self.modelcodes= ['PS','PX','BS','BO','LS','LX','LO','LK','TS','TX']
        self.npars=[4,6,7,12,7,9,12,14,10,12]
        self.logposs=[[0,1,3],
                 [1,3],
                 [0,1,6],
                 [0,1,6],
                 [0,1,4,5],
                 [0,1,4,5],
                 [0,1,4,5],
                 [0,1,4,5],
                 [0,1,4,5,7,8],
                 [0,1,4,5,7,8]]
        self.parnames = [['u0','tE','t0','rho'],
                    ['u0','tE','t0','rho','piN','piE'],
                    ['tE','FR','u01','u02','t0','t02','rho'],
                    ['tE','FR','u01','u02','t0','t02','rho','piN','piE','gamma1','gamma2','gammaz'],
                    ['s','q','u0','alpha','rho','tE','t0'],
                    ['s','q','u0','alpha','rho','tE','t0','piN','piE'],
                    ['s','q','u0','alpha','rho','tE','t0','piN','piE','gamma1','gamma2','gammaz'],
                    ['s','q','u0','alpha','rho','tE','t0','piN','piE','gamma1','gamma2','gammaz','sz_s','a_s3d'],
                    ['s','q','u0','alpha','rho','tE','t0','s2','q2','beta'],
                    ['s','q','u0','alpha','rho','tE','t0','s2','q2','beta','piN','piE']]
        self.astroparnames = ['muS_Dec','muS_RA','piS','thetaE']
        if(colors != None):
            self.colors = colors
        else:
            self.colors = ['blue','red','green','darkorange','magenta','cyan','gray','teal','maroon','gold','lime','darkviolet']
        self.sourcecolor = 'pink'
        self.causticcolor = (0.5,0.5,0.5)
        if(satellitecolors != None):
            self.satellitecolors = colors
        else:
            self.satellitecolors = ['k-','r-','g-','b-','y-']
        self.legendlocation = 'best'
        self.astrometric = False
        self.nlinpars = 2
        self.parstring = ''
        if(animate):
            self.animateplots(interval = interval)
        else:
            self.readdata()
            self.readoptions()
            self.readparameters()
            self.calculate()
            self.showall()
        
    # Reading data from LCToFit.txt
    def readdata(self):
        if(self.eventname== None):
            self.lightcurves = []
            self.telescopes = []            
            self.limbdarkenings = []   
            self.nfil = 0
            self.npoints = 0
        else:
            os.chdir(self.eventname)
            with open('LCToFit.txt') as f:
                self.npoints=int(f.readline())
                ofil=0
                data=[[]]
                while True:
                    line=f.readline()
                    if(line == ''):
                        break
                    chunks=line.split(' ')
                    self.nfil=int(chunks[0])
                    if(self.nfil != ofil):
                        data.append([])
                        ofil=self.nfil
                    data[self.nfil].append([float(chunks[1]),float(chunks[2]),float(chunks[3]),int(chunks[4]), float(chunks[5]),float(chunks[6]),float(chunks[7]),float(chunks[8])])
                    if(data[-1][-1][-3]>0):
                        self.astrometric = True
                self.nfil +=1
            if(self.astrometric):
                self.nlinpars = 4
                for i in range(len(self.npars)):
                    self.npars[i] += 4
                    self.parnames[i] += self.astroparnames
            self.lightcurves = [ [np.array([dl[0] for dl in d]),np.array([dl[1] for dl in d]),np.array([dl[2] for dl in d]) , np.array([dl[4] for dl in d]),np.array([dl[5] for dl in d]),np.array([dl[6] for dl in d]),np.array([dl[7] for dl in d]),d[0][3]] for d in data]
            with open('FilterToData.txt') as f:
                self.telescopes = f.readlines()
                for i in range(0,self.nfil):
                    self.telescopes[i] = self.telescopes[i][0:self.telescopes[i].index('.')]
                self.telescopes = [tel.split('_')[0] for tel in self.telescopes]
            while(len(self.colors)<self.nfil):
                self.colors.extend(self.colors)
            if(os.path.exists(self.eventname + '/Data/LimbDarkening.txt')):
                with open(self.eventname + '/Data/LimbDarkening.txt') as f:
                    lines = f.readlines()
                    self.limbdarkenings = [float(ld) for ld in lines]
            else:
                self.limbdarkenings = [0 for t in self.telescopes]

    # Reading options from LevMar.ini
    def readoptions(self):
        if(os.path.exists(self.eventname + '/ini/LevMar.ini')):
            with open(self.eventname + '/ini/LevMar.ini') as f:
                lines = f.readlines()
                for line in lines:
                    chunks = line.split()
                    if(chunks[0] == 'turn_off_secondary_source' and chunks[2] == 'True'):
                        self.vbm.turn_off_secondary_source = True
                    elif(chunks[0] == 'turn_off_secondary_lens' and chunks[2] == 'True'):
                        self.vbm.turn_off_secondary_lens = True
                    elif(chunks[0] == 'mass_luminosity_exponent'):
                        self.vbm.mass_luminosity_exponent = float(chunks[2])
                    elif(chunks[0] == 'mass_radius_exponent'):
                        self.vbm.mass_radius_exponent = float(chunks[2])
                    elif(chunks[0] == 'lens_mass_luminosity_exponent'):
                        self.vbm.lens_mass_luminosity_exponent = float(chunks[2])
                        
    # Reading model parameters
    def readparameters(self):
        self.modnumber = self.modelcodes.index(self.model[0:2])
        if(self.eventname != None):
            os.chdir(self.eventname)
        if(self.parameters == []):
            with open(self.modelfile) as f:
                lines = f.readlines()
                if(not lines[0][0].isnumeric()):
                    lines.pop(0)
                self.maxlines = len(lines)
                if(self.line>=self.maxlines):
                    return False
                else:
                    line = lines[self.line]
                    chunks = line.split(' ')
                    values = [float(v) for v in chunks]
                    self.pars = values[0:self.npars[self.modnumber]]
                    self.parsprint=self.pars[:]
                    for i in self.logposs[self.modnumber]:
                        self.pars[i] = math.log(self.pars[i])
                    self.blends = np.array([values[self.npars[self.modnumber]+i*self.nlinpars] for i in range(0,self.nfil)])
                    self.sources = np.array([values[self.npars[self.modnumber]+i*self.nlinpars+1] for i in range(0,self.nfil)])
                    if(self.astrometric):
                        self.L0Dec = np.array([values[self.npars[self.modnumber]+i*self.nlinpars+2] for i in range(0,self.nfil)])
                        self.L0RA = np.array([values[self.npars[self.modnumber]+i*self.nlinpars+3] for i in range(0,self.nfil)])
                    self.blendings = self.blends/(self.sources+1.e-12*self.blends)
                    self.baselines = -2.5*np.log10(self.blends+self.sources)
                    self.chi2=values[-1]
                    #Errors
                    if(not self.animate):
                        line = lines[self.line +1 ]
                        chunks = line.split(' ')
                        values = [float(v) for v in chunks]
                        self.parerrs = values[0:self.npars[self.modnumber]]
                        self.blendingerrs = np.array([values[self.npars[self.modnumber]+i*self.nlinpars] for i in range(0,self.nfil)])
                        self.baselineerrs = np.array([values[self.npars[self.modnumber]+i*self.nlinpars+1] for i in range(0,self.nfil)])
                    return True
        else:
            self.pars = self.parameters[0:self.npars[self.modnumber]]
            self.parsprint=self.pars[:] 
            for i in self.logposs[self.modnumber]:
                self.pars[i] = math.log(self.pars[i])
            self.blends = []
            self.sources = []
            self.L0Dec = []
            self.L0RA = []
            self.chi2 = 0
            for i in range(0,self.nfil):
                lc0 = self.lightcurves[i]
                lcarr = np.array(lc0[0:3])
                lctran=np.transpose(lcarr)
                lc = np.transpose(lctran)                
                self.t = lc[0]
                self.vbm.satellite = lc0[7]
                self.lightcurve()
                f = np.array(self.results[0])
                sumf = (f/(lc[2]*lc[2])).sum()
                sumf2 = (f*f/(lc[2]*lc[2])).sum()
                sumsigma =(1/(lc[2]*lc[2])).sum()
                sumy = (lc[1]/(lc[2]*lc[2])).sum()
                sumfy = (f*lc[1]/(lc[2]*lc[2])).sum()
                sumy2 = (lc[1]*lc[1]/(lc[2]*lc[2])).sum()
                p1=sumf*sumf-sumf2*sumsigma + 1.0e-50;
                self.blends.append((sumf*sumfy-sumf2*sumy)/p1)
                self.sources.append((sumf*sumy-sumsigma*sumfy)/p1)
                self.chi2 += ((lc[1] - self.blends[-1] - self.sources[-1]*f)**2/(lc[2]*lc[2])).sum()
                L0N = L0E = 0
                if(self.astrometric and self.lightcurves[i][4][0]>0):
                    c12 =  self.vbm.CombineCentroids(self.results, self.blends[-1]/(self.sources[-1]+1.e-12*self.blends[-1]))
                    c1 = np.array(c12[0])
                    c2 = np.array(c12[1])
                    cc0 = self.lightcurves[i]
                    ccarr = np.array(cc0[3:7])
                    cctran=np.transpose(ccarr)
                    cc = np.transpose(cctran)
                    sumcN = (cc[0]/(cc[1]*cc[1])).sum()
                    sumc1 = (c1/(cc[1]*cc[1])).sum()
                    sumsigmaN = (1/(cc[1]*cc[1])).sum()
                    sumcE = (cc[2]/(cc[3]*cc[3])).sum()
                    sumc2 = (c2/(cc[3]*cc[3])).sum()
                    sumsigmaE = (1/(cc[3]*cc[3])).sum()
                    L0N = (sumcN-sumc1)/sumsigmaN
                    L0E = (sumcE-sumc2)/sumsigmaE
                    chia = ((cc[0] - c1 - L0N)**2/(cc[1]*cc[1])).sum() + ((cc[2] - c2 - L0E)**2/(cc[3]*cc[3])).sum()      
                    self.chi2 += chia
                self.L0Dec.append(L0N)
                self.L0RA.append(L0E)
                    
            self.blends =np.array(self.blends)
            self.sources =np.array(self.sources)
            self.blendings = self.blends/(self.sources+1.e-12*self.blends)
            self.baselines = -2.5*np.log10(self.blends+self.sources)
                
    def lightcurve(self):
        if(self.modnumber == 1 or self.modnumber == 3 or self.modnumber > 4):
            self.vbm.SetObjectCoordinates(glob.glob('Data/*.coordinates')[0],self.satellitedir)
            self.vbm.parallaxsystem = 1
        if(self.modnumber == 0):
            self.results = self.vbm.ESPLLightCurve(self.pars,self.t)
        elif(self.modnumber == 1):
            if(self.astrometric):
                self.results = self.vbm.ESPLAstroLightCurve(self.pars,self.t)
            else:
                self.results = self.vbm.ESPLLightCurveParallax(self.pars,self.t)
        elif(self.modnumber == 2):
            self.results = self.vbm.BinSourceExtLightCurve(self.pars,self.t)
        elif(self.modnumber == 3):
            if(self.astrometric):
                self.results = self.vbm.BinSourceAstroLightCurveXallarap(self.pars,self.t)
            else:
                self.results = self.vbm.BinSourceExtLightCurveXallarap(self.pars,self.t)
        elif(self.modnumber == 4):
            self.results = self.vbm.BinaryLightCurve(self.pars,self.t)
        elif(self.modnumber == 5):
            if(self.astrometric):
                self.results = self.vbm.BinaryAstroLightCurve(self.pars,self.t)
            else:
                self.results = self.vbm.BinaryLightCurveParallax(self.pars,self.t)
        elif(self.modnumber == 6):
            if(self.astrometric):
                self.results = self.vbm.BinaryAstroLightCurveOrbital(self.pars,self.t)
            else:
                self.results = self.vbm.BinaryLightCurveOrbital(self.pars,self.t)
        elif(self.modnumber == 7):
            if(self.astrometric):
                self.results = self.vbm.BinaryAstroLightCurveKepler(self.pars,self.t)
            else:
                self.results = self.vbm.BinaryLightCurveKepler(self.pars,self.t)
        elif(self.modnumber == 8):
            self.results = self.vbm.TripleLightCurve(self.pars,self.t)
        elif(self.modnumber == 9):
            if(self.astrometric):
                self.results = self.vbm.TripleAstroLightCurve(self.pars,self.t)
            else:
                self.results = self.vbm.TripleLightCurveParallax(self.pars,self.t)
    
    def calculate(self):
        # Light curve calculation
        t0i = self.parnames[self.modnumber].index('t0')
        tEi = self.parnames[self.modnumber].index('tE')
        rhoi = self.parnames[self.modnumber].index('rho')
        self.rho = self.parsprint[rhoi]
        if(self.tmin == ''):
            self.tmin = self.pars[t0i]-2*self.parsprint[tEi]
        if(self.tmax == ''):
            self.tmax = self.pars[t0i]+2*self.parsprint[tEi]
            
        self.lctimes=[]
        self.lcmags=[]
        self.lcerrs=[]
        self.centroids=[]
        self.satellites = []
        for i in range(0,self.nfil):
            lc0 = self.lightcurves[i]
            lcarr = np.array(lc0[:7])
            lctran=np.transpose(lcarr)
            lcsel = [x for x in lctran if(x[0]<self.tmax and x[0]>self.tmin and ((x[1]-self.blends[i])/(self.sources[i]+1.e-12*self.blends[i])*self.sources[self.referencephot] +self.blends[self.referencephot]>0 or x[2] < 0))]
            lc = np.transpose(lcsel)
            if(len(lc)>0):
                self.lctimes.append(lc[0])
                if(lc[2][0]>0):
                    self.lcmags.append(np.array([-2.5*math.log10((y-self.blends[i])/(self.sources[i]+1.e-12*self.blends[i])*self.sources[self.referencephot]+self.blends[self.referencephot]) for y in lc[1]]))
                    self.lcerrs.append(lc[2]/lc[1]*2.5/math.log(10.0))
                else:
                    self.lcmags.append([])
                    self.lcerrs.append([])
                self.centroids.append([np.array(lc[3]),np.array(lc[4]),np.array(lc[5]),np.array(lc[6])])
            else:
                self.lctimes.append(np.array([]))
                self.lcmags.append(np.array([]))
                self.lcerrs.append(np.array([]))
                self.centroids.append([])
            self.satellites.append(lc0[7])

        
        self.t0 = np.linspace(self.tmin,self.tmax,self.timesteps)
        self.usedsatellites = list(set(self.satellites))
        if(self.eventname == None):
            self.usedsatellites = [0]
            self.referencephot = 0
            self.sources = [1]
            self.blends = [0]
        while(len(self.satellitecolors)<len(self.usedsatellites)):
            self.satellitecolors.extend(self.satellitecolors)
        self.minmag = 1000
        self.maxmag = -1000
        self.maxy1 = -1000
        self.maxy2 = -1000
        self.miny1 = 1000
        self.miny2 = 1000
        self.magnitudes = []
        self.magnifications = []
        self.trajectories = []
        self.c1s = []
        self.c2s = []
        self.c1l = []
        self.c2l = []
        for satellite in self.usedsatellites:   
            self.t =self.t0[:]
#            for i in range(self.nfil):
#                if(self.satellites[i] == satellite):
#                    self.t = np.concatenate((self.t,self.lctimes[i]))                
            self.t = np.sort(self.t)
            self.vbm.satellite = satellite
            self.vbm.a1 = self.limbdarkenings[self.referencephot]
            self.lightcurve()        
            self.mags = [-2.5*math.log10(max(self.sources[self.referencephot]*yi+self.blends[self.referencephot], 1.e-100)) for yi in self.results[0]]
            if(self.astrometric):
                self.magnifications.append(self.results[0])
                self.c1s.append(self.results[1])
                self.c2s.append(self.results[2])
                self.c1l.append(self.results[3])
                self.c2l.append(self.results[4])
                self.y1 = self.results[5]
                self.y2 = self.results[6]
            else:
                self.y1 = self.results[1]
                self.y2 = self.results[2]
            self.magnitudes.append([self.t,self.mags])
            self.trajectories.append([self.y1,self.y2])
            minmag = min(self.mags)
            maxmag = max(self.mags)
            self.minmag = min(self.minmag, minmag)
            self.maxmag = max(self.maxmag, maxmag)            
            miny = min(self.y1)
            maxy = max(self.y1)
            self.miny1 = min(self.miny1, miny)
            self.maxy1 = max(self.maxy1, maxy)            
            miny = min(self.y2)
            maxy = max(self.y2)
            self.miny2 = min(self.miny2, miny)
            self.maxy2 = max(self.maxy2, maxy)            
        margin = (self.maxmag-self.minmag)*0.1
        self.minmag -= margin
        self.maxmag += margin
        if(self.inputmagmin != ''):
            self.minmag = self.inputmagmin
        if(self.inputmagmax != ''):
            self.maxmag = self.inputmagmax

        self.lcress = []
        for i in range(self.nfil):
            self.t = self.lctimes[i]                
            self.vbm.satellite = self.satellites[i]
            self.vbm.a1 = self.limbdarkenings[i]
            if(len(self.lcmags[i]))>0:
                self.lightcurve()
                self.mags = [-2.5*math.log10(max(self.sources[self.referencephot]*yi+self.blends[self.referencephot],1.e-100)) for yi in self.results[0]]
                ress = self.mags-self.lcmags[i]
                self.lcress.append(ress)
            else:
                self.lcress.append([])

    # Caustic plot preparation   
        self.rancau = max([self.maxy1-self.miny1,self.maxy2 - self.miny2])
        if(self.modnumber>7):
            self.caus = self.vbm.Multicaustics()
        elif(self.modnumber>3):
            self.caus = self.vbm.Caustics(self.parsprint[0],self.parsprint[1])
        else:
            self.caus = np.array([[[0,1,0,-1,0],[1,0,-1,0,1]]])*self.rancau*0.001        
        # tE=parsprint[parnames[modnumber].index('tE')]

    def printparameters(self):
        self.parstring = ''
        table =[['chi2', str(self.chi2)]]
        for i in range(self.npars[self.modnumber]):
            if((not self.animate) and len(self.parameters)==0):
                table.append([self.parnames[self.modnumber][i], self.approx(i)])
#                    self.parstring = self.parstring + self.parnames[self.modnumber][i] + ' = ' + self.approx(i) + '\n'  #+ str(self.parsprint[i]) +  ' +- ' + str(self.parerrs[i]) + '\n'
            else:
                table.append([self.parnames[self.modnumber][i], str(self.parsprint[i])])
#                self.parstring = self.parstring + self.parnames[self.modnumber][i] + ' = ' + str(self.parsprint[i]) + '\n'
        self.parstring = self.parstring + tabulate(table, headers='firstrow', tablefmt='fancy_grid')
        self.parstring = self.parstring + '\n\n'
        
        if((not self.animate) and len(self.parameters)==0):
            table = [[self.telescopes[i], self.approxbaselines(i), self.approxblendings(i)] for i in range(self.nfil)]
            table.insert(0,['telescope','baseline', 'blending'])
        else:
            table = [[t,base, bl] for t,base,bl in zip(self.telescopes,self.baselines,self.blendings)]
            table.insert(0,['telescope','baseline', 'blending'])
        for i in range(self.nfil,0,-1):
            if(len(self.lcmags[i-1])==0):
                del(table[i])
        self.parstring = self.parstring + tabulate(table, headers='firstrow', tablefmt='fancy_grid')
        print(self.parstring)

    def fexp(self, f):
        return int(math.floor(math.log10(abs(f)))) if f != 0 else 0
    
    def approx(self, i):
        exerr= self.fexp(self.parerrs[i])-1
        expars= self.fexp(self.parsprint[i])-1
        return f'{self.parsprint[i]:.{max(0,-exerr,-expars)}f}' + ' +- ' + f'{self.parerrs[i]:.{max(0,-exerr,-expars)}f}'

    def approxbaselines(self, i):
        exerr= self.fexp(self.baselineerrs[i])-1
        return f'{self.baselines[i]:.{max(0,-exerr)}f}' + ' +- ' + f'{self.baselineerrs[i]:.{max(0,-exerr)}f}'

    def approxblendings(self, i):
        exerr= self.fexp(self.blendingerrs[i])-1
        return f'{self.blendings[i]:.{max(0,-exerr)}f}' + ' +- ' + f'{self.blendingerrs[i]:.{max(0,-exerr)}f}'

    def axeslightcurve(self,ax):
        for i in range(0,self.nfil):
            if(len(self.lcmags[i])>0):
                ax.errorbar(self.lctimes[i],self.lcmags[i],yerr=self.lcerrs[i],color=self.colors[i],fmt='.',label=self.telescopes[i])
        for i in range(len(self.usedsatellites)):
            ax.plot(self.magnitudes[i][0],self.magnitudes[i][1],self.satellitecolors[i],linewidth=0.5)
        ax.set_ylabel(self.maglabel)
        ax.set_ylim([self.maxmag,self.minmag])
        ax.set_xlim([self.tmin,self.tmax])
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        if(self.eventname != None):
            ax.legend(loc=self.legendlocation)

    def axesresiduals(self,ax):
        ax.plot((self.tmin,self.tmax),(0,0),'k-',linewidth=0.5)
        for i in range(0,self.nfil):
            if(len(self.lcress[i])>0):
                ax.errorbar(self.lctimes[i],self.lcress[i],yerr=self.lcerrs[i],color=self.colors[i],fmt='.')
        ax.set_xlabel(self.tlabel)
        ax.set_ylabel(self.reslabel)
        ax.set_ylim([-0.1,0.1])
        ax.set_xlim([self.tmin,self.tmax])
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

    def showlightcurve(self):
        plt.figure()
        fig, axs =plt.subplots(2,1,figsize=(12,6), gridspec_kw={'height_ratios': [5, 1]},sharex='col')
        self.axeslightcurve(axs[0])
        self.axesresiduals(axs[1])
        plt.subplots_adjust(hspace=0.1)
        self.figure = fig

    def showicon(self,outputfile):
        plt.figure()
        fig, ax =plt.subplots(1,1,figsize=(6,6))
        for i in range(0,self.nfil):
            ax.errorbar(self.lctimes[i],self.lcmags[i],yerr=self.lcerrs[i],color=self.colors[i],fmt='.',label=self.telescopes[i])
        for i in range(len(self.usedsatellites)):
            ax.plot(self.magnitudes[i][0],self.magnitudes[i][1],self.satellitecolors[i],linewidth=2.5)
        ax.set_ylim([self.maxmag,self.minmag])
        ax.set_xlim([self.tmin,self.tmax])
        ax.tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False) 
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(6)  # change width
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        self.figure = fig
        self.figure.savefig(outputfile)
        plt.close()
  
    def showcaustics(self):
        plt.figure()
        fig, ax =plt.subplots(figsize=[5,5])
        self.axescaustics(ax)
        self.figure = fig

    def axescaustics(self,ax):
        for cau in self.caus:
            ax.plot(cau[0],cau[1],color = self.causticcolor)
        ax.set_xlabel('y1')
        ax.set_ylabel('y2')
        ax.set_xlim([(self.miny1+self.maxy1-self.rancau)/2,(self.miny1+self.maxy1+self.rancau)/2])
        ax.set_ylim([(self.miny2+self.maxy2-self.rancau)/2,(self.miny2+self.maxy2+self.rancau)/2])
        for i in range(len(self.usedsatellites)):
            ax.plot(self.trajectories[i][0],self.trajectories[i][1],self.satellitecolors[i])
            n = len(self.trajectories[i][0])
            arrlength = self.rancau*0.01
            i0 = n-10
            dir = np.array([self.trajectories[i][0][i0+1]-self.trajectories[i][0][i0],self.trajectories[i][1][i0+1]-self.trajectories[i][1][i0]])
            dir = dir/(math.sqrt(dir.dot(dir)))
            dirort = np.array([-dir[1], dir[0]])
            p0 = np.array([self.trajectories[i][0][i0],self.trajectories[i][1][i0]])
            p1 = (-dir +dirort)*arrlength + p0
            p2 = (-dir -dirort)*arrlength + p0
            ax.plot((p0[0],p1[0]),(p0[1],p1[1]),self.satellitecolors[i])
            ax.plot((p0[0],p2[0]),(p0[1],p2[1]),self.satellitecolors[i])
            i0 = 20
            dir = np.array([self.trajectories[i][0][i0+1]-self.trajectories[i][0][i0],self.trajectories[i][1][i0+1]-self.trajectories[i][1][i0]])
            dir = dir/(math.sqrt(dir.dot(dir)))
            dirort = np.array([-dir[1], dir[0]])
            p0 = np.array([self.trajectories[i][0][i0],self.trajectories[i][1][i0]])
            p1 = (-dir +dirort)*arrlength + p0
            p2 = (-dir -dirort)*arrlength + p0
            ax.plot((p0[0],p1[0]),(p0[1],p1[1]),self.satellitecolors[i])
            ax.plot((p0[0],p2[0]),(p0[1],p2[1]),self.satellitecolors[i])
            i0 = int(n*0.5)
            p0 = np.array([self.trajectories[i][0][i0],self.trajectories[i][1][i0]])
            ax.add_patch(Circle(p0,self.rho,color = self.sourcecolor))
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

    def showall(self):
        if(self.eventname != None):
            fig, axes = plt.subplot_mosaic([['tl','right'],['bl','right']],figsize=(12,5.5), gridspec_kw={'height_ratios': [5, 1]})
            self.axesresiduals(axes['bl'])
            axes['tl'].xaxis.set_ticklabels([])
            fig.suptitle(self.model);
        else:
            fig, axes = plt.subplot_mosaic([['tl','right']],figsize=(12,5.5))
        self.axeslightcurve(axes['tl'])
        self.axescaustics(axes['right'])
        plt.subplots_adjust(hspace=0.1)
        plt.show()
        self.figure = fig
        if(self.printpars):
            self.printparameters()

    def calculate_astrometry(self, i):
        isat = self.satellites[i]
        mags = np.array(self.magnifications[isat][1])
        c1s = np.array(self.c1s[isat])
        c2s = np.array(self.c2s[isat])
        c1l = np.array(self.c1l[isat])
        c2l = np.array(self.c2l[isat])
        self.centroidi = np.array([(mags*c1s+self.blendings[i]*c1l)/(mags+self.blendings[i])+self.L0Dec[i],\
                             (mags*c2s+self.blendings[i]*c2l)/(mags+self.blendings[i])+self.L0RA[i]])
        maxDec = max(self.centroidi[0])
        minDec = min(self.centroidi[0])
        maxRA = max(self.centroidi[1])
        minRA = min(self.centroidi[1])
        meanerrDec = np.median(self.centroids[i][1])
        meanerrRA = np.median(self.centroids[i][3])
        self.ran = max(maxDec-minDec+3*meanerrDec,maxRA-minRA+3*meanerrRA)*0.55
        self.cenDec = 0.5*(maxDec+minDec)
        self.cenRA = 0.5*(maxRA+minRA)

    def showastrometry(self, i = None):
        plt.figure()
        fig, ax =plt.subplots(figsize=[5,5])
        self.axesastrometry(ax, i)
        self.figure = fig
    
    def axesastrometry(self,ax, i):
        if(i==None):
            for j in range(self.nfil):
                if(len(self.centroids[j])>0 and self.centroids[j][1][0]>0):
                    i=j
                    break
            if(i==None):
                print('No astrometric data!')
                return
        if(self.centroids[i][1][0]>0):          
            self.calculate_astrometry(i)
            ax.set_ylabel('Dec')
            ax.set_xlabel('RA')
            ax.set_ylim([self.cenDec-self.ran,self.cenDec+self.ran])
            ax.set_xlim([self.cenRA+self.ran,self.cenRA-self.ran])
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.errorbar(self.centroids[i][2],self.centroids[i][0],xerr=self.centroids[i][3],yerr=self.centroids[i][1],color=self.colors[i],fmt='.',label=self.telescopes[i])
            ax.plot(self.centroidi[1],self.centroidi[0],self.satellitecolors[self.satellites[i]],linewidth=1,zorder = 1.e10)
            if(self.eventname != None):
                ax.legend(loc=self.legendlocation)
        else:
            print('No astrometry for this dataset')

    def showastrometryRA(self, i = None):
        plt.figure()
        fig, ax =plt.subplots(figsize=[7,5])
        self.axesastrometryRA(ax, i)
        self.figure = fig
    
    def axesastrometryRA(self,ax, i):
        if(i==None):
            for j in range(self.nfil):
                if(len(self.centroids[j])>0 and self.centroids[j][1][0]>0):
                    i=j
                    break
            if(i==None):
                print('No astrometric data!')
                return
        if(self.centroids[i][1][0]>0):          
            self.calculate_astrometry(i)
            ax.set_ylabel('RA')
            ax.set_xlabel(self.tlabel)
            ax.set_ylim([self.cenRA-self.ran,self.cenRA+self.ran])
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.errorbar(self.lctimes[i], self.centroids[i][2],yerr=self.centroids[i][3],color=self.colors[i],fmt='.',label=self.telescopes[i])
            ax.plot(self.magnitudes[self.satellites[i]][0],self.centroidi[1],self.satellitecolors[self.satellites[i]],linewidth=1,zorder = 1.e10)
            if(self.eventname != None):
                ax.legend(loc=self.legendlocation)
        else:
            print('No astrometry for this dataset')
   
    def showastrometryDec(self, i = None):
        plt.figure()
        fig, ax =plt.subplots(figsize=[7,5])
        self.axesastrometryDec(ax, i)
        self.figure = fig

    def axesastrometryDec(self,ax, i):
        if(i==None):
            for j in range(self.nfil):
                if(len(self.centroids[j])>0 and self.centroids[j][1][0]>0):
                    i=j
                    break
            if(i==None):
                print('No astrometric data!')
                return
        if(self.centroids[i][1][0]>0):          
            self.calculate_astrometry(i)
            ax.set_ylabel('Dec')
            ax.set_xlabel(self.tlabel)
            ax.set_ylim([self.cenDec-self.ran,self.cenDec+self.ran])
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.errorbar(self.lctimes[i], self.centroids[i][0],yerr=self.centroids[i][1],color=self.colors[i],fmt='.',label=self.telescopes[i])
            ax.plot(self.magnitudes[self.satellites[i]][0],self.centroidi[0],self.satellitecolors[self.satellites[i]],linewidth=1,zorder = 1.e10)
            if(self.eventname != None):
                ax.legend(loc=self.legendlocation)
        else:
            print('No astrometry for this dataset')
    
    def update(self, frame):
        self.im.set_array(self.images_array[frame])
        return self.im,
    
    def animateplots(self, interval = 1000):
        if(os.path.exists('tmpsteps')):
            shutil.rmtree('tmpsteps')
        os.mkdir('tmpsteps')
        os.chdir('tmpsteps')
        dirocco = os.getcwd()
        self.line = 0
        self.readdata()
        with open(self.modelfile) as f:
            lines = f.readlines()
            self.maxlines = len(lines)
        pbar = tqdm(total = self.maxlines,desc = 'Frames',file=sys.stdout, colour='GREEN', smoothing = 0)
        for frame in range(self.maxlines):
            self.line = frame
            self.readparameters()
            self.fig, self.axes = plt.subplot_mosaic([['tl','right'],['bl','right']],figsize=(12,6), gridspec_kw={'height_ratios': [5, 1]})
            self.fig.suptitle(self.model + ' - frame: ' + str(self.line));
            self.calculate()
            self.axeslightcurve(self.axes['tl'])
            self.axes['tl'].xaxis.set_ticklabels([])
            self.axesresiduals(self.axes['bl'])
            self.axescaustics(self.axes['right'])
            _ = plt.subplots_adjust(hspace=0.1)
            os.chdir(dirocco)
            plt.savefig(str(self.line) + '.png', bbox_inches = 'tight', dpi = 300)
            plt.close(self.fig)
            pbar.update(1)
        pbar.close()
        self.images_array = []
        for i in range(self.maxlines):
            image = Image.open(str(i) + '.png')
            self.images_array.append(image)
        plt.close(self.fig)
        self.fig, ax = plt.subplots()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        self.fig.subplots_adjust(bottom=0, top=1, left=0, right=1)
        self.im = ax.imshow(self.images_array[0], animated=True)
        self.animation_fig = animation.FuncAnimation(self.fig, self.update, frames=len(self.images_array), interval=interval, blit=True,repeat = False)
        os.chdir(dirocco)
        os.chdir('..')
        self.animation_fig.save('ani.gif',dpi = 150)        
        plt.close(self.fig)
        
def plotchain(eventname, model, par1, par2):
    chains = []
    filenames = glob.glob(eventname+ '/PreModels/' + model + '-step*')
    for fil in filenames:
        with open(fil) as f:
            chainstrings = f.readlines()
        chainlist = []
        for st in chainstrings:
            chunks = st.split(' ')
            chainlist.append([float(v) for v in chunks])
        chain = np.array(chainlist)
        chains.append(chain)
    colors = ['blue','red','green','darkorange','magenta','cyan','gray','teal','maroon','gold','lime','darkviolet']
    while(len(colors)<len(chains)):
        colors.extend(colors)
    fig, ax = plt.subplots()
#    ax.set_xlabel('s')
#    ax.set_ylabel('q')
    for i in range(len(filenames)):
        chain = chains[i]
        x = chain.transpose()[par1]
        y = chain.transpose()[par2]
        ax.plot(x,y,color = colors[i])
        ax.scatter(x[-1], y[-1],s=20,color = colors[i])
        
def orbital_elements(modelfile):
    with open(modelfile) as f:
        line=f.readline().split()
        parsall=[float(chunk) for chunk in line]
        orbitalparameters = None
        if(os.path.basename(modelfile)[0]=='L'):
            if(os.path.basename(modelfile)[1]=='O'):
                parameters = parsall[0:12]
                w1 = parameters[9]
                w2 = parameters[10]
                w3 = parameters[11]
                s = parameters[0]
                t0 = parameters[6]
                alpha = parameters[3]
                calpha = math.cos(alpha)
                salpha = math.sin(alpha)
                w13 = w1 * w1 + w3 * w3
                w123 = math.sqrt(w13 + w2 * w2)
                w13 = math.sqrt(w13)
                if (w13 > 1.e-8):
                    if(w3 < 1.e-8): 
                        w3 = 1.e-8
                    w = w3 * w123 / w13
                    s3d = s*w13/w3
                    inc = math.acos(w2 * w3 / w13 / w123)
                    Om = math.atan2(w1*w2/w13,w13)
                    phi0 = math.atan2(-w1 / w3, w13 / w123)
                else:
                    w = w2
                    inc = 0.0
                    Om = math.pi/2
                    phi0 = -Om
                orbitalparameters = {'T': 2*math.pi/w, 'a': s3d, 'e': 0, 'inc': inc, 'OM': Om, 'om': 0, 'phi0': phi0,'epoch': t0-phi0*w}
            elif(os.path.basename(modelfile)[1]=='K'):
                parameters = parsall[0:14]
                w1 = parameters[9]
                w2 = parameters[10]
                w3 = parameters[11]
                szs = parameters[12]
                ar = parameters[13]
                s = parameters[0]
                t0 = parameters[6]
                alpha = parameters[3]
                calpha = math.cos(alpha)
                salpha = math.sin(alpha)
                smix = 1 + szs * szs
                sqsmix = math.sqrt(smix)
                w22 = w2 * w2
                w11 = w1 * w1
                w33 = w3 * w3
                w12 = w11 + w22
                w23 = w22 + w33
                wt2 = w12 + w33
                
                szs2 = szs * szs
                ar2 = ar * ar
                arm1 = ar - 1
                arm2 = 2 * ar - 1
                n = math.sqrt(wt2 / arm2 / smix) / ar
                Z0 = -szs * w2
                Z1 = szs * w1 - w3
                Z2 = w2
                h = math.sqrt(Z0 * Z0 + Z1 * Z1 + Z2 * Z2)
                Z0 /= h
                Z1 /= h
                Z2 /= h
                X0 = -ar * w11 + arm1 * w22 - arm2 * szs * w1 * w3 + arm1 * w33
                X1 = -arm2 * w2 * (w1 + szs * w3)
                X2 = arm1 * szs * w12 - arm2 * w1 * w3 - ar * szs * w33
                e = math.sqrt(X0 * X0 + X1 * X1 + X2 * X2)
                X0 /= e
                X1 /= e
                X2 /= e
                e /= ar * sqsmix * wt2
                a = ar * s * math.sqrt(smix)
                Y0 = Z1 * X2 - Z2 * X1;
                Y1 = Z2 * X0 - Z0 * X2;
                Y2 = Z0 * X1 - Z1 * X0;
                
                cosnu = (X0 + X2 * szs) / sqsmix
                sinnu = Y0 + Y2 * szs
        
                co1EE0 = cosnu + e
                co2EE0 = 1 + e * cosnu
                cosE = co1EE0 / co2EE0
                EE0 = math.acos(cosE)
                EE0 *= np.sign(sinnu)
                sinE = math.sqrt(1 - cosE * cosE) * np.sign(sinnu);
                co1tperi = e * sinE;
                tperi = t0 - (EE0 - co1tperi) / n;
                inc = math.acos(Z2)
                Om = math.atan2(Z0,-Z1)
                om = math.acos((X1*Z0-X0*Z1)/math.sin(inc))*np.sign(X2)
                phi0 = math.acos(cosnu)*np.sign(sinnu)       
                orbitalparameters = {'T': 2*math.pi/n, 'a': a, 'e': e, 'inc': inc, 'OM': Om, 'om': om, 'phi0': phi0, 'epoch': tperi}
    return orbitalparameters
