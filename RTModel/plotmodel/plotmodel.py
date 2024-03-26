import VBBinaryLensing
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


class plotmodel:
    def __init__(self, eventname,model = '', tmin = '', tmax = '', referencephot = 0, timesteps = 300, \
                 modelfile = None, parameters = [], line = 0,printpars = True, animate = False,interval = 1000, satellitedir = '.'):
        self.satellitedir = satellitedir
        self.parameters = parameters
        filin=inspect.getfile(VBBinaryLensing)
        self.filout= os.path.dirname(filin) + '/data/ESPL.tbl'
        self.eventname = eventname
        self.model = model
        self.tmin = tmin
        self.tmax = tmax
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
        self.vbbl = VBBinaryLensing.VBBinaryLensing()
        # General information on models
        self.modelcodes= ['PS','PX','BS','BO','LS','LX','LO']
        self.npars=[4,6,7,10,7,9,12]
        self.logposs=[[0,1,3],
                 [1,3],
                 [0,1,6],
                 [2,3,9],
                 [0,1,4,5],
                 [0,1,4,5],
                 [0,1,4,5]]
        self.parnames = [['u0','tE','t0','rho'],
                    ['u0','tE','t0','rho','piN','piE'],
                    ['tE','FR','u01','u02','t0','t02','rho'],
                    ['u0','t0','tE','rho','xi1','xi2','om','inc','phi','qs'],
                    ['s','q','u0','alpha','rho','tE','t0'],
                    ['s','q','u0','alpha','rho','tE','t0','piN','piE'],
                    ['s','q','u0','alpha','rho','tE','t0','piN','piE','gamma1','gamma2','gammaz']]
        self.colors = ['blue','red','green','darkorange','magenta','cyan','gray','teal','maroon','gold','lime','darkviolet']
        self.sourcecolor = 'pink'
        self.causticcolor = (0.5,0.5,0.5)
        self.satellitecolors = ['k-','r-','g-','b-','y-']
        self.legendlocation = 'best'
        if(animate):
            self.animateplots(interval = interval)
        else:
            self.readdata()
            self.readparameters()
            self.calculate()
            self.showall()
        
    # Reading data from LCToFit.txt
    def readdata(self):
        if(self.eventname== None):
            self.lightcurves = []
            self.telescopes = []
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
                    data[self.nfil].append([float(chunks[1]),float(chunks[2]),float(chunks[3]),int(chunks[4])])
                self.nfil +=1
            self.lightcurves = [ [np.array([dl[0] for dl in d]),np.array([dl[1] for dl in d]),np.array([dl[2] for dl in d]),d[0][3]] for d in data]
            with open('FilterToData.txt') as f:
                self.telescopes = f.readlines()
                for i in range(0,self.nfil):
                    self.telescopes[i] = self.telescopes[i][0:self.telescopes[i].index('.')]
            while(len(self.colors)<self.nfil):
                self.colors.extend(self.colors)

    # Reading model parameters
    def readparameters(self):
        self.modnumber = self.modelcodes.index(self.model[0:2])
        if(self.eventname != None):
            os.chdir(self.eventname)
        if(self.parameters == []):
            with open(self.modelfile) as f:
                lines = f.readlines()
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
                    self.blends = np.array([values[self.npars[self.modnumber]+i*2] for i in range(0,self.nfil)])
                    self.sources = np.array([values[self.npars[self.modnumber]+i*2+1] for i in range(0,self.nfil)])
                    self.chi2=values[-1]
                    #Errors
                    if(not self.animate):
                        line = lines[self.line +1 ]
                        chunks = line.split(' ')
                        values = [float(v) for v in chunks]
                        self.parerrs = values[0:self.npars[self.modnumber]]
                        self.blenderrs = np.array([values[self.npars[self.modnumber]+i*2] for i in range(0,self.nfil)])
                        self.sourceerrs = np.array([values[self.npars[self.modnumber]+i*2+1] for i in range(0,self.nfil)])
                    return True
        else:
            self.pars = self.parameters[0:self.npars[self.modnumber]]
            self.parsprint=self.pars[:] 
            for i in self.logposs[self.modnumber]:
                self.pars[i] = math.log(self.pars[i])
            self.blends = []
            self.sources = []
            self.chi2 = 0
            for i in range(0,self.nfil):
                lc0 = self.lightcurves[i]
                lcarr = np.array(lc0[0:3])
                lctran=np.transpose(lcarr)
                lc = np.transpose(lctran)                
                self.t = lc[0]
                self.vbbl.satellite = lc0[3]
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
                self.chi2 += sumy2+(sumfy*sumfy*sumsigma+sumf2*sumy*sumy-2*sumf*sumy*sumfy)/p1
                
    def lightcurve(self):
        if(self.modnumber < 4):
            self.vbbl.LoadESPLTable(self.filout)
        if(self.modnumber == 1 or self.modnumber > 4):
            self.vbbl.SetObjectCoordinates(glob.glob('Data/*.coordinates')[0],self.satellitedir)
            self.vbbl.parallaxsystem = 1
        if(self.modnumber == 0):
            self.results = self.vbbl.ESPLLightCurve(self.pars,self.t)
        elif(self.modnumber == 1):
            self.results = self.vbbl.ESPLLightCurveParallax(self.pars,self.t)
        elif(self.modnumber == 2):
            self.results = self.vbbl.BinSourceExtLightCurve(self.pars,self.t)
        elif(self.modnumber == 3):
            self.results = self.vbbl.BinSourceSingleLensXallarap(self.pars,self.t)
        elif(self.modnumber == 4):
            self.results = self.vbbl.BinaryLightCurve(self.pars,self.t)
        elif(self.modnumber == 5):
            self.results = self.vbbl.BinaryLightCurveParallax(self.pars,self.t)
        elif(self.modnumber == 6):
            self.results = self.vbbl.BinaryLightCurveOrbital(self.pars,self.t)

    
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
        self.satellites = []
        for i in range(0,self.nfil):
            lc0 = self.lightcurves[i]
            lcarr = np.array(lc0[0:3])
            lctran=np.transpose(lcarr)
            lcsel = [x for x in lctran if(x[0]<self.tmax and x[0]>self.tmin and (x[1]-self.blends[i])/self.sources[i]*self.sources[self.referencephot] +self.blends[self.referencephot]>0)]
            lc = np.transpose(lcsel)
            if(len(lc)>0):
                self.lctimes.append(lc[0])
                self.lcmags.append(np.array([-2.5*math.log10((y-self.blends[i])/self.sources[i]*self.sources[self.referencephot]+self.blends[self.referencephot]) for y in lc[1]]))
                self.lcerrs.append(lc[2]/lc[1]*2.5/math.log(10.0))
            else:
                self.lctimes.append(np.array([]))
                self.lcmags.append(np.array([]))
                self.lcerrs.append(np.array([]))
            self.satellites.append(lc0[3])

        
        self.t0 = np.linspace(self.tmin,self.tmax,self.timesteps)
        self.usedsatellites = list(set(self.satellites))
        while(len(self.satellitecolors)<len(self.usedsatellites)):
            self.satellitecolors.extend(self.satellitecolors)
        self.minmag = 1000
        self.maxmag = -1000
        self.maxy1 = -1000
        self.maxy2 = -1000
        self.miny1 = 1000
        self.miny2 = 1000
        self.magnifications = []
        self.trajectories = []
        for satellite in self.usedsatellites:   
            self.t =self.t0[:]
            for i in range(self.nfil):
                if(self.satellites[i] == satellite):
                    self.t = np.concatenate((self.t,self.lctimes[i]))                
            self.t = np.sort(self.t)
            self.vbbl.satellite = satellite
            self.lightcurve()        
            self.mags = [-2.5*math.log10(self.sources[self.referencephot]*yi+self.blends[self.referencephot]) for yi in self.results[0]]
            self.y1 = self.results[1]
            self.y2 = self.results[2]
            self.magnifications.append([self.t,self.mags])
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

        self.lcress = []
        for i in range(self.nfil):
            ress = []
            isat = np.where(np.array(self.usedsatellites) == self.satellites[i])[0][0]
            for j in range(0,len(self.lctimes[i])):
                ip = np.where(self.magnifications[isat][0] == self.lctimes[i][j])[0][0]
                ress.append(self.magnifications[isat][1][ip]-self.lcmags[i][j])
            self.lcress.append(ress)

    # Caustic plot preparation   
        self.rancau = max([self.maxy1-self.miny1,self.maxy2 - self.miny2])
        if(self.modnumber>3):
            self.caus = self.vbbl.Caustics(self.parsprint[0],self.parsprint[1])
        else:
            self.caus = np.array([[[0,1,0,-1,0],[1,0,-1,0,1]]])*self.rancau*0.001        
        # tE=parsprint[parnames[modnumber].index('tE')]

    def printparameters(self):
        print('Parameters')
        for i in range(self.npars[self.modnumber]):
            if((not self.animate) and len(self.parameters)==0):
                    print(self.parnames[self.modnumber][i],' = ', self.parsprint[i], ' +- ', self.parerrs[i])
            else:
                print(self.parnames[self.modnumber][i],' = ', self.parsprint[i])
        print()
        print('blending = ', np.array(self.blends)/np.array(self.sources))
        print('baseline = ', -2.5*np.log10(np.array(self.blends)+np.array(self.sources)))
        print('chi2 =', self.chi2)

    def axeslightcurve(self,ax):
        for i in range(0,self.nfil):
            ax.errorbar(self.lctimes[i],self.lcmags[i],yerr=self.lcerrs[i],color=self.colors[i],fmt='.',label=self.telescopes[i])
        for i in range(len(self.usedsatellites)):
            ax.plot(self.magnifications[i][0],self.magnifications[i][1],self.satellitecolors[i],linewidth=0.5)
        ax.set_ylabel('mag')
        ax.set_ylim([self.maxmag,self.minmag])
        ax.set_xlim([self.tmin,self.tmax])
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.legend(loc=self.legendlocation)

    def axesresiduals(self,ax):
        ax.plot((self.tmin,self.tmax),(0,0),'k-',linewidth=0.5)
        for i in range(0,self.nfil):
            ax.errorbar(self.lctimes[i],self.lcress[i],yerr=self.lcerrs[i],color=self.colors[i],fmt='.')
        ax.set_xlabel('t')
        ax.set_ylabel('Res')
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

    
    def showcaustics(self):
        plt.figure()
        fig, ax =plt.subplots(figsize=[5,5])
        self.axescaustics(ax)

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
        fig, axes = plt.subplot_mosaic([['tl','right'],['bl','right']],figsize=(12,5.5), gridspec_kw={'height_ratios': [5, 1]})
        fig.suptitle(self.model);
        self.axeslightcurve(axes['tl'])
        axes['tl'].xaxis.set_ticklabels([])
        self.axesresiduals(axes['bl'])
        self.axescaustics(axes['right'])
        plt.subplots_adjust(hspace=0.1)
        plt.show()
        if(self.printpars):
            self.printparameters()
    
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
    filenames = glob.glob(eventname+ '/PreModels/' + model + '/*step*')
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
    





