
import subprocess
import os
import sys
import glob
import time
import inspect
from tqdm import tqdm
import shutil


class RTModel:
    def __init__(self, event = None):
        # Directory preliminaries
        print('*********************')
        print('****   RTModel   ****')
        print('*********************')
        self.pathtoRTM = inspect.getfile(RTModel)
        self.bindir = os.path.dirname(self.pathtoRTM) + '/bin/'
        if os.name =='nt':
            self.readerexe = 'Reader.exe'
            self.initcondexe = 'InitCond.exe'
            self.levmarexe = 'LevMar.exe'
            self.modelselectorexe = 'ModelSelector.exe'
            self.finalizerexe = 'Finalizer.exe'
        else:
            self.readerexe = 'Reader'
            self.initcondexe = 'InitCond'
            self.levmarexe = 'LevMar'
            self.modelselectorexe = 'ModelSelector'
            self.finalizerexe = 'Finalizer'
        if(event == None):
            self.eventname = os.getcwd()
        else:
            self.eventname = os.path.realpath(event)
            print("Event name: " + self.eventname)
        self.inidir = "ini"
        self.modelcodes = ['PS', 'PX', 'BS', 'BO', 'LS', 'LX', 'LO']
        self.endphase =  len(self.modelcodes)*2+3
        self.eventinifile = 'event.ini'        
        self.nprocessors = os.cpu_count()
        print('Number of processors: {}'.format(self.nprocessors))
        self.config_Reader()
        self.config_InitCond()
        self.config_LevMar()
        self.config_ModelSelector()
        self.satellitedir = '.'

    def set_processors(self, nprocessors):
        self.nprocessors = nprocessors

    def set_event(self, event):
        self.eventname = os.path.realpath(event)

    def set_satellite_dir(self, satellitedir):
        self.satellitedir = satellitedir

    def config_Reader(self, tau = 0.1, binning = 4000, otherseasons = 1, renormalize = 1, thresholdoutliers = 10):
        self.Reader_tau= tau # conventional correlation time for consecutive points
        self.Reader_binning = binning # maximum number of points left after re-binning
        self.Reader_otherseasons = otherseasons # How to use other seasons (0 = Yes, 1 = decrease significance, 2 = remove)
        self.Reader_renormalize = renormalize # Re-normalize error bars if non-zero
        self.Reader_thresholdoutliers = thresholdoutliers # Threshold in sigmas for removing outliers
        
    def Reader(self):
        if(not os.path.exists(self.eventname + '/' + self.inidir)):
            os.makedirs(self.eventname + '/' + self.inidir)
        with open(self.eventname + '/' + self.inidir + '/Reader.ini','w') as f:
            f.write('tau = ' + str(self.Reader_tau) + '\n')
            f.write('binning = ' + str(self.Reader_binning) + '\n')
            f.write('otherseasons = ' + str(self.Reader_otherseasons) + '\n')
            f.write('renormalize = ' + str(self.Reader_renormalize) + '\n')
            f.write('thresholdoutliers = ' + str(self.Reader_thresholdoutliers) + '\n')
        print('- Launching: Reader')
        print('  Pre-processing data...')
        completedprocess = subprocess.run([self.bindir+self.readerexe,self.eventname], cwd = self.bindir, shell = False, stdout=subprocess.DEVNULL)
        if(completedprocess.returncode != 0):
            print('! Error in pre-processing. Please check your data!')            
            self.done = True
        else:
            print('  OK')

    def config_InitCond(self, npeaks = 2, peakthreshold = 10.0, oldmodels = 4, override = None, nostatic = False, onlyorbital = False, usesatellite = 0):
        self.InitCond_npeaks = npeaks # Number of peaks in the observed light curve to be considered for setting initial conditions.
        self.InitCond_peakthreshold = peakthreshold # Number of sigmas necessary for a deviation to be identified as a maximum or a minimum.
        self.InitCond_oldmodels = oldmodels # Maximum number of old models to include in new run as initial conditions
        self.InitCond_override = override # Override peak identification and manually set peak times
        self.InitCond_nostatic = nostatic or onlyorbital # No static models will be calculated.
        self.InitCond_noparallax = onlyorbital; # Only orbital motion models will be calculated.
        self.InitCond_usesatellite = usesatellite; # Satellite to be used for initial conditions. Ground telescopes by default.
        
    def InitCond(self):
        if(not os.path.exists(self.eventname + '/' + self.inidir)):
            os.makedirs(self.eventname + '/' + self.inidir)
        with open(self.eventname + '/' + self.inidir + '/InitCond.ini','w') as f:
            f.write('npeaks = ' + str(self.InitCond_npeaks) + '\n')
            f.write('peakthreshold = ' + str(self.InitCond_peakthreshold) + '\n')
            f.write('oldmodels = ' + str(self.InitCond_oldmodels) + '\n')
            f.write('usesatellite = ' + str(self.InitCond_usesatellite) + '\n')
            if(self.InitCond_nostatic):            
                f.write('nostatic = 1\n')
            if(self.InitCond_noparallax):            
                f.write('noparallax = 1\n')
            if(self.InitCond_override != None):
                f.write('override = ' + str(self.InitCond_override[0])+ ' ' + str(self.InitCond_override[1]) + '\n')            
        print('- Launching: InitCond')
        print('  Setting initial conditions...')
        completedprocess = subprocess.run([self.bindir+self.initcondexe,self.eventname], cwd = self.bindir, shell = False, stdout=subprocess.DEVNULL)
        if(completedprocess.returncode != 0):
            print('! Error in setting initial conditions!')
            self.done = True
        else:
            print('  OK')

    def config_LevMar(self, nfits = 5, timelimit = 600.0, maxsteps = 50, bumperpower = 2.0):
        self.LevMar_nfits = nfits # Number of models to be calculated from the same initial condition using the bumper method
        self.LevMar_maxsteps = maxsteps # Maximum number of steps in each fit
        self.LevMar_timelimit = timelimit # Maximum time in seconds for total execution
        self.LevMar_bumperpower = bumperpower # Repulsion factor of bumpers
    
    def LevMar(self,strmodel):
        if(not os.path.exists(self.eventname + '/' + self.inidir)):
            os.makedirs(self.eventname + '/' + self.inidir)
        with open(self.eventname + '/' + self.inidir + '/LevMar.ini','w') as f:
            f.write('nfits = ' + str(self.LevMar_nfits) + '\n')
            f.write('maxsteps = ' + str(self.LevMar_maxsteps) + '\n')
            f.write('timelimit = ' + str(self.LevMar_timelimit) + '\n')
            f.write('bumperpower = ' + str(self.LevMar_bumperpower) + '\n')
        print('- Launching: LevMar')
        print('  Fitting ' + strmodel + ' ...')
        completedprocess = subprocess.run([self.bindir+self.levmarexe,self.eventname, strmodel,self.satellitedir], cwd = self.bindir, shell = False, stdout=subprocess.DEVNULL)
        if(completedprocess.returncode != 0):
            print('! Error in fit!')
            self.done = True
        else:
            print('  OK')
  
    def launch_fits(self,modelcode):
        if(not os.path.exists(self.eventname + '/' + self.inidir)):
            os.makedirs(self.eventname + '/' + self.inidir)
        with open(self.eventname + '/' + self.inidir + '/LevMar.ini','w') as f:
            f.write('nfits = ' + str(self.LevMar_nfits) + '\n')
            f.write('maxsteps = ' + str(self.LevMar_maxsteps) + '\n')
            f.write('timelimit = ' + str(self.LevMar_timelimit) + '\n')
            f.write('bumperpower = ' + str(self.LevMar_bumperpower) + '\n')
        stringfits = {'PS' : '- Single-lens-Single-source fits',
                      'PX' : '- Single-lens-Single-source fits with parallax',
                      'BS' : '- Single-lens-Binary-source fits',
                      'BO' : '- Single-lens-Binary-source fits with xallarap',
                      'LS' : '- Binary-lens-Single-source fits',
                      'LX' : '- Binary-lens-Single-source fits with parallax',
                      'LO' : '- Binary-lens-Single-source fits with orbital motion'}       
        print(stringfits[modelcode])
        initcondfile = self.eventname + '/InitCond/' + 'InitCond'+ modelcode + '.txt'
        if(os.path.exists(initcondfile)):
            with open(self.eventname + '/InitCond/' + 'InitCond'+ modelcode + '.txt') as f:
                line = f.readline().split()
                npeaks = int(line[0])
                ninitconds = int(line[1])   
            processes = []
            procnumbers = []
            iinitcond = 0
            finitcond = 0
            finitcondold = -1       
            pbar = tqdm(total = ninitconds,desc = 'Fits completed',file=sys.stdout, colour='GREEN', smoothing = 0)
            while(finitcond < ninitconds):
                i=0
                while i < len(processes):
                    if(processes[i].poll() != None):
                        processes.pop(i)
                        procnumbers.pop(i)
                        finitcond += 1
                    else:
                        i += 1
                while(iinitcond < ninitconds and len(processes) < self.nprocessors):
                    strmodel =  modelcode + '{:0>4}'.format(str(iinitcond))
                    if(glob.glob(self.eventname +'/PreModels/' + strmodel + '/t' + strmodel + '.dat')==[]):
                        processes.append(subprocess.Popen([self.bindir+self.levmarexe,self.eventname, strmodel,self.satellitedir], cwd = self.bindir, shell = False, stdout=subprocess.DEVNULL))
                        procnumbers.append(iinitcond)
                    else:
                        finitcond += 1
                    iinitcond += 1
                if(finitcond != finitcondold):
                    #print('  Fits launched: {}; completed: {}/{}'.format(iinitcond, finitcond, ninitconds))
                    pbar.update(finitcond - max(finitcondold,0))
                    finitcondold =finitcond
                time.sleep(0.1)
            pbar.close()
        else:
            print('- No initial conditions for this class')
 
    def config_ModelSelector(self, sigmasoverlap = 3.0, sigmachisquare = 1.0, maxmodels = 10):
        self.ModelSelector_sigmasoverlap = sigmasoverlap # factor multiplying the inverse covariance in search for superpositions (models are incompatible if farther than sigmasoverlap*sigma)
        self.ModelSelector_sigmachisquare = sigmachisquare # number of sigmas in the chi square distribution for accepting alternative models after the best one
        self.ModelSelector_maxmodels = maxmodels # maximum number of models returned
    
    def ModelSelector(self, modelcode):
        if(not os.path.exists(self.eventname + '/' + self.inidir)):
            os.makedirs(self.eventname + '/' + self.inidir)
        with open(self.eventname + '/' + self.inidir + '/ModelSelector.ini','w') as f:
            f.write('sigmasoverlap = ' + str(self.ModelSelector_sigmasoverlap) + '\n')
            f.write('sigmachisquare = ' + str(self.ModelSelector_sigmachisquare) + '\n')
            f.write('maxmodels = ' + str(self.ModelSelector_maxmodels) + '\n')
        stringmodels = {'PS' : '- Selecting models for Single-lens-Single-source fits',
                        'PX' : '- Selecting models for Single-lens-Single-source fits with parallax',
                        'BS' : '- Selecting models for Single-lens-Binary-source fits',
                        'BO' : '- Selecting models for Single-lens-Binary-source fits with xallarap',
                        'LS' : '- Selecting models for Binary-lens-Single-source fits',
                        'LX' : '- Selecting models for Binary-lens-Single-source fits with parallax',
                        'LO' : '- Selecting models for Binary-lens-Single-source fits with orbital motion'}
        print(stringmodels[modelcode])
        completedprocess = subprocess.run([self.bindir+self.modelselectorexe,self.eventname, modelcode], cwd = self.bindir, shell = False, stdout=subprocess.DEVNULL)
        if(completedprocess.returncode != 0):
            print('! Error in model selection!')
            self.done = True
        else:
            print('  OK')
            
    def Finalizer(self):
        print('- Launching: Finalizer')
        print('  Making final assessment for this event')
        completedprocess = subprocess.run([self.bindir+self.finalizerexe,self.eventname], cwd = self.bindir, shell = False, stdout=subprocess.DEVNULL)
        if(completedprocess.returncode != 0):
            print('! Error in finalization. Maybe there are problems with models')
            self.done = True
        else:
            with open(self.eventname + '/Nature.txt') as f:
                for line in f.readlines():
                    print("  " + line,end='')
                print("  OK")  

    def run(self, event = None):
        phase =0
        if(event!= None):
            self.eventname = os.path.realpath(event)
        self.done = False        
        while not(self.done):
            print("o " + time.asctime())
            # Check that event directory exists
            if phase == 0:
                if(os.path.exists(self.eventname + '/Data')):
                    print('- Analyzing event: ',self.eventname)
                    phase = 1
                else:
                    print('! Event data for ' + self.eventname + ' not found !')
                    self.done = True
            # Launch Reader
            elif phase == 1:
                self.Reader()
                phase = 2
            # Launch InitCond
            elif phase == 2:
                self.InitCond()
                phase = 3
            # Launch Finalizer
            elif phase == self.endphase:
                self.Finalizer()
                phase += 1
            # Conclude analysis
            elif phase > self.endphase:
                print("- Analysis of " + self.eventname + " successfully completed!")
                self.done = True
            # Launch LevMar for next class
            elif phase%2 == 1:
                self.launch_fits(self.modelcodes[phase//2-1]) 
                phase += 1            
            # Launch ModelSelector for this class
            else:
                self.ModelSelector(self.modelcodes[phase//2-2])
                phase += 1     
                
    def archive_run(self, destination = None):
        olddir = os.getcwd()
        os.chdir(self.eventname)
        if(os.path.exists('LCToFit.txt')):
            previousrunslist = glob.glob('run-*')
            previousrunslist.sort()
            if(destination == None):
                if(len(previousrunslist)>0):
                    lastrun = int(previousrunslist[-1].split('-')[-1])
                else:
                    lastrun = 0
                rundir = 'run-' + str(lastrun+1). zfill(4)
            else:
                rundir = destination
            alllist = glob.glob('*')
            alllist.remove('Data')
            filelist = list(set(alllist) - set(previousrunslist))
            os.mkdir(rundir)
            shutil.copytree('Data',rundir+'/Data')
            for nam in filelist:
                shutil.move(nam,rundir)
        os.chdir(olddir)
