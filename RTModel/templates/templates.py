import site
from pathlib import Path
import RTModel.plotmodel as plm
import shutil
import math
import numpy as np

template_library_path = None
site_packages_directories = site.getsitepackages()
site_packages_directories.append(site.getusersitepackages())
for site_packages_directory in site_packages_directories:
    template_library_path = Path(site_packages_directory).joinpath('RTModel/data/TemplateLibrary.txt')
if template_library_path is None:
    raise FileNotFoundError(f'RTModel binary directory not found. Searched {site_packages_directories} '
                            f'site-packages directories.')

def clone_default_library(destination):
    shutil.copy(template_library_path, destination)
    return

def load_library(source = None):
    templates = None
    if(source == None):
        source = template_library_path
    with open(source) as f:
        lines = f.readlines()
        templates = []
        for i in range(1,len(lines)):
            template = [float(par) for par in lines[i].split()]
            templates.append(template)
    return templates

def save_library(destination, templates):
    with open(destination, 'w') as f:
        f.write(str(len(templates)) + '\n')
        for temp in templates:
            for num in temp:
                f.write(str(num) + ' ')
            f.write('\n')
    return

def show_template(parameters, tmin = -3, tmax = +3, tstep = 0.001, accuracy = 0.01):
    logs = math.log(parameters[0])
    logq = math.log(parameters[1])
    u0 = parameters[2]
    alpha = parameters[3]
    logrho = math.log(parameters[4])
    logtE = 0
    t0 = 0
    print('s: ' + str(parameters[0]) + '  q: ' + str(parameters[1]) + '  u0: ' + str(parameters[2]) + '  alpha: ' + str(parameters[3]) + '  rho: '+ str(parameters[4]))
    parnew = parameters[0:5] + [1,0]
    pl = plm.plotmodel(None, model = 'LS', parameters = parnew, tmin = tmin, tmax = tmax, timesteps = math.floor((tmax-tmin)/tstep+1) ,accuracy = accuracy, printpars = False)
    times = pl.t
    mags = pl.results[0]

    peaks =[]
    lastmin = mags[0]
    lastmax = -1
    timmax = -1
    derivative = -1
    i = 1
    for i in range(1,len(mags)):
        if(derivative == -1):
            if(mags[i]<mags[i-1]):
                lastmin = mags[i]
            elif(mags[i]>lastmin + 5*accuracy):
                derivative = +1
                lastmax = mags[i]
        else:
            if(mags[i]>mags[i-1]):
                lastmax = mags[i]
                timmax = times[i]
            elif(mags[i]<lastmax - 5*accuracy):
                derivative = -1
                lastmin = mags[i]
                peaks.append(timmax)
    print('peaks:  ' + str(peaks))
    rettemps = []
    for i in range(len(peaks)):
        for j in range(i+1,len(peaks)):
            rettemps.append(parameters[0:5] + [int(peaks[i]/tstep+0.5)*tstep, int(peaks[j]/tstep+0.5)*tstep])
    return rettemps
