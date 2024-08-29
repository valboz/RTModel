import site
from pathlib import Path
import VBBinaryLensing
import matplotlib.pyplot as plt
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
    times = np.arange(tmin, tmax, tstep)
    fig, ax = plt.subplots()
    vbbl = VBBinaryLensing.VBBinaryLensing()
    vbbl.Tol = accuracy
    results = vbbl.BinaryLightCurve([logs, logq, u0, alpha, logrho, logtE, t0], times)
    mags = results[0]
    ax.plot(times, mags)
    ax.set_xlabel('t/tE')
    ax.set_ylabel('mag')

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
            elif(mags[i]>lastmin + 5*vbbl.Tol):
                derivative = +1
                lastmax = mags[i]
        else:
            if(mags[i]>mags[i-1]):
                lastmax = mags[i]
                timmax = times[i]
            elif(mags[i]<lastmax - 5*vbbl.Tol):
                derivative = -1
                lastmin = mags[i]
                peaks.append(timmax)
    print('peaks:  ' + str(peaks))
    rettemps = []
    for i in range(len(peaks)):
        for j in range(i+1,len(peaks)):
            rettemps.append(parameters[0:5] + [math.floor(peaks[i]/tstep+0.5)*tstep, math.floor(peaks[j]/tstep+0.5)*tstep])
    return rettemps