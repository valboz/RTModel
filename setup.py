import sys
import setuptools
from setuptools import setup
from distutils.ccompiler import new_compiler
import os

compiler = new_compiler()
compiler.add_include_dir('RTModel/include')
if(os.name == 'nt'):
    opts = ['/std:c++17']
    opts2 = []
else:
    opts = ['-std=c++17', '-fpermissive']
    opts2 = ['-lstdc++', '-lm']
    
objects = compiler.compile(['RTModel/lib/Reader.cpp'],extra_postargs = opts)
compiler.link_executable(objects, 'RTModel/bin/Reader',extra_postargs = opts2)
objects = compiler.compile(['RTModel/lib/InitCond.cpp'],extra_postargs = opts)
compiler.link_executable(objects, 'RTModel/bin/InitCond',extra_postargs = opts2)
objects = compiler.compile(['RTModel/lib/LevMar.cpp','RTModel/lib/LevMarFit.cpp','RTModel/lib/bumper.cpp','RTModel/lib/VBBinaryLensingLibrary.cpp'],extra_postargs = opts)
compiler.link_executable(objects, 'RTModel/bin/LevMar',extra_postargs = opts2)
objects = compiler.compile(['RTModel/lib/ModelSelector.cpp','RTModel/lib/bumper.cpp'],extra_postargs = opts)
compiler.link_executable(objects, 'RTModel/bin/ModelSelector',extra_postargs = opts2)
objects = compiler.compile(['RTModel/lib/Finalizer.cpp','RTModel/lib/bumper.cpp'],extra_postargs = opts)
compiler.link_executable(objects, 'RTModel/bin/Finalizer',extra_postargs = opts2)

setup(
    packages=['RTModel', 'RTModel.plotmodel']
)
