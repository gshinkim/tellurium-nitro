r"""
####################################################################################################
                            tellurium 2.2.13.1
-+++++++++++++++++-         Python Environment for Modeling and Simulating Biological Systems
 .+++++++++++++++.
  .+++++++++++++.           Homepage:      http://tellurium.analogmachine.org/
-//++++++++++++/.   -:/-`   Documentation: https://tellurium.readthedocs.io/en/latest/index.html
.----:+++++++/.++  .++++/   Forum:         https://groups.google.com/forum/#!forum/tellurium-discuss
      :+++++:  .+:` .--++   Bug reports:   https://github.com/sys-bio/tellurium/issues
       -+++-    ./+:-://.   Repository:    https://github.com/sys-bio/tellurium
        .+.       `...`

SED-ML simulation experiments: http://www.sed-ml.org/
    # Change back to the original (with 'getName') when libsedml is fixed
    sedmlDoc: L1V1 
    inputType:      'SEDML_FILE'
    workingDir:     '/Users/gio/Desktop/TelluriumNitro/tellurium_nitro/tests/testdata/sedml/sed-ml'
    saveOutputs:    'False'
    outputDir:      'None'
    plottingEngine: '<MatplotlibEngine>'

macOS-26.6.2-arm64-arm-64bit
python 3.12.14 (main, Aug 12 2026, 13:57:54) [Clang 21.0.0 (clang-2100.1.1.101)]
####################################################################################################
"""
import tellurium_nitro as te
from roadrunner import Config
from tellurium_nitro.sedml.mathml import *
from tellurium_nitro.sedml.tesedml import process_trace, terminate_trace, fix_endpoints

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
try:
    import libsedml
except ImportError:
    import tesedml as libsedml
import pandas
import os.path
Config.LOADSBMLOPTIONS_RECOMPILE = True

workingDir = r'/Users/gio/Desktop/TelluriumNitro/tellurium_nitro/tests/testdata/sedml/sed-ml'

# --------------------------------------------------------
# Models
# --------------------------------------------------------
# Model <Application0>
Application0 = te.loadSBMLModel(os.path.join(workingDir, r'../models/app2sim.xml'))
# Model <Application0_0>
Application0_0 = te.loadSBMLModel(os.path.join(workingDir, r'../models/app2sim.xml'))
# /sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='s1'] 10.0
Application0_0['init([s1])'] = 10.0
if 'init([s1])' in Application0_0.getInitialAssignmentIds():
     Application0_0.removeInitialAssignment('init([s1])')
# Need to save the 'reset state' of the model to the current state.
temp_sbml = Application0_0.getCurrentSBML()
Application0_0 = te.loadSBMLModel(temp_sbml)




# --------------------------------------------------------
# Tasks
# --------------------------------------------------------
# Task <task_0_0>
Application0.resetAll()
# Task: <task_0_0>
task_0_0 = [None]
Application0.setIntegrator('cvode')
if Application0.conservedMoietyAnalysis == True: Application0.conservedMoietyAnalysis = False
Application0.timeCourseSelections = ['[s1]', 'time', '[s0]']
task_0_0[0] = Application0.simulate(start=0.0, end=20.0, steps=1000)
# Task <task_0_1>
Application0_0.resetAll()
# Task: <task_0_1>
task_0_1 = [None]
Application0_0.setIntegrator('cvode')
if Application0_0.conservedMoietyAnalysis == True: Application0_0.conservedMoietyAnalysis = False
Application0_0.timeCourseSelections = ['[s1]', 'time', '[s0]']
task_0_1[0] = Application0_0.simulate(start=0.0, end=30.0, steps=1000)

# --------------------------------------------------------
# DataGenerators
# --------------------------------------------------------
# DataGenerator <time_task_0_0>
__var__t = np.column_stack([sim['time'] for sim in task_0_0])
if len(__var__t.shape) == 1:
     __var__t.shape += (1,)
time_task_0_0 = __var__t
# DataGenerator <dataGen_task_0_0_s0>
__var__s0 = np.column_stack([sim['[s0]'] for sim in task_0_0])
if len(__var__s0.shape) == 1:
     __var__s0.shape += (1,)
dataGen_task_0_0_s0 = __var__s0
# DataGenerator <dataGen_task_0_0_s1>
__var__s1 = np.column_stack([sim['[s1]'] for sim in task_0_0])
if len(__var__s1.shape) == 1:
     __var__s1.shape += (1,)
dataGen_task_0_0_s1 = __var__s1
# DataGenerator <time_task_0_1>
__var__t = np.column_stack([sim['time'] for sim in task_0_1])
if len(__var__t.shape) == 1:
     __var__t.shape += (1,)
time_task_0_1 = __var__t
# DataGenerator <dataGen_task_0_1_s0>
__var__s0 = np.column_stack([sim['[s0]'] for sim in task_0_1])
if len(__var__s0.shape) == 1:
     __var__s0.shape += (1,)
dataGen_task_0_1_s0 = __var__s0
# DataGenerator <dataGen_task_0_1_s1>
__var__s1 = np.column_stack([sim['[s1]'] for sim in task_0_1])
if len(__var__s1.shape) == 1:
     __var__s1.shape += (1,)
dataGen_task_0_1_s1 = __var__s1

# --------------------------------------------------------
# Outputs
# --------------------------------------------------------
# Output <plot2d_Simulation0>
_stacked = False
if _stacked:
    tefig = te.getPlottingEngine().newStackedFigure(title='Application0plots', xlabel='time_task_0_0', ylabel='', xlim=(None, None), ylim=(None, None))
else:
    tefig = te.nextFigure(title='Application0plots', xlabel='time_task_0_0', ylabel='', xlim=(None, None), ylim=(None, None))

for k in range(time_task_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_0'
    tefig.addXYDataset(time_task_0_0[:,k], dataGen_task_0_0_s0[:,k], color='#1f77b4', tag='tag0', **extra_args)
for k in range(time_task_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_1'
    tefig.addXYDataset(time_task_0_0[:,k], dataGen_task_0_0_s1[:,k], color='#ff7f0e', tag='tag1', **extra_args)
if te.tiledFigure():

    if te.tiledFigure().renderIfExhausted():

        te.clearTiledFigure()

else:

    fig = tefig.render()

# Output <plot2d_Simulation1>
_stacked = False
if _stacked:
    tefig = te.getPlottingEngine().newStackedFigure(title='Application0plots', xlabel='time_task_0_1', ylabel='', xlim=(None, None), ylim=(None, None))
else:
    tefig = te.nextFigure(title='Application0plots', xlabel='time_task_0_1', ylabel='', xlim=(None, None), ylim=(None, None))

for k in range(time_task_0_1.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_0'
    tefig.addXYDataset(time_task_0_1[:,k], dataGen_task_0_1_s0[:,k], color='#1f77b4', tag='tag0', **extra_args)
for k in range(time_task_0_1.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_1'
    tefig.addXYDataset(time_task_0_1[:,k], dataGen_task_0_1_s1[:,k], color='#ff7f0e', tag='tag1', **extra_args)
if te.tiledFigure():

    if te.tiledFigure().renderIfExhausted():

        te.clearTiledFigure()

else:

    fig = tefig.render()

####################################################################################################