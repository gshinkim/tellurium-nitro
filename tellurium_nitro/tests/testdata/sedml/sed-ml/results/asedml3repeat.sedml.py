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
    workingDir:     '/Users/gio/Desktop/tellurium-nitro/.claude/worktrees/agent-aeb060fff99a59918/tellurium_nitro/tests/testdata/sedml/sed-ml'
    saveOutputs:    'False'
    outputDir:      'None'
    plottingEngine: '<MatplotlibEngine>'

macOS-26.6.2-arm64-arm-64bit-Mach-O
python 3.14.7 (main, Aug  5 2026, 10:29:49) [Clang 21.0.0 (clang-2100.1.1.101)]
####################################################################################################
"""
import tellurium_nitro as te
from tellurium_nitro.nitro import Config
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

workingDir = r'/Users/gio/Desktop/tellurium-nitro/.claude/worktrees/agent-aeb060fff99a59918/tellurium_nitro/tests/testdata/sedml/sed-ml'

# --------------------------------------------------------
# Models
# --------------------------------------------------------
# Model <Application0>
Application0 = te.loadSBMLModel(os.path.join(workingDir, r'../models/asedml3repeat.xml'))




# --------------------------------------------------------
# Tasks
# --------------------------------------------------------
# Task <task_0_0>
# not part of any DataGenerator: task_0_0
# Task <repeatedTask_0_0>
Application0.resetAll()

repeatedTask_0_0 = []
# Task: <task_0_0>
task_0_0 = [None]
Application0.setIntegrator('cvode')
if Application0.conservedMoietyAnalysis == True: Application0.conservedMoietyAnalysis = False
__range__range_0_0_s1_init_uM = np.linspace(start=10.0, stop=30.0, num=4)
for __k__range_0_0_s1_init_uM, __value__range_0_0_s1_init_uM in enumerate(__range__range_0_0_s1_init_uM):
    Application0.resetAll()
    Application0['init([s1])'] = __value__range_0_0_s1_init_uM
    Application0.timeCourseSelections = ['[s0]', 'time', '[s1]']
    task_0_0[0] = Application0.simulate(start=0.0, end=30.0, steps=1000)

    repeatedTask_0_0.extend(task_0_0)

# --------------------------------------------------------
# DataGenerators
# --------------------------------------------------------
# DataGenerator <time_repeatedTask_0_0>
__var__t = np.column_stack([sim['time'] for sim in repeatedTask_0_0])
if len(__var__t.shape) == 1:
     __var__t.shape += (1,)
time_repeatedTask_0_0 = __var__t
# DataGenerator <dataGen_repeatedTask_0_0_s0>
__var__s0 = np.column_stack([sim['[s0]'] for sim in repeatedTask_0_0])
if len(__var__s0.shape) == 1:
     __var__s0.shape += (1,)
dataGen_repeatedTask_0_0_s0 = __var__s0
# DataGenerator <dataGen_repeatedTask_0_0_s1>
__var__s1 = np.column_stack([sim['[s1]'] for sim in repeatedTask_0_0])
if len(__var__s1.shape) == 1:
     __var__s1.shape += (1,)
dataGen_repeatedTask_0_0_s1 = __var__s1

# --------------------------------------------------------
# Outputs
# --------------------------------------------------------
# Output <plot2d_Simulation0>
_stacked = False
if _stacked:
    tefig = te.getPlottingEngine().newStackedFigure(title='Application0plots', xlabel='time_repeatedTask_0_0', ylabel='', xlim=(None, None), ylim=(None, None))
else:
    tefig = te.nextFigure(title='Application0plots', xlabel='time_repeatedTask_0_0', ylabel='', xlim=(None, None), ylim=(None, None))

for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_0'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_s0[:,k], color='#1f77b4', tag='tag0', **extra_args)
for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_1'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_s1[:,k], color='#ff7f0e', tag='tag1', **extra_args)
if te.tiledFigure():

    if te.tiledFigure().renderIfExhausted():

        te.clearTiledFigure()

else:

    fig = tefig.render()

####################################################################################################