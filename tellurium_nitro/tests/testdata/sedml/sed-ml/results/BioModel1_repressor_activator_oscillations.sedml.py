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
    workingDir:     '/Users/gio/Desktop/tellurium-nitro/tellurium_nitro/tests/testdata/sedml/sed-ml'
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

workingDir = r'/Users/gio/Desktop/tellurium-nitro/tellurium_nitro/tests/testdata/sedml/sed-ml'

# --------------------------------------------------------
# Models
# --------------------------------------------------------
# Model <repressor_activator_oscillations>
repressor_activator_oscillations = te.loadSBMLModel(os.path.join(workingDir, r'../models/BioModel1_repressor_activator_oscillations.xml'))




# --------------------------------------------------------
# Tasks
# --------------------------------------------------------
# Task <task_0_0>
# not part of any DataGenerator: task_0_0
# Task <repeatedTask_0_0>
repressor_activator_oscillations.resetAll()

repeatedTask_0_0 = []
# Task: <task_0_0>
task_0_0 = [None]
repressor_activator_oscillations.setIntegrator('cvode')
if repressor_activator_oscillations.conservedMoietyAnalysis == True: repressor_activator_oscillations.conservedMoietyAnalysis = False
__range__range_0_0_common_delta_A = np.linspace(start=0.2, stop=1.0, num=4)
for __k__range_0_0_common_delta_A, __value__range_0_0_common_delta_A in enumerate(__range__range_0_0_common_delta_A):
    repressor_activator_oscillations.resetAll()
    repressor_activator_oscillations['common_delta_A'] = __value__range_0_0_common_delta_A
    repressor_activator_oscillations.timeCourseSelections = ['[A]', '[C]', '[PrmR_bound]', '[PrmA]', 'time', '[mRNA_A_]', '[PrmA_bound]', '[PrmR]', '[R]', '[mRNA_R]']
    task_0_0[0] = repressor_activator_oscillations.simulate(start=0.0, end=200.0, steps=400)

    repeatedTask_0_0.extend(task_0_0)

# --------------------------------------------------------
# DataGenerators
# --------------------------------------------------------
# DataGenerator <time_repeatedTask_0_0>
__var__t = np.column_stack([sim['time'] for sim in repeatedTask_0_0])
if len(__var__t.shape) == 1:
     __var__t.shape += (1,)
time_repeatedTask_0_0 = __var__t
# DataGenerator <dataGen_repeatedTask_0_0_mRNA_R>
__var__mRNA_R = np.column_stack([sim['[mRNA_R]'] for sim in repeatedTask_0_0])
if len(__var__mRNA_R.shape) == 1:
     __var__mRNA_R.shape += (1,)
dataGen_repeatedTask_0_0_mRNA_R = __var__mRNA_R
# DataGenerator <dataGen_repeatedTask_0_0_A>
__var__A = np.column_stack([sim['[A]'] for sim in repeatedTask_0_0])
if len(__var__A.shape) == 1:
     __var__A.shape += (1,)
dataGen_repeatedTask_0_0_A = __var__A
# DataGenerator <dataGen_repeatedTask_0_0_R>
__var__R = np.column_stack([sim['[R]'] for sim in repeatedTask_0_0])
if len(__var__R.shape) == 1:
     __var__R.shape += (1,)
dataGen_repeatedTask_0_0_R = __var__R
# DataGenerator <dataGen_repeatedTask_0_0_PrmA>
__var__PrmA = np.column_stack([sim['[PrmA]'] for sim in repeatedTask_0_0])
if len(__var__PrmA.shape) == 1:
     __var__PrmA.shape += (1,)
dataGen_repeatedTask_0_0_PrmA = __var__PrmA
# DataGenerator <dataGen_repeatedTask_0_0_PrmR>
__var__PrmR = np.column_stack([sim['[PrmR]'] for sim in repeatedTask_0_0])
if len(__var__PrmR.shape) == 1:
     __var__PrmR.shape += (1,)
dataGen_repeatedTask_0_0_PrmR = __var__PrmR
# DataGenerator <dataGen_repeatedTask_0_0_C>
__var__C = np.column_stack([sim['[C]'] for sim in repeatedTask_0_0])
if len(__var__C.shape) == 1:
     __var__C.shape += (1,)
dataGen_repeatedTask_0_0_C = __var__C
# DataGenerator <dataGen_repeatedTask_0_0_PrmA_bound>
__var__PrmA_bound = np.column_stack([sim['[PrmA_bound]'] for sim in repeatedTask_0_0])
if len(__var__PrmA_bound.shape) == 1:
     __var__PrmA_bound.shape += (1,)
dataGen_repeatedTask_0_0_PrmA_bound = __var__PrmA_bound
# DataGenerator <dataGen_repeatedTask_0_0_PrmR_bound>
__var__PrmR_bound = np.column_stack([sim['[PrmR_bound]'] for sim in repeatedTask_0_0])
if len(__var__PrmR_bound.shape) == 1:
     __var__PrmR_bound.shape += (1,)
dataGen_repeatedTask_0_0_PrmR_bound = __var__PrmR_bound
# DataGenerator <dataGen_repeatedTask_0_0_mRNA_A_>
__var__mRNA_A_ = np.column_stack([sim['[mRNA_A_]'] for sim in repeatedTask_0_0])
if len(__var__mRNA_A_.shape) == 1:
     __var__mRNA_A_.shape += (1,)
dataGen_repeatedTask_0_0_mRNA_A_ = __var__mRNA_A_

# --------------------------------------------------------
# Outputs
# --------------------------------------------------------
# Output <plot2d_scan_for_delta_A>
_stacked = False
if _stacked:
    tefig = te.getPlottingEngine().newStackedFigure(title='repressor_activator_oscillationsplots', xlabel='time_repeatedTask_0_0', ylabel='', xlim=(None, None), ylim=(None, None))
else:
    tefig = te.nextFigure(title='repressor_activator_oscillationsplots', xlabel='time_repeatedTask_0_0', ylabel='', xlim=(None, None), ylim=(None, None))

for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_0'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_mRNA_R[:,k], color='#1f77b4', tag='tag0', **extra_args)
for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_1'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_A[:,k], color='#ff7f0e', tag='tag1', **extra_args)
for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_2'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_R[:,k], color='#2ca02c', tag='tag2', **extra_args)
for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_3'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_PrmA[:,k], color='#d62728', tag='tag3', **extra_args)
for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_4'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_PrmR[:,k], color='#9467bd', tag='tag4', **extra_args)
for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_5'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_C[:,k], color='#8c564b', tag='tag5', **extra_args)
for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_6'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_PrmA_bound[:,k], color='#e377c2', tag='tag6', **extra_args)
for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_7'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_PrmR_bound[:,k], color='#7f7f7f', tag='tag7', **extra_args)
for k in range(time_repeatedTask_0_0.shape[1]):
    extra_args = {}
    if k == 0:
        extra_args['name'] = 'curve_8'
    tefig.addXYDataset(time_repeatedTask_0_0[:,k], dataGen_repeatedTask_0_0_mRNA_A_[:,k], color='#bcbd22', tag='tag8', **extra_args)
if te.tiledFigure():

    if te.tiledFigure().renderIfExhausted():

        te.clearTiledFigure()

else:

    fig = tefig.render()

####################################################################################################