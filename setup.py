# -*- coding: utf-8 -*-
###################################
# tellurium setup script
#
# build with 
# python -m build --wheel
###################################

from setuptools import setup
import os
try:
    with open(os.path.join(os.path.dirname(__file__), 'VERSION.txt'), 'r') as f:
        version = f.read().rstrip()
except IOError:
    # fallback location
    with open(os.path.join(os.path.dirname(__file__), 'tellurium_nitro/VERSION.txt'), 'r') as f:
        version = f.read().rstrip()

setup(name='tellurium_nitro',
      version=version,
      author='J. Kyle Medley, Kiri Choi, Matthias König, Lucian Smith, Herbert M. Sauro',
      description='Tellurium: An biological modeling environment for Python',
      url='http://tellurium.analogmachine.org/',
      long_description=open('README.md', "r").read(),
      long_description_content_type='text/markdown',
      packages=[
          'tellurium_nitro',
          'tellurium_nitro.analysis',
          'tellurium_nitro.notebooks',
          'tellurium_nitro.plotting',
          'tellurium_nitro.roadrunner',
          'tellurium_nitro.sedml',
          'tellurium_nitro.teconverters',
          'tellurium_nitro.teio',
          'tellurium_nitro.utils',
          #'tellurium_nitro.optimization',
          'tellurium_nitro.visualization',
          #'tellurium_nitro.tests',
          #'tellurium_nitro.tests.testdata',
      ],
      package_data={
          "tellurium_nitro": ["*.txt"],
          "tellurium_nitro.sedml": ["templates/*.template"],
          "tellurium_nitro.roadrunner": ["test_models/*.xml"],
      },
      install_requires=[
          # general
          'numpy>=1.23',
          'scipy>=1.5.1',
          'matplotlib>=2.0.2',
          'pandas>=2.0',
          # SBW-derived
          'libroadrunner_nitro>=0.2',  # the Rust engine and libroadrunner's Python API
          'antimony>=3.1.0',
          'phrasedml>=1.3',
          # standards
          'python-libnuml>=1.0.0',
          'python-libsedml>=2.0.17',
          'python-libcombine>=0.2.2',
          # misc
          'appdirs>=1.4.3',
          'jinja2>=3.0.0',
          'plotly>=2.0.12',
          'requests',
          # Jupyter / IPython
          'jupyter-client>=5.1.0',
          'jupyter-core>=4.3.0',
          'ipython',
          'ipykernel>=4.6.1',
          # testing
          'pytest',
          ],
      extras_require={
        # te.gpu(): CUDA or Apple MPS for te.controlCoefficients
        'gpu': [
          'torch>=2.0',
          ],
        ':(python_version < "3.10" and platform_machine != "arm64")': [
          'sbml2matlab>=0.9.1',
          ],
        ':(platform_machine != "arm64")': [
          'phrasedml>=1.0.9',
          ],
          },
      )
