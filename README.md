# Tellurium Turbo <img width="50" height="50" alt="image" src="https://github.com/user-attachments/assets/5f0099d5-27bd-4d8f-b6da-2039999c0911" />


A fork of [sys-bio/tellurium](https://github.com/sys-bio/tellurium) that adds parallel, optionally GPU-accelerated
metabolic control analysis. Everything else in tellurium is unchanged; the original README follows below.

## What's new

### Quick start: using the GPU

Call `te.gpu(True)` **once** per session. After that, every `te.controlCoefficients` call uses the GPU
automatically; you don't pass anything else.

```python
import tellurium as te, numpy as np

if __name__ == '__main__':             # needed in .py scripts (uses all CPU cores); not needed in Jupyter
    te.gpu(True)                       # once: picks NVIDIA (CUDA) or Apple (MPS), else CPU

    r = te.loada('''
      J0: $X0 -> S1; k1*X0
      J1: S1 -> S2;  k2*S1
      J2: S2 -> $X1; k3*S2
      X0 = 10; k1 = 1; k2 = 2; k3 = 0.5
    ''')

    res = te.controlCoefficients(r, {'k1': np.linspace(0.5, 3, 500)})   # 500 parameter sets in one call
    print(res['fcc'][0])               # flux control coefficients for the first parameter set

    te.gpu(False)                      # optional: back to CPU
```

Only `te.controlCoefficients` uses the GPU. Regular roadrunner calls (`r.simulate()`, `r.steadyState()`,
`r.getScaledFluxControlCoefficientMatrix()`, ...) run exactly as before.

### Full API

```python
import tellurium as te, numpy as np

te.gpu("true")        # finds an NVIDIA (CUDA) or Apple Silicon (MPS) GPU, falls back to CPU
te.getGPU()           # -> "cuda", "mps" or "cpu"

r = te.loada(model)
res = te.controlCoefficients(r, {'k1': np.linspace(0.5, 3, 1000)})
res['fcc']            # (1000, reactions, reactions)  scaled flux control coefficients
res['ccc']            # (1000, species, reactions)    scaled concentration control coefficients
res['elasticities']   # (1000, reactions, species)    scaled elasticities
res['ufcc'], res['uccc']        # unscaled versions
res['species'], res['reactions'], res['ok']   # ids, and False where no steady state was found (NaN rows)
```

- `samples` can be a dict of equal-length arrays (zipped, not a grid), a list of dicts, or a pandas DataFrame.
  `None` evaluates the model's current values once.
- Steady states and elasticities run in parallel CPU worker processes (all cores when there are 64+ samples;
  override with `processes=`). Scripts need an `if __name__ == '__main__':` guard; notebooks don't.
- Control coefficients are computed as `C = -L (N_R E L)^-1 N_R` and `F = I + E C` in one batched solve,
  on the GPU when `te.gpu` is on. Results match roadrunner's own matrices to ~1e-15 on CPU.
  Apple MPS only supports float32, so results there agree to ~1e-5.
- GPU support needs PyTorch: `pip install -e ".[gpu]"`. Without it `te.gpu` warns and stays on CPU.

## Benchmark report

Apple Silicon Mac, wall clock. Task: scaled flux + concentration control coefficients for every parameter set.

| Workload | Regular tellurium | Tellurium Turbo | Speedup |
|---|---|---|---|
| 200 parameter sets, 40 species, 1 CPU core | 321 s¹ | 105.5 s | **3.0×** |
| 200 parameter sets, 40 species, all CPU cores | 321 s¹ | 20.0 s | **16.1×** |
| 200 parameter sets, 40 species, all CPU cores + MPS GPU | 321 s¹ | 22.8 s | **14.1×** |
| Matrix algebra only, 2000 sets, 200 species, MPS GPU | 3.86 s² | 3.01 s | **1.3×** |

¹ Loop of `steadyState()` + `getScaledFluxControlCoefficientMatrix()` + `getScaledConcentrationControlCoefficientMatrix()`,
timed on 20 sets and extrapolated to 200.
² Regular tellurium has no batched equivalent; baseline is Tellurium Turbo's own CPU (numpy) path.

**Takeaway:** ~95% of the time is roadrunner's steady-state solve, which only runs on CPU, so the big win is the
algorithm (steady state + elasticities once, instead of roadrunner's per-call MCA) and using every core.
The GPU only accelerates the matrix algebra: it slightly *slows down* small models (transfer overhead) and helps
for large models with many samples. CUDA was not benchmarked.

## Changed files

- `tellurium/analysis/gpu.py` (new): `gpu()`, `getGPU()`, `controlCoefficients()`
- `tellurium/__init__.py`: exports `te.gpu`, `te.getGPU`, `te.controlCoefficients`
- `tellurium/tests/test_gpu.py` (new): checks results against roadrunner for CPU serial, CPU multi-process and GPU
- `tellurium/visualization/sbmldiagram.py`: removed the Python 2 `six` import that made `import tellurium` fail
  when `six` isn't installed
- `setup.py`: new optional extra `gpu` (`torch>=2.0`)

---

<h1><img title="tellurium logo" src="./docs/images/tellurium_logo.png" height="50" />&nbsp;&nbsp;tellurium</h1>

<!-- These badges come from shield.io, zenodo or travis -->

<table style="width:100%">
  <tr>
    <!--<td><img alt="Github version", src="https://travis-ci.com/sys-bio/tellurium.svg?branch=master"></td> -->
    <td><img alt="Read the Docs" src="https://img.shields.io/readthedocs/tellurium"></td>
    <td><a href="https://badge.fury.io/gh/sys-bio%2Ftellurium"><img src="https://badge.fury.io/gh/sys-bio%2Ftellurium.svg" alt="GitHub version" height="18"></a></td>
  </tr>
</table> 

 <table style="width:100%">
  <tr>
    <td><img alt="Licence", src="https://img.shields.io/hexpm/l/tellurium"</td>
    <td><img alt="PyPI - Downloads", src="https://img.shields.io/pypi/dm/tellurium"></td>
    <td><img alt="Funding", src="https://img.shields.io/badge/Funding-NIGMS%20(GM123032)-blue"></td>
    <td><img alt="Funding", src="https://img.shields.io/badge/Funding-NIBIB%20(EB028887)-blue"></td>
    <td><a href="https://badge.fury.io/py/tellurium"><img src="https://badge.fury.io/py/tellurium.svg" alt="PyPI version" height="18"></a> </td>
   </tr>
</table> 

 <table style="width:100%">
  <tr>
    <td><a href="https://doi.org/10.5281/zenodo.2548944"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.2548944.svg" alt="DOI"></a></td>
  </tr>
</table>

Copyright 2014-2026
Kiri Choi, J Kyle Medley, Matthias König, Kaylene Stocking, Caroline Cannistra, Michal Galdzicki, Ciaran Welsh, Lucian Smith, Adel Heydarabadipour, and Herbert Sauro

## What is it?

Tellurium is a comprehensive Python accessible SBML-compliant modeling platform for biochemical reaction systems. 

## Introduction

For installation instructions scroll down a few sections.

Tellurium is a python environment for reproducible dynamical modeling of biological networks. 
Tellurium provides the interface code to convert between standard formats and uses powerful 
libraries without requiring technical expertise, allowing teh user to focus on what’s important: 
building better models. Tellurium also provides first-class support for exchangeability via 
[COMBINE archives](http://co.mbine.org/documents/archive) and SBML, allowing you to share your models 
and simulations with many other tools.

Tellurium combines state-of-the-art scientific Python libraries, such 
as [NumPy](http://www.numpy.org/) and [SciPy](http://www.scipy.org/), 
and includes special-purpose systems biology Python tools. Out of the box, 
Tellurium includes the high performance simulation engine [libroadrunner](https://github.com/sys-bio/roadrunner), the model descrption language [antimony](http://antimony.sourceforge.net/), [phrasedml](http://phrasedml.sf.net/), [libsbml](http://sbml.org/Software/libSBML), and [libsedml](https://github.com/fbergmann/libSEDML).

The Tellurium (and libroadrunner project) project is funded from the NIH/NIGMS (GM081070) and NIH/NIBIB U24EB028887.

## Documentation 
* General: http://tellurium.readthedocs.org/en/latest/
* API: http://tellurium.readthedocs.io/en/latest/_apidoc/tellurium.html
* Webpage: http://tellurium.analogmachine.org/

## Usage Example

```{python}
import tellurium as te

r = te.loada('''
      S1 -> S2; k1*S1
      S1 = 10; S2 = 0
      k1 = 0.1
''')

result = r.simulate(0, 40, 500) 
r.plot()
```

<img src="./docs/images/tellurium-front-page-image.png" width=50% height=50%>

## Installation Instructions

Tellurium can be installed as a Python package, with a front-end (Notebook, IDE), or in a cloud environment (Google Colab):


A. [Python Pip](#Install-via-pip)

B. [IDE front-end](#Front-end-installation) based on [Spyder](https://www.spyder-ide.org/)  

C. [Notebook front-end](#front-end-2-tellurium-notebook)  

D. [Google Colab](#google-colab)  

We recommend first-time users choose one of the front-ends, while developers looking to integrate Tellurium use the pip package. All python desktop IDE graphical front-ends should work. We know for certain that telluirum will work from Spyder, pyCharm and Visual Studio Code.

The IDE front-end provides a MATLAB like experience with a code editor and Python console. 
The notebook front-end provides a notebook interface similar to [Jupyter](http://jupyter.org/), 
and features notebook cells for inline OMEX, a human-readable representation of COMBINE archives.

For any issues with installation and setup, please check [FAQ](https://github.com/sys-bio/tellurium/wiki/FAQ) or [contact us](#contact-us). 

-------

### Install via pip
[![PyPI version](https://badge.fury.io/py/tellurium.svg)](https://badge.fury.io/py/tellurium)

If you have a version of Python in an environment you're comfortable with, you can add Tellurium with pip:

```
pip install tellurium
```	

This will install on all platforms. For Colab see the end of this screen. 

For those using Anaconda, we currently do not have a separate Anaconda version of Tellurium, so you'll need to install it via pip, as above.

:exclamation:**For detailed instructions on how to setup Tellurium on Anaconda distributions, 
see [this page](https://github.com/sys-bio/tellurium/wiki/FAQ#i-would-like-to-use-tellurium-on-anaconda-what-should-i-do).** 


For developers, the latest stable version from the repository can be installed via
```
pip install git+https://github.com/sys-bio/tellurium.git
```

-------

### Front-end installation

<img align="left" width="32px" id="windows" src="https://raw.githubusercontent.com/wiki/sys-bio/tellurium/img/windows.png">
<h4>Windows</h4>
<br style="clear:both"/>

#### Approach 1: Spyder for Tellurium IDE
The Tellurium Spyder installers are tested with Windows 10 and 11, and come with Python 3.12.

:exclamation:**Attention to those upgrading to Tellurium Spyder IDE from an older version :** 
We strongly suggest you to completely remove the older version of Tellurium Spyder IDE prior to installing the latest version.

1. Download Tellurium Spyder for Windows:
    * [Python 3.12](https://sourceforge.net/projects/pytellurium/files/Tellurium-2.4/2.4.2/Tellurium-2.4.2-Python-3.12-win64-setup.exe/download)
2. Double-click the installer to start the installation
3. Follow the instructions

If you wish to use a different version of Python, you'll need to install Tellurium from the windows command line by using the command: `pip install tellurium`

NOTE: Installation requires administrative rights. It is recommended to accept the default settings.

#### Approach 2: Thonny IDE for Python
Tellurium can be used with [Thonny](https://en.wikipedia.org/wiki/Thonny), a free, open-source IDE for Python that is designed for beginners.

1. [Download the latest version of Thonny for Windows](https://thonny.org)
2. Double-click the installer `*.exe` to start the installation
3. Follow the instructions
4. Launch Thonny
5. Go to Tools > Manage packages...
6. Type `tellurium` in the search box and press the "Search on PyPI" button
7. Once the `tellurium` package is found in "Search results", click on it, and then press the "Install" button. Wait until the installation process is completed.

&nbsp;

<img align="left" width="32px" id="mac-osx" src="https://raw.githubusercontent.com/wiki/sys-bio/tellurium/img/macos.png">
<h4>Mac OS X</h4>
<br style="clear:both"/>

For Mac users there a a number of options. By far, the easiest is to use the self-cotbained Thonny Python IDE.

#### Approach 1: Thonny IDE for Python
Tellurium can be used with [Thonny](https://en.wikipedia.org/wiki/Thonny), a free, open-source IDE for Python that is designed for beginners.

1. [Download the latest version of Thonny for Mac OS X](https://thonny.org)
2. Double-click the installer `*.pkg` to start the installation
3. Follow the instructions
4. Launch Thonny from Spotlight or directly from your Applications folder
5. Go to Tools > Manage packages...
6. Type `tellurium` in the search box and press the "Search on PyPI" button
7. Once the `tellurium` package is found in "Search results", click on it, and then press the "Install" button. Wait until the installation process is completed.

#### Approach 2: Spyder for Tellurium IDE

1. [Download the Anaconda distribution for Mac OS X.](https://www.anaconda.com/download). 

2. [Follow the installer instructions.](https://docs.anaconda.com/anaconda/install/linux). The installer will install Spyder as part of the installation process. If you customize the installation, you will have a chance to select which packages to install - ensure Spyder is selected.
3. When prompted, add Anaconda to your PATH (optional, but this will make the following steps easier).
4. [Open a Terminal](https://www.wikihow.com/Open-a-Terminal-Window-in-Mac) and run the command:
```
conda install msgpack-python
```
Next, install Tellurium itself:
```
pip install tellurium
```
If this fails, try using the default absolute path to Anaconda, which will probably be something like:
```
/Users/<your-user>/opt/anaconda3/bin/pip install tellurium
```
5. Launch the Anaconda Navigator via Launchpad. In Anaconda Navigator, run Spyder and then try `import tellurium` within the Spyder editor or console.

#### Approach 3: Spyder for Tellurium IDE (only for Mac OS X Legacy)

Currently, only a Legacy version of the Mac Spyder IDE is available. 

1. [Download Tellurium IDE for Mac OS X 10.10 or later](https://github.com/sys-bio/tellurium/releases/download/1.3.5-rc3/Tellurium-1.3.5-Spyder-2.3.8-OSX.dmg)
2. Double-click the .dmg file to open a new window
3. Double-click the Spyder icon


**For those who wish to use Tellurium with Spyder IDE on Mac OS X, we also suggest that you install through [Anaconda](https://www.anaconda.com/) by 
following the instructions on [this page](https://github.com/sys-bio/tellurium/wiki/FAQ#i-would-like-to-use-tellurium-on-anaconda-what-should-i-do).**

&nbsp;

<img align="left" width="32px" id="redhat" src="https://raw.githubusercontent.com/wiki/sys-bio/tellurium/img/redhat.png">
<h4>Linux (RedHat)</h3>
<br style="clear:both"/>

<img align="left" width="32px" id="debian" src="https://raw.githubusercontent.com/wiki/sys-bio/tellurium/img/debian.png">
<h4>Linux (Debian)</h4>
<br style="clear:both"/>

**For those who wish to use Tellurium with Spyder IDE on Linux, we also suggest that you install through [Anaconda](https://www.anaconda.com/) by 
following the instructions on [this page](https://github.com/sys-bio/tellurium/wiki/FAQ#i-would-like-to-use-tellurium-on-anaconda-what-should-i-do).**

:question:Looking for old releases? [You can find them here.](https://sourceforge.net/projects/pytellurium/files/notebook) To install an old release over a new release, you will need to [manually wipe the data directory.](https://github.com/sys-bio/tellurium/wiki/FAQ#how-do-i-uninstall-tellurium-notebook)

-------

### Google Colab
Tellurium can be used entirely in a browser using [Google Colab](https://colab.research.google.com/).  Ideally, it will work with the following steps (updated May 2023)
  
1. run a cell with ```!pip install -q tellurium```
2. (Very important) Restart the runtime (From the menu: 'Runtime / Restart session')
3. Test by typing ```import tellurium as te```
4. te.getVersionInfo()  # Check for versions if you want to 
   Try a model
5. r = te.loada ('S1 -> S2; k1*S1; k1 = 0.1; S1 = 10'); r.simulate(); r.plot()
  
The Python version behind Colab changes periodically, so what worked one day may stop working the next, but the following Colab notebook worked when used at ICSB 2022:  https://colab.research.google.com/drive/1wddLftHNhetbozZY29r2HRkzQLl1F_fs#scrollTo=l1bCgW46-adR and will hopefully be instructive.

-------
  
## Citing

If you use Tellurium in your research, we would appreciate following citations in any works you publish:

Medley et al. (2018). ["Tellurium notebooks—An environment for reproducible dynamical modeling in systems biology."](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006220) *PLoS Computational Biology*, 14(6), e1006220.

Choi et al. (2018). ["Tellurium: An extensible python-based modeling environment for systems and synthetic biology."](https://www.sciencedirect.com/science/article/pii/S0303264718301254) *Biosystems*, 171, 74-79.

## Contact Us

For general questions, requesting help, or reporting bugs, feel free to use the [GitHub issue tracker](https://github.com/sys-bio/tellurium/issues).  You can also post to the  [Tellurium-discuss](https://groups.google.com/d/forum/tellurium-discuss) mailing list.

## Legal

The source code for the Tellurium Python package is hosted at https://github.com/sys-bio/tellurium and is licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0). Tellurium uses third-party dependencies which may be licensed under different terms. Consult the documentation for the respective third-party packages for more details.

TELLURIUM AND ALL SOFTWARE BUNDLED WITH TELLURIUM (HEREAFTER "THE SOFTWARE") IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THE SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
