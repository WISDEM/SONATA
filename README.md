[![Build Status](https://github.com/WISDEM/SONATA/workflows/CI_SONATA/badge.svg?branch=master)](https://github.com/WISDEM/SONATA/actions)

# SONATA

## Background
SONATA is a toolbox for Multidiciplinary Rotor Blade Design Environment for Structural Optimization and Aeroelastic Analysis. SONATA has originally been developed at the Helicopter Technology Institute of the Technical University of Munich (TUM), Germany. The original repository is available at https://gitlab.lrz.de/HTMWTUM/SONATA. The original work was supported by the German Federal Ministry for Economic Affairs and Energy through the German Aviation Research Program LuFo V-2 and the Austrian Research Promotion Agency through the Austrian Research Program TAKE OFF in the project VARI-SPEED.

SONATA has been adapted to wind energy applications thanks to work performed at the National Renewable Energy Laboratory ([NREL](https://www.nrel.gov)) in Colorado, USA and funded by the US Department of Energy, Wind Energy Technology Office under the Big Adaptive Rotor program. This repository is managed by Pietro Bortolotti, researcher in the systems engineering group at NREL.


## Installation
SONATA can be run on mac and linux machines. No Windows installation is supported at the moment. We make use of Anaconda, which is a commonly used package manager for python. Download and install the latest anaconda version [here](https://docs.anaconda.com/anaconda/install/)

At NREL (and possibly at other institutes), first disconnect from vpn client during installation in order to avoid remote server error when trying to retrieve URLs for installation.

First setup an anaconda environment, here named sonata-env, activate it, and add the pythonocc library (v7.4.1).
You should do this in the folder that you want to clone SONATA to.

```
conda config --add channels conda-forge
conda config --add channels tpaviot
git clone git@github.com:ptrbortolotti/SONATA.git
cd SONATA
```
If you have a mac with a newer chip, run the following:
```
CONDA_SUBDIR=osx-64 conda env create --name sonata-env -f environment.yaml
conda activate sonata-env
conda config --env --set subdir osx-64 # run with sonata-env active.
cd ..
```
Otherwise, you should be able to run:
```
conda env create --name sonata-env -f environment.yaml
conda activate sonata-env
cd ..
```

Next, download the solvers VABS (commercial, use wine to run it on mac/linux systems) or in the same conda environment compile ANBA4 (open-source)

```
git clone git@github.com:ANBA4/anba4.git # (or git clone https://github.com/ANBA4/anba4.git)
cd anba4
pip install -e .
cd ..
```

Go back to where you cloned SONATA and type:

```
cd SONATA
pip install -e .
```

Done! now check your installation trying running an example

## Usage

Navigate to examples/0_beams and run the example

```
cd examples/0_beams
python 0_SONATA_init_box_beam_HT_antisym_layup_15_6_SI_SmithChopra91.py
```

Next try running the section at 30% along the blade span of the [IEA15MW refence wind turbine](https://github.com/IEAWindTask37/IEA-15-240-RWT)
```
cd ../examples/1_IEA15MW
python 1_sonata_IEA15.py
```


## Publications:

**Feil, R., Pflumm, T., Bortolotti, P., Morandini, M.:** A cross-sectional aeroelastic analysis and structural optimization tool for slender composite structures. Composite Structures Volume 253, 1 December 2020, 112755.[[link]](https://www.sciencedirect.com/science/article/pii/S0263822320326817)

**Pflumm, T., Garre, W., Hajek, M.:** A Preprocessor for Parametric Composite Rotor Blade Cross-Sections, 44th European Rotorcraft Forum, Delft, The Netherlands, 2018  [[pdf]](docs/Pflumm,%20T.%20-%20A%20Preprocessor%20for%20Parametric%20Composite%20Rotor%20Blade%20Cross-Sections%20(2018,%20ERF).pdf) [[more…\]](https://mediatum.ub.tum.de/604993?query=Pflumm&show_id=1455385) [[BibTeX\]](https://mediatum.ub.tum.de/export/1455385/bibtex)

**Pflumm, T., Rex, W., Hajek, M.:** Propagation of Material and Manufacturing Uncertainties in Composite Helicopter Rotor Blades, 45th European Rotorcraft Forum, Warsaw, Poland, 2019 [[more…\]](https://mediatum.ub.tum.de/1520025) [BibTeX\]](https://mediatum.ub.tum.de/export/1520025/bibtex)
