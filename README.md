[![Build Status](https://github.com/NLRWindSystems/SONATA/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/NLRWindSystems/SONATA/actions/workflows/python-package-conda.yml)

# SONATA

## Background
SONATA is a toolbox for Multidiciplinary Rotor Blade Design Environment for Structural Optimization and Aeroelastic Analysis. SONATA has originally been developed at the Helicopter Technology Institute of the Technical University of Munich (TUM), Germany. The original repository is available at https://gitlab.lrz.de/HTMWTUM/SONATA. The original work was supported by the German Federal Ministry for Economic Affairs and Energy through the German Aviation Research Program LuFo V-2 and the Austrian Research Promotion Agency through the Austrian Research Program TAKE OFF in the project VARI-SPEED.

SONATA has been adapted to wind energy applications thanks to work performed at the National Renewable Energy Laboratory ([NREL](https://www.nrel.gov)), recently renamed National Laboratory of the Rockies (NLR), located in Golden, Colorado, USA. Work at NREL/NLR is funded by the US Department of Energy, Wind Energy Technology Office under the Big Adaptive Rotor program, and the STability and Aeroelastic Behavior of Large wind turbinEs (STABLE) project.

This repository is managed by Pietro Bortolotti and Justin Porter, researchers in the wind energy systems group at NLR.


## Part of the WETO Stack

SONATA is primarily developed with the support of the U.S. Department of Energy and is part of the [WETO Software Stack](https://nrel.github.io/WETOStack). For more information and other integrated modeling software, see:
- [Portfolio Overview](https://nrel.github.io/WETOStack/portfolio_analysis/overview.html)
- [Entry Guide](https://nrel.github.io/WETOStack/_static/entry_guide/index.html)
- [Systems Engineering Workshop](https://nrel.github.io/WETOStack/workshops/user_workshops_2024.html#systems-engineering)



## Installation
SONATA can be run on mac and linux machines. No Windows installation is supported at the moment. We make use of Anaconda, which is a commonly used package manager for Python. Download and install the latest anaconda version [here](https://docs.anaconda.com/anaconda/install/)

First setup an anaconda environment, here named sonata-env, activate it, and add the pythonocc library (v7.4.1).
You should do this in the folder that you want to clone SONATA to.

```
conda config --add channels conda-forge
git clone git@github.com:NLRWindSystems/SONATA.git
cd SONATA
git checkout wr1_b3

conda env create --name sonata-env -f environment.yaml
conda activate sonata-env
cd ..
```

Next, in the same conda environment compile ANBA4 (open-source)

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

To improve readability of git blame run (ignores commits involved with linting):
```
git config blame.ignoreRevsFile .git-blame-ignore-revs
```
If you are contributing code to SONATA, you also need to install pre-commit
```
pre-commit install
```

Done! now check your installation trying running an example

## Usage

Navigate to examples/0_beams and run the example

```
cd examples/0_beams
python 0_SONATA_init_box_beam_HT_antisym_layup_15_6_SI_SmithChopra91.py
```

Next try modeling the blade of the [IEA-15MW reference wind turbine](https://github.com/IEAWindSystems/IEA-15-240-RWT).
```
cd ../examples/1_IEA15MW
python 1_sonata_IEA15.py
```

Note that two examples use two different structures of the input files. Example 0 uses the original yaml structure from TUM, whereas the second file uses the wind turbine ontology coded in [windIO](https://github.com/IEAWindSystems/windIO). windIO is developed within IEA Wind Task 55 REFWIND and its documentation lives [here](https://ieawindsystems.github.io/windIO/main/index.html).


### Structural Damping Example

The IEA 22MW example includes instructions on how to estimate structural damping via the modal strain energy approach.
Detailed instructions can be found [here](examples/2_IEA22MW/README.md). Further details are provided in publications including Porter et al. (2025).

### Notes

1. You must define a TE anchor for SONATA. SONATA does not autodefine the TE anchor as implied by WindIO.
2. Layer definition in SONATA only supports `start_nd_arc` and `end_nd_arc`. Other tools can generate these quantities from the WindIO definition.


## Publications:

**Feil, R., Pflumm, T., Bortolotti, P., Morandini, M.:** A cross-sectional aeroelastic analysis and structural optimization tool for slender composite structures. Composite Structures Volume 253, 1 December 2020, 112755.[[link]](https://www.sciencedirect.com/science/article/pii/S0263822320326817)

**Pflumm, T., Garre, W., Hajek, M.:** A Preprocessor for Parametric Composite Rotor Blade Cross-Sections, 44th European Rotorcraft Forum, Delft, The Netherlands, 2018  [[pdf]](docs/Pflumm,%20T.%20-%20A%20Preprocessor%20for%20Parametric%20Composite%20Rotor%20Blade%20Cross-Sections%20(2018,%20ERF).pdf) [[more…\]](https://mediatum.ub.tum.de/604993?query=Pflumm&show_id=1455385) [[BibTeX\]](https://mediatum.ub.tum.de/export/1455385/bibtex)

**Pflumm, T., Rex, W., Hajek, M.:** Propagation of Material and Manufacturing Uncertainties in Composite Helicopter Rotor Blades, 45th European Rotorcraft Forum, Warsaw, Poland, 2019 [[more…\]](https://mediatum.ub.tum.de/1520025) [BibTeX\]](https://mediatum.ub.tum.de/export/1520025/bibtex)

**Porter, J. H., Mace, T., Bortolotti, P., et al.:** Prediction of structural damping in a composite structure from coupon tests, 2025, Preprint. [[link]](https://doi.org/10.2139/ssrn.5408061)
