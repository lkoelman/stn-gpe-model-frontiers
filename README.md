# Introduction

This is the code for biological neural network model published in the paper
[Beta-Band Resonance and Intrinsic Oscillations in a Biophysically Detailed Model of the Subthalamic Nucleus-Globus Pallidus Network](https://www.frontiersin.org/articles/10.3389/fncom.2019.00077/full).

The main branch is the original code written for python 2.7. The updated code for Python 3.X is maintained on the branch [python3](https://github.com/lkoelman/stn-gpe-model-frontiers/tree/python3), with help of [@DavidCromp](https://github.com/DavidCromp).

# Installation

These installation instructions are for the Conda Python ditribution.

```sh
# BluePyOpt for various simulation tools
git clone https://github.com/BlueBrain/BluePyOpt.git
cd BluePyOpt && git checkout 1456941abe425b4a20fb084eab1cb6415ccfe2b8
pip install -e .
cd ..

# PyNN for network simulation (patched version)
git clone -b lkmn-multicomp https://github.com/lkoelman/PyNN.git
cd PyNN/pyNN/neuron/nmodl && nrnivmodl
cd ../../.. & pip install -e .
cd ..

# Neo electrophysiology data formats (patched version with MATLAB annotation support)
pip uninstall neo
git clone -b lkmn-dev https://github.com/lkoelman/python-neo.git
cd python-neo 
pip install -e .
cd ..

# Tools for LFP simulation
git clone https://github.com/lkoelman/LFPsim.git
cd LFPsim/lfpsim && nrnivmodl && cd ..
pip install -e ./LFPsim


# This network model
git clone https://github.com/lkoelman/stn-gpe-model-frontiers
cd stn-gpe-model-frontiers
python setup.py develop # or: pip install -e path_to_repo
python setup.py mechanisms # Build NMODL files automatically
```

NOTE: the custom PyNN classes in this repository have only been tested with:
- PyNN 0.9.2 (commit id `4b80343ba17afd688a92b3a0755ea73179c9daab`)
- BluePyOpt version 1.5.35 (commit id `1456941abe425b4a20fb084eab1cb6415ccfe2b8`).


# Running the model

The network model can be found in the subdirectory `bgcellmodels/models/network/LuNetStnGpe`

The commands for running the model can be found in the docstring on top, in the file `bgcellmodels/models/network/LuNetStnGpe/model_parameterized.py`. All possible command-line arguments can be found at the bottom of this file.

To run using MPI, use the `mpirun` or `mpiexec` command:

`mpirun -n 6 python model_parameterized.py --scale 0.5 --dur 3000 --dnorm --nolfp --seed 888 --transient-period 0.0 --write-interval 3000 -id CALDNORM --outdir ~/storage/LuNetStnGpe --config configs/syn-V18_template.json 2>&1 | tee CALDNORM.log`

To run from the IPython shell, use the `%run` magic function. This lets you inspect the model interactively after running and is similar to `ipython -i <script> -- <args>`:

`%run model_parameterized.py --scale 0.5 --dd --dur 100 --seed 888 --transient-period 0.0 --write-interval 1000 --no-gui -id test1 --outdir ~/storage --config configs/syn-V18_template.json`

Scripts to run the model on SLURM or PBS-based computing clusters are also provided:
- `LuNetStnGpe/batchjobs/slurm_generate_jobs.sh` will submit SLURM jobs with various parameter combinations
	- it submits the job script `LuNetStnGpe/batchjobs/slurm_runjob_netsim.sh`
- `LuNetStnGpe/batchjobs/pbs_generate_jobs.sh` will submit PBS (OpenPBS/Torque) jobs with various parameter combinations
	- it submits the job script `LuNetStnGpe/batchjobs/pbs_runjob_bgnetmodel.sh`

## BluePyOpt optimizations

```sh
# For Parallel BluePyOpt:
# In another terminal session, create ipyparallel instances in new environment
source activate neuro
cd ~/workspace/bgcellmodels
ipcluster start
```


# Install Supporting Tools (optional)

## NEURON Syntax Definitions

NEURON Syntax definitions for Hoc and NMODL languages [for Sublime Text](https://github.com/jordan-g/NEURON-for-Sublime-Text) and [for VS Code](https://github.com/imatlopez/vscode-neuron).


## nb_conda

Allows you to select a conda environment as Jupyter kernel.

```sh
conda install nb_conda
```

To enable an environment as a kernel, install the module ipykernel in it:

```sh
conda install -n my_environment ipykernel
```

## nbstripout

This facilitates working with git and Jupyter notebooks by stripping output cells before commits. This avoids committing large binary outputs like images embedded into notebooks.

```bash
pip install --upgrade nbstripout
cd bgcellmodels
nbstripout --install
```

## Jupyter extensions

```bash
pip install jupyter_contrib_nbextensions
jupyter contrib nbextension install --user
```

Now extensions can be enabled via the tab 'Nbextensions' in Jupyter. For example, 'Table of Contents' will add a handy navigation widget when working with notebooks.

--------------------------------------------------------------------------------

