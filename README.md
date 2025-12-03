[![run_test](https://github.com/plariono/tutorial_acts_june25/actions/workflows/run_test.yml/badge.svg)](https://github.com/plariono/tutorial_acts_june25/actions/workflows/run_test.yml)



# ACTS Tutorial for ALICE 3 studies

## Instructions how to run the full chain simulation and reconstruction.

Full ACTS documentation is located [here](https://acts.readthedocs.io/en/latest/).

# ACTS lxplus installation: 
```
lrnvmbp14@macbookpro ~ % ssh [your username]@lxplus.cern.ch
[plariono@lxplus986 plariono]$ cd /afs/cern.ch/work/[your username]
[plariono@lxplus986 plariono]$ mkdir actsdir && cd actsdir
[plariono@lxplus986 actsdir]$ git clone git@github.com:acts-project/acts.git source && cd source
```
Then the dependencies can be satisfied via an LCG release in the following way: 
`source /cvmfs/sft.cern.ch/lcg/views/<lcg_release>/<lcg_platform>/setup.sh`

```
[plariono@lxplus986 actsdir]$ source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh
[plariono@lxplus986 actsdir]$ mkdir acts
[plariono@lxplus986 actsdir]$ cmake -S source -B acts -DACTS_BUILD_PLUGIN_GEANT4=on -DACTS_BUILD_PLUGIN_ROOT=on -DACTS_BUILD_FATRAS=on -DACTS_BUILD_FATRAS_GEANT4=on -DACTS_BUILD_EXAMPLES_GEANT4=on -DACTS_BUILD_EXAMPLES_PYTHON_BINDINGS=on -DACTS_BUILD_ANALYSIS_APPS=on -DACTS_BUILD_EXAMPLES_PYTHIA8=on
[plariono@lxplus986 actsdir]$ cmake -—build acts -j$nproc
```
# ACTS aliBuild installation:
```
aliBuild init
git clone git@github.com:AliceO2Group/acts.git
aliBuild -d build ACTS
```
## Fixing potential issues with the aliBuild installation (M. Faggin)

Issue: `ImportError: libHepMC3.so.4: cannot open shared object file: No such file or directory`

Fix: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/path-to-HepMC3/3.3.0-44/lib"`

Issue: Missing or not found Boost package

Fix: Install Boost or edit the $PATH: `export PATH=$PATH:"/path-to-boost/v1.83.0-alice2-local1/include"`

When the installation is completed you can pull the tutorial repository and place it in ```actsdir```.

# Material maps:
The material-map is provided in .zip format. The full-chain job will unzip it automatically at the first run. 

## IRIS 4 geometry

https://cernbox.cern.ch/s/IcBXKkxEgaxvyLt

After that the file must be placed in ```actsdir/tutorial_acts_june25/geom/iris4```. 

## Geometry version Oct24, IRIS pure cylindrical
https://cernbox.cern.ch/s/PpjPUjx1F3HACFC

After that the file must be placed in ```actsdir/tutorial_acts_june25/geom/geom_oct24```. 

# Running the full chain simulation and reconstruction on lxplus

```
cd tutorial_acts_june25
source ../source/CI/setup_cvmfs_lcg.sh
source ../acts/python/setup.sh
python3 full_chain_alice3_tutorial.py --usePythia -n1000
```

The output files will be located in ```actsdir/reco_output_pythia```.

# Get the $\it{p}_{\rm{T}}$ resolution plot

Located in ```tutorial_acts_june25```, run:
```
root -l -q 'getPtResolution.C+g("reco_output_pythia")'
```

This will create the folders ```Plots``` and ```treeoutput``` in ```reco_output_pythia``` which contain the figure and the ROOT with the $\it{p}_{\rm{T}}$ resolution.

# Get the reconstruction efficiency 

The reconstruction efficiency histogram can be found in ```reco_output_pythia/performance_finding_ambi.root``` as a TEfficiency object:
``` trackeff_vs_pT```,``` trackeff_vs_eta```, ``` trackeff_vs_phi```.

The binning can be changed my modifying ```source/Examples/Framework/src/Validation/EffPlotTool.cpp```.
