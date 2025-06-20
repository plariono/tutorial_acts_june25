# tutorial_acts_june25
ACTS Tutorial June 25

## Instructions how to run the full chain simulation and reconstruction.

Full ACTS documentation is located [here](https://acts.readthedocs.io/en/latest/).

# ACTS lxplus installation: 
```
lrnvmbp14@macbookpro ~ % ssh [your username]@lxplus.cern.ch
[plariono@lxplus986 plariono]$ cd /afs/cern.ch/work/[your username]
[plariono@lxplus986 plariono]$ mkdir actsdir && cd actsdir
[plariono@lxplus986 actsdir]$ git clone git@github.com:acts-project/acts.git source && cd source
[plariono@lxplus986 source]$ git fetch origin —-tags
[plariono@lxplus986 source]$ git switch tags/v40.1.0 -c mybranch-v40
[plariono@lxplus986 actsdir]$ source source/CI/setup_cvmfs_lcg.sh
[plariono@lxplus986 actsdir]$ mkdir acts
[plariono@lxplus986 actsdir]$ cmake -S source -B acts -DACTS_BUILD_PLUGIN_GEANT4=on -DACTS_BUILD_PLUGIN_TGEO=on -DACTS_BUILD_FATRAS=on -DACTS_BUILD_FATRAS_GEANT4=on -DACTS_BUILD_EXAMPLES_GEANT4=on -DACTS_BUILD_EXAMPLES_PYTHON_BINDINGS=on -DACTS_BUILD_ANALYSIS_APPS=on -DACTS_BUILD_EXAMPLES_PYTHIA8=on
[plariono@lxplus986 actsdir]$ cmake —build acts -j$nproc
```
# ACTS alibuild installation:
```
aliBuild init
git clone git@github.com:AliceO2Group/acts.git
aliBuild -d build ACTS
```
When the installation is completed you can pull the tutorial repository and place it in ```actsdir```.

# IMPORTANT: 
due to file size limitations the material-map.json file must be downloaded separately from:
https://cernbox.cern.ch/s/t8ZTuM5vZuehFkY

After that the file must be placed in ```actsdir/tutorial_acts_june25/geom/iris_segm_june25```. 

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