# LArRecoND

Standalone Pandora application for developing and running DUNE ND reconstruction.

## Building Pandora with LArRecoND

The [scripts](scripts) directory contains example build and environment setup files.

The [build](scripts/build) folder contains a recipe for building LArRecoND with
all of the required [Pandora](https://github.com/PandoraPFA) packages, using the
versions defined in [tags.sh](scripts/tags.sh). This just requires 3 things to
either be installed or available (i.e. via `spack` or
[CVMFS](https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html)):
[ROOT](https://root.cern/install), [CMake](https://cmake.org/install/) and a
modern C++ compiler.
Depending on the OS, [Eigen](https://gitlab.com/libeigen/eigen) may also be
setup, if not provided by the system.

There is 3 main ways to build LArRecoND and its dependencies, depending on the
OS you are using, and where that OS is running (e.g. on a Fermilab gpvm, or on a
personal machine).

All commands below assume that you have cloned the LArRecoND repository and are
in the top-level directory of the repository (i.e. some `$MY_TEST_AREA/`
directory, that contains `LArRecoND` inside it).

### Alma9 (FNAL GPVM Native Build)

Pandora supports being built natively on Alma9, such that no container is required.
All dependencies are supplied via `spack`, so the build is relatively straight forward:

```sh
source LArRecoND/scripts/setup/Alma9_FNAL.sh
./LArRecoND/scripts/build/build_al9.sh
```

This will build a full build of Pandora, including DL support.

### SL7 (FNAL GPVM Container Build)

If instead, you need to run in an SL7 container still, you can follow these
instructions, and get your dependencies from `ups` and `cvmfs`:

```sh
source LArRecoND/scripts/setup/ContainerSL7_FNAL.sh
source LArRecoND/scripts/setup/SL7_FNAL.sh
./LArRecoND/scripts/build/build_sl7.sh
```

There also exists more specialised build scripts for SL7 inside the build folder,
to enable building with LibTorch support, and also EDepsim support.

### Personal Machine (Linux or MacOS)

Finally, if you want a local development build, you can use the following script
to build Pandora natively locally. This works by grabbing a pre-built version of
LibTorch, which may not be compatible with every single linux version. If that occurs,
you may need to alter the download link to a version of LibTorch that is
compatible with your system.

```sh
./LArRecoND/scripts/build/build.sh
```

This will build a full build of Pandora, including DL support.

### LArMachineLearningData

Various neutrino algorithms need to use MicroBooNE/SBND and DUNE training files from the
[LArMachineLearningData](https://github.com/PandoraPFA/LArMachineLearningData) package. These need to be
downloaded from CERNBox using the
[download.sh](https://github.com/PandoraPFA/LArMachineLearningData/blob/master/download.sh) script:

These are downloaded by default by the build scripts, but if you need to
download them manually, you can do so by running the following commands:

```Shell
cd $MY_TEST_AREA/LArMachineLearningData
source download.sh dunend
source download.sh dune lbl
source download.sh uboone
```

## Running LArRecoND

For each new terminal/interactive session, make sure the environment is setup by first running either the
[tags.sh](scripts/tags.sh), [Alma9_FNAL.sh](scripts/Alma9_FNAL.sh) or [SL7_FNAL.sh](scripts/SL7_FNAL.sh) scripts, where
the optional MyTestAreaDirPath parameter sets the $MY_TEST_AREA environment variable (which defaults to the current
working directory if this is not provided):

```Shell
source tags.sh MyTestAreaDirPath
```

or

```Shell
source Alma9_FNAL.sh
source tags.sh MyTestAreaDirPath
```

or

```Shell
source ContainerSL7_FNAL.sh
source SL7_FNAL.sh
source tags.sh MyTestAreaDirPath
```

The LArRecoND software is run by using the `PandoraInterface` executable, which is created from the
[PandoraInterface.cxx](test/PandoraInterface.cxx) main program that is steered with xml files from the
[settings](settings) directory. This application uses the neutrino reconstruction methods defined in the
[LArContent](https://github.com/PandoraPFA/LArContent) package, as well as the 3D algorithms listed in
[LArNDContent.cc](src/LArNDContent.cc) along with the [MasterThreeDAlgorithm.cc](src/MasterThreeDAlgorithm.cc) class.

If everything has been built correctly, running

```Shell
cd $MY_TEST_AREA/LArRecoND
./build/PandoraInterface -h
```

will list all available (required and optional) run options.

If you get runtime warnings about missing parameter files, then make sure they are downloaded in the
[LArMachineLearningData](#larmachinelearningdata) directory and their relative locations are specified
by the $FW_SEARCH_PATH environment variable, which is set by the [tags.sh](scripts/tags.sh) script.

### Converting hdf5 files from FLOW to ROOT

Scripts exist to extract the necessary components from the FLOW files, which are in hdf5 format, and store them as a ROOT
TTree read by `LArRecoND`. Running this conversion does **NOT** require one to have `LArRecoND` installed and relies on
Python and ROOT. It does however, require some non-standard Python packages, and thus the following instructions make use
of a virtual environment. First, source the relevant setup script from `LArRecoND`, e.g. assuming Alma9 and a Fermilab gpvm:

```Shell
source scripts/Alma9_FNAL.sh
```

then

```Shell
git clone https://github.com/lbl-neutrino/h5flow.git
python3 -m venv pandora.venv
source pandora.venv/bin/activate
cd h5flow
pip3 install .

pip3 install uproot
```

Then to run the conversion, this takes place in 2 steps:
- hdf5 to ROOT: make a first ROOT file from the hdf5 file, but it is in a
"subevent" format which may have more than 1 entry per event
- ROOT to ROOT: take the "subevent" level ROOT file and package it with one entry per event

**STEP 1:**
This hdf5-to-ROOT conversion works as follows:
```
python h5_to_root_ndlarflow.py FileList IsData IsFinalHits OutName
-- Parameters
FileList    [REQUIRED]:                                         comma separated set of files to convert - note it will be one output
IsData      [OPTIONAL, DEFAULT = 0, is MC]:                     1 = Data, otherwise = MC (Monte Carlo simulation)
IsFinalHits [OPTIONAL, DEFAULT = 0, prompt hits]:               1 = use "final" hits, otherwise = "prompt"
OutName     [OPTIONAL, DEFAULT = input[0]+"_hits_uproot.root"]: string for an output file name if you want to override. Note that default writes to current directory.
```

**STEP 2:**
Then the ROOT-to-ROOT conversion is a ROOT macro that can be run in compiled form. It has the following options:

```
rootToRootConversion(
  const bool isMC=true,           ---> whether this is MC [true] or data [false]
  const std::string fname,        ---> the input file name
  const std::string outname )     ---> the output file name
```

so for example

```Shell
root -l -b -q rootToRootConversion.C++\(true,\"myfile_hits_uproot.root\",\"myfile_hits.root\"\)
```

### Geometry files

The ND-LAr geometry needs to be provided as a ROOT file containing the
[TGeoManager](https://root.cern.ch/doc/master/classTGeoManager.html) object, specified by the `-g` run parameter.
These can be created from [GDML](https://gdml.web.cern.ch/GDML/) files using ROOT, where the first two calls
set the length units to cm:

```Shell
root -l
TGeoManager::LockDefaultUnits(kFALSE);
TGeoManager::SetDefaultUnits(TGeoManager::kRootUnits);
TGeoManager::Import("GeometryFile.gdml");
gGeoManager->Export("GeometryFile.root");
.q
```

The GDML files for the 2x2 ArgonCube-based prototype geometry are available in the
[2x2_sim/geometry](https://github.com/DUNE/2x2_sim/tree/develop/geometry) repository,
and the current recommended file to use is
[Merged2x2MINERvA_v4_withRock.gdml](https://github.com/DUNE/2x2_sim/blob/develop/geometry/Merged2x2MINERvA_v4/Merged2x2MINERvA_v4_withRock.gdml).


### 2x2 data

The following example can be used to run LArRecoND 3D-clustering and 2D-projection neutrino reconstruction
algorithms (without deep learning vertexing) for the first 10 events from a 2x2 data file:

```Shell
cd $MY_TEST_AREA/LArRecoND
./build/PandoraInterface -i settings/PandoraSettings_LArRecoND_ThreeD.xml \
-r AllHitsNu -e Input2x2Data.root -g Geometry2x2.root -n 10 -N
```

where the mandatory settings `-i`, `-r`, `-e` and `-g` specify the xml steering run file, reconstruction hit method,
the input data ROOT file containing the hits and the geometry ROOT file, respectively. The `-n` option sets the
number of events (in this case 10) while `-N` prints out event information.

The input hit data ROOT file uses the default [SpacePoint](include/LArSP.h) format, which needs to be previously
converted from the original HDF5 format (see relevant section above):

```Shell
python ndlarflow/h5_to_root_ndlarflow.py inputHDF5File 1 outputDir
```

where the first argument is the input HDF5 file, the second integer specifies data (1) and the final argument
is the location of the output directory which will store the equivalent ROOT file.

To use deep learning vertexing (DLVtx), make sure LArRecoND and LArContent is first built with LibTorch enabled, then use
the [PandoraSettings_LArRecoND_ThreeD_DLVtx.xml](settings/PandoraSettings_LArRecoND_ThreeD_DLVtx.xml) settings file.
You can tell if the DL vertexing is running if you see messages such as
`Loaded the TorchScript model PandoraNetworkDataFileName.pt` when the first event is getting processed.

### 2x2 simulation

The following example can be used to run LArRecoND 3D-clustering and 2D-projection neutrino reconstruction
algorithms (without deep learning vertexing) for the first 10 events from a 2x2 Monte Carlo (MC) file:

```Shell
cd $MY_TEST_AREA/LArRecoND
./build/PandoraInterface -i settings/PandoraSettings_LArRecoND_ThreeD.xml \
-r AllHitsNu -e Input2x2MC.root -g Geometry2x2.root -f SPMC -n 10 -N
```

where the mandatory settings `-i`, `-r`, `-e` and `-g` specify the xml steering run file, reconstruction hit method,
the input MC ROOT file containing the hits and the geometry ROOT file, respectively. The `-f SPMC` option sets the
input to use the [SpacePoint MC](include/LArSPMC.h) format, which stores all of the MC truth information;
this is not done by the default `-f SP` format option (the ROOT data structures are different).
The `-n` option sets the number of events (in this case 10) while `-N` prints out event information.

The input MC ROOT file needs to be previously converted from the original HDF5 format (see relevant section above):

```Shell
python ndlarflow/h5_to_root_ndlarflow.py inputHDF5File 0 outputDir
```

where the first argument is the input HDF5 file, the second integer specifies MC (0) and the final argument
is the location of the output directory which will store the equivalent ROOT file.

To use deep learning vertexing (DLVtx), make sure LArRecoND and LArContent is first built with LibTorch enabled, then use
the [PandoraSettings_LArRecoND_ThreeD_DLVtx.xml](settings/PandoraSettings_LArRecoND_ThreeD_DLVtx.xml) settings file.

### Cheating

It is possible to cheat various steps of the reconstruction, in which MC truth information is used to create "perfectly" reconstructed objects.
These can then be compared with the full reconstruction objects to see which algorithms need further improvement.
There are two cheating workflows available: fully-cheated and partially-cheated reconstruction. The fully-cheated workflow is enabled using
the [PandoraSettings_LArRecoND_ThreeD_Cheated.xml](settings/PandoraSettings_LArRecoND_ThreeD_Cheated.xml) settings file, while the
partially-cheated workflow that only cheats rock muon tracks can be run using
[PandoraSettings_LArRecoND_ThreeD_PartialCheated.xml](settings/PandoraSettings_LArRecoND_ThreeD_PartialCheated.xml).
These workflows require the `SPMC` input format, since all the MC information needs to be available; they won't work for the data `SP` format.

### edep-sim

The following example can be used to run LArRecoND 3D-clustering and 2D-projection neutrino reconstruction
algorithms (without deep learning vertexing) for the first 10 events from an edep-sim MC file:

```Shell
cd $MY_TEST_AREA/LArRecoND
./build/PandoraInterface -i settings/PandoraSettings_LArRecoND_ThreeD.xml \
-r AllHitsNu -e EDepSimMC.root -g EDepSimMC.root -f EDepSim -n 10 -N
```

where the mandatory settings `-i`, `-r`, `-e` and `-g` specify the xml steering run file, reconstruction hit method,
the input ROOT file containing the hits and the geometry ROOT file, respectively. The `-f EDepSim` option sets the
input to use the [edep-sim format](https://github.com/ClarkMcGrew/edep-sim/tree/master/io), which also stores all of the
available MC truth information. The `-n` option sets the number of events (in this case 10) while `-N` prints out
event information. Usually, the TGeoManager geometry information is stored in the event input ROOT file, so the same
filename should be used for both the `-e` and `-g` options if this is indeed the case.

To use deep learning vertexing (DLVtx), make sure LArRecoND and LArContent is first built with LibTorch enabled, then use
the [PandoraSettings_LArRecoND_ThreeD_DLVtx.xml](settings/PandoraSettings_LArRecoND_ThreeD_DLVtx.xml) settings file.

### Event displays

Pandora uses ROOT's [TEve](https://root.cern/doc/master/group__TEve.html) module for event displays in monitoring algorithms such as
[LArVisualMonitoring](https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArMonitoring/VisualMonitoringAlgorithm.cc#L364)
and [LArHierarchyMonitoring](https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArMonitoring/HierarchyMonitoringAlgorithm.cc#L479).
Calling the `LArVisualMonitoring` and/or `LArHierarchyMonitoring` algorithms at specific places in the xml run file,
for example [PandoraSettings_LArRecoND_ThreeD.xml](settings/PandoraSettings_LArRecoND_ThreeD.xml), will run the
event display at those moments in the reconstruction algorithm flow. To disable the event display (e.g. to run in
batch jobs or if there are display problems with ROOT), comment out or remove the visual monitoring calls in the
xml run file, or set the global `IsMonitoringEnabled` variable to false (which also disables the ROOT output from
the hierarchy validation tools):

```xml
    <IsMonitoringEnabled>false</IsMonitoringEnabled>
```

### Hierarchy Tools validation and analysis output

The [HierarchyAnalysisAlgorithm.cc](src/HierarchyAnalysisAlgorithm.cc) class uses
[Hierarchy Tools](https://github.com/PandoraPFA/Documentation/blob/master/Hierarchy_Tools/Hierarchy_Tools_Overview.pdf)
to create an output ROOT file that contains summary information about the reconstructed Particle Flow Objects (PFOs)
and their best-matched MC particles. The hierarchy structure and logic is implemented by LArContent's
[LArHierarchyHelper](https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArHelpers/LArHierarchyHelper.h)
class. The hierarchy analysis algorithm is enabled using the following example xml settings:

```xml
   <algorithm type = "LArHierarchyAnalysis">
        <CaloHitListName>CaloHitList2D</CaloHitListName>
        <PfoListName>RecreatedPfos</PfoListName>
        <AnalysisFileName>LArRecoND.root</AnalysisFileName>
        <AnalysisTreeName>LArRecoND</AnalysisTreeName>
        <FoldToPrimaries>false</FoldToPrimaries>
        <FoldToLeadingShowers>false</FoldToLeadingShowers>
	<FoldDynamic>true</FoldDynamic>
        <MinPurity>0.5</MinPurity>
        <MinCompleteness>0.1</MinCompleteness>
        <MinRecoHits>15</MinRecoHits>
        <MinRecoHitsPerView>5</MinRecoHitsPerView>
        <MinRecoGoodViews>2</MinRecoGoodViews>
        <RemoveRecoNeutrons>true</RemoveRecoNeutrons>
    </algorithm>

```

This creates the [TTree](https://root.cern.ch/doc/master/classTTree.html) `LArRecoND` in the output ROOT file `LArRecoND.root`
using the PFOs stored in Pandora's `RecreatedPfos` list along with the list of hits named `CaloHitList2D` (currently the
hierarchy tools can only use the 2D views). Here, the hierarchy building and matching requires minimum quality selection criteria
and removes neutrons. There are also several folding options available, which have the following (mutually exclusive)
initialization order: `FoldToPrimaries` folds the hierarchy down to the primary particles (who's direct parent is the neutrino),
else the `FoldToLeadingShowers` folds the hierarchy to the main cluster showers, otherwise the recommended option `FoldDynamic`
attempts to keep both primary and secondary interactions in the hierarchy, only folding elastic-type scatters to the relevant
parent or child particles. If none of these options are set, then `FoldDynamic` is enabled by default with all the other options
turned off. If all of these options are set to false, then no folding is done and the full hierarchy tree is kept.

The hierarchy analysis algorithm sets the event number by incrementing the number of times the `Run()` function is called
(0 to N-1 for N events). If the `-e` input file contains event numbers that are not contiguous, then the following xml
parameter settings (which must be added to the previous ones) need to be included to set the event numbers correctly,
which include settings for storing the run numbers and trigger timing, as well as associating the unique and local MC particle
IDs for CAF truth matching:

```xml
    <algorithm type = "LArHierarchyAnalysis">
        <EventFileName>LArRecoNDInput.root</EventFileName>
        <EventTreeName>events</EventTreeName>
        <EventLeafName>event</EventLeafName>
	<RunLeafName>run</RunLeafName>
	<SubRunLeafName>subrun</SubRunLeafName>
	<UnixTimeLeafName>unix_ts</UnixTimeLeafName>
	<UnixTimeUsecLeafName>unix_ts_usec</UnixTimeUsecLeafName>
	<StartTimeLeafName>event_start_t</StartTimeLeafName>
	<EndTimeLeafName>event_end_t</EndTimeLeafName>
        <MCIdLeafName>mcp_id</MCIdLeafName>
        <MCLocalIdLeafName>mcp_idLocal</MCLocalIdLeafName>
        <EventsToSkip>0</EventsToSkip>
    </algorithm>
```

Here, `EventFileName` needs to match the input file name specified by the `-e` run parameter, `EventTreeName` defines what TTree
contains the event numbers (which defaults to `events`) and `EventLeafName` defines the name of the event number variable
(which defaults to `event`). This workaround is needed since it is currently not possible to pass event (and run) information
between Pandora algorithms. By default, no events are skipped, but if the `-s` run option is used, then `EventsToSkip` must be equal
to this integer to ensure that the correct event numbers are found. The additional parameters `RunLeafName` and `SubRunLeafName`
define the run and sub-run numbers, `UnixTimeLeafName` defines the trigger unix time variable name (units in seconds), while
`StartTimeLeafName` and `EndTimeLeafName` define the decimal start and end time variables (units in ticks = 0.1 microseconds),
respectively. If these are not provided, then they will all be set to zero. Finally, `MCIdLeafName` and `MCLocalIdLeafName`
define the variable names used for the unique and local MC truth matching IDs. The unique IDs are used internally by Pandora's
Hierarchy Tools while the local ones are used by the CAF truth matching. A map is used to retrieve the local ID values given
the unique ID numbers. If these MC ID variable names are not provided, then only the unique MC IDs are used (the local IDs are
set equal to them), meaning that the CAF truth matching will be incomplete.

The xml settings files [PandoraSettings_LArRecoND_ThreeD.xml](settings/PandoraSettings_LArRecoND_ThreeD.xml) and
[PandoraSettings_LArRecoND_ThreeD_DLVtx.xml](settings/PandoraSettings_LArRecoND_ThreeD_DLVtx.xml) contain
(commented out) examples of using LArContent's
[LArHierarchyValidation](https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.h)
algorithm, which can be used to create event and particle-level ROOT output files that contain much more detail of the hierarchy than the
above LArRecoND analysis algorithm, such as complete lists of all possible reco-MC matches.
The xml files also contain (commented out) examples of using the MicroBooNE validation algorithm
[LArNeutrinoEventValidation](https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArMonitoring/NeutrinoEventValidationAlgorithm.h),
which only works for events containing single neutrino interactions (with cosmic rays).


## Fermigrid jobs

The template python script [createFNALJobs.py](scripts/createFNALJobs.py) can be used to submit LArRecoND jobs
(in SL7 containers) on the [Fermigrid](https://dune.github.io/computing-basics/07-grid-job-submission/index.html)
batch system.

It has example settings for submitting MiniRun4 or edep-sim reconstruction jobs, and it should be
relatively straightforward to extend or modify it to deal with other event samples. It uses objects
to define the setup, geometry and reconstruction parameters, which change depending on the number and
format of the sample input files. The required Pandora packages and xml steering files are stored in a tarball
by the script, which is then used by each reconstruction job, along with a copy of the input data file.
The job directories use the `/pnfs/dune/scratch/users/$USER` area, which is visible to all of the batch nodes.
Each job copies the tarball and input file to its own temporary directory area, extracts the tarball, then runs
the LArRecoND executable and copies the output files to the appropriate `/pnfs/dune/scratch/users/$USER` area.
It is recommended that the completed job output files are first copied from the scratch area to a directory in
the `/exp/dune/data/users/$USER` data area before they are used for analysis.

As an example, the following will submit reconstruction jobs for a sample of MiniRun4 input files:

```Shell
python scripts/createFNALJobs.py --option MiniRun4
source runJobs_MiniRun4.sh
```

The job run file created by the python script depends on the sample option and the number of input files.
Separate job run files are also made for each sample, which can be sourced individually to split up the
job submission process.

## GitHub Actions CI

There is a basic CI workflow for LArRecoND that is run on GitHub Actions, which
builds the software and runs a few basic tests.

Its main goals are two things: Catch build and basic runtime errors (i.e. the
software builds and reconstructs things), and to validate the reconstruction
output against reference files, to catch any changes in the reconstruction
performance. This second point is only done at low stats though, so should not
be relied on for detailed performance validation, but it is a good first check
to catch any major changes in the reconstruction output.

It does a few things:

 - It builds the software using the `build.sh` script, which builds Pandora and
 LArRecoND with all dependencies.
 - It runs a few basic tests, which include running the `PandoraInterface`
 executable with a few different settings files and input files, and checking
 that the output ROOT files are created and look vaguely valid.
 - It compares the output ROOT files against 2 refence ROOT files: one dynamically
 built on the master branch, and another for a given tagged release. This
 comparison is doing using the `validation/ciRecoValidation.py` script.
 - Finally, the result of this check, both a table of the raw max, min, average
 etc, is commented on any new PRs, and a full image of the output plots is
 uploaded to the PR for inspection.
   - The data used for the plots is also uploaded to the CI, for manual inspection if needed.


## License and Copyright

LArRecoND is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

Copyright (C), LArRecoND Authors

LArRecoND is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
