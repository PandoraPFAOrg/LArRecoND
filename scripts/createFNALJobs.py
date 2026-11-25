# Script to submit LArRecoND jobs using Fermigrid

import argparse
import os
import platform
import random
import ROOT
import sys

# Initial setup parameters
class setupPars(object):

    def __init__(self, sample, minSample, nSamples, inputDir, firstEvt,
                 nEvtJob, eventTree, dataFormat):

        # Sample name
        self.sample = sample
        # First sample number, usually 0 or 1
        self.minSample = minSample
        # Number of sample input files
        self.nSamples = nSamples

        # Input file directory (these will be copied to the job scratch area)
        # relative to /exp/dune/data/users/$USER
        user = os.getenv('USER')
        self.inputDir = '/exp/dune/data/users/{0}/{1}'.format(user, inputDir)

        # First event number, usually 0 or 1
        self.firstEvt = firstEvt
        # Number of events per job
        self.nEvtJob = nEvtJob
        # Event tree name
        self.eventTree = eventTree
        # Data format
        self.dataFormat = dataFormat


# Input file parameters
class inputPars(object):

    def __init__(self, inFileName, inFileCopy, eventTree, dataFormat = 'EDepSim'):

        # Input ROOT file containing the events (-e)
        self.inFileName = inFileName
        self.inFileCopy = inFileCopy
        
        # TTree object containing the hits (-k)
        self.eventTree = eventTree
        # Data format (-f): EDepSim (default), SED, SP, SPMC
        self.dataFormat = dataFormat

        # Create copy if it doesn't exist
        if os.path.exists(inFileName) and not os.path.exists(inFileCopy):
            print('Copying {0} to {1}'.format(inFileName, inFileCopy))
            os.system('cp {0} {1}'.format(inFileName, inFileCopy))

        # Get number of events from the copy (which the jobs will use)
        self.N = 0
        if os.path.exists(inFileCopy):
            rootFile = ROOT.TFile.Open(inFileCopy, 'read')
            if rootFile.GetListOfKeys().Contains(self.eventTree):
                events = rootFile.Get(self.eventTree)
                self.N = events.GetEntries()
            rootFile.Close()


# Event numbers
class eventPars(object):

    def __init__(self, startEvt, endEvt):

        # Start, end and number of events (-s and -n)
        self.startEvt = startEvt
        self.endEvt = endEvt
        self.NEvents = self.endEvt - self.startEvt + 1


# LArRecoND run settings
class recoPars(object):

    def __init__(self, xmlFile, projection = 'both',
                 recoOption = 'AllHitsSliceNu', useHierarchy = True):

        # Xml run settings file (-i), relative to LArRecoND/settings.
        # Make sure these files don't call any monitoring algorithms, since
        # the event display can't run in batch mode and it will stop the job
        self.xmlFile = 'LArRecoND/settings/{0}'.format(xmlFile)

        # 3D, LArTPC or both projections (-j)
        self.projection = projection
        # Reco option (-r):
        # Full, AllHitsCR, AllHitsNu, CRRemHitsSliceCR,
        # CRRemHitsSliceNu, AllHitsSliceCR, AllHitsSliceNu
        self.recoOption = recoOption

        # Hierarchy tools or Validation
        self.useHierarchy = useHierarchy
        

# Geometry parameters
class geomPars(object):

    def __init__(self, physVol, sensDet, useModular = True,
                 geomFile='', geomTree=''):

        # Physical Geant4 mother volume containing LAr regions (-v)
        self.physVol = physVol
        # Sensitive detector name containing the energy hits (-d)
        self.sensDet = sensDet
        # Use modular geometry for TPC active volumes (-M)
        self.useModular = useModular
        # Geometry could be stored in a separate ROOT file (-g).
        # Location should be relative to the setupPars.inputDir directory.
        # If not set, geometry is stored in the event ROOT file
        self.geomFile = geomFile
        # TGeoManager name (-t)
        self.geomTree = geomTree


# Common job parameters
class jobPars(object):

    def __init__(self, sample, inputDir, nEvtJob, runtime = '8h', memory = '2000MB'):

        # Sample name, e.g. MiniRun4 or numu_all
        self.sample = sample

        # Directory containing the Pandora ROOT input files
        self.inputDir = inputDir

        # Number of events per job
        self.nEvtJob = nEvtJob

        # Job run time and memory
        self.runtime = runtime
        self.memory = memory

        # Pandora release (home directory) containing LArRecoND & LArContent etc.
        # This should correspond to the MY_TEST_AREA environment variable (trailing / removed)
        myTestArea = os.getenv('MY_TEST_AREA').rstrip('/')
        print('myTestArea = {0}'.format(myTestArea))

        # Set Pandora release name based on the final pandoraDir word
        pandoraName = myTestArea.split('/')[-1]

        # User scratch base directory for job input and output
        user = os.getenv('USER')
        scratchBase = '/pnfs/dune/scratch/users/{0}/{1}'.format(user, pandoraName)
        #self.scratchBase = 'testDir/{0}'.format(pandoraName)
        self.scratchSample = '{0}/{1}'.format(scratchBase, sample)

        # Create scratch area if it doesn't exist
        if not os.path.exists(self.scratchSample):
            print('Creating {0}'.format(self.scratchSample))
            os.makedirs(self.scratchSample)
            os.chmod(self.scratchSample, 0o744)

        # Tarball name to store Pandora packages. This is then used by all jobs
        self.pandoraTarName = '{0}Install.tar.gz'.format(pandoraName)

        # Required packages for LArRecoND executable (shared libraries & xml files).
        # These make up the Pandora tarball
        pandoraPkgs = ['PandoraSDK/lib', 'LArContent/lib', 'LArMachineLearningData',
                       'PandoraMonitoring/lib', 'edep-sim/install', 'LArRecoND/lib',
                       'LArRecoND/bin', 'LArRecoND/settings']

        # Create the tar file in the scratch area for the jobs
        self.pandoraTarFile = '{0}/{1}'.format(scratchBase, self.pandoraTarName)

        if not os.path.exists(self.pandoraTarFile):
            print('Creating {0}'.format(self.pandoraTarFile))
            # List of required install packages
            pgkList = ' '.join([str(pkg) for pkg in pandoraPkgs])
            tarCmd = 'tar -czf {0} -C {1} {2}'.format(self.pandoraTarFile, myTestArea, pgkList)
            print('tarCmd = {0}'.format(tarCmd))
            os.system(tarCmd)
        else:
            # Tar file has already been created
            print('Using tarball {0}'.format(self.pandoraTarFile))

        # Batch machine temporary directory location
        self.batchDir = '$_CONDOR_SCRATCH_DIR'

        # LArRecoND executable (inside tarball)
        self.recoExe = 'LArRecoND/bin/PandoraInterface'

        # Singularity container for DUNE environment
        self.singularity = '/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest'

        # Check cvmfs availability (DUNE and fifeuser), otherwise setup scripts don't work
        self.cvmfs = '\'(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)\''
        


def getInFileLabel(sample, iS):

    # Return the expected input name label, e.g. LArEDepSim_sample_1.root
    label = 'LArEDepSim_{0}_{1}.root'.format(sample, iS)

    if sample == 'MiniRun4':
        # Run number needs leading zeros (up to 5)
        runNo = '{:05d}'.format(iS)
        label = 'MiniRun4_1E19_RHC.flow.{0}.FLOWTestMergedhits.root'.format(runNo)

    return label


def createJobScript(jobScript, jobDir, ePars, gPars, iPars, jPars, rPars):

    # If jobScript exists, delete it
    if os.path.exists(jobScript):
        os.remove(jobScript)
    jobFile = open(jobScript, 'w')

    # Setup job environment
    jobFile.write('#!/bin/bash\n')
    jobFile.write('date\n')
    jobFile.write('source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh\n')
    jobFile.write('setup ifdhc\n')

    # For libxxhash.so
    jobFile.write('setup dune_oslibs v1_0_0\n')
    # Compiler
    jobFile.write('setup gcc v9_3_0\n')
    # ROOT
    jobFile.write('setup root v6_22_08d -q e20:p392:prof\n')
    # Geant4 (for edep-sim)
    jobFile.write('setup geant4 v4_10_6_p01c -q e19:prof\n')
    # LibTorch for deep learning
    jobFile.write('setup libtorch v1_6_0d -q e20\n')
    
    # Limit ifdh copying attempts to avoid stalled jobs
    jobFile.write('export IFDH_CP_MAXRETRIES=0\n')

    # Print working directory
    jobFile.write('echo Batch dir is {0}\n'.format(jPars.batchDir))
    
    # Copy Pandora release tarball from pnfs scratch area to batch machine local area
    batchTarFile = '{0}/{1}'.format(jPars.batchDir, jPars.pandoraTarName)
    jobFile.write('ifdh cp {0} {1}\n'.format(jPars.pandoraTarFile, batchTarFile))
    jobFile.write('tar -xzf {0} -C {1}\n\n'.format(batchTarFile, jPars.batchDir))

    # Specify edep-sim header file location for LArRecoND executable (for ROOT-based shared-libraries)
    edepSimHDir = '$_CONDOR_SCRATCH_DIR/edep-sim/install/include/EDepSim'
    jobFile.write('export ROOT_INCLUDE_PATH={0}:$ROOT_INCLUDE_PATH\n'.format(edepSimHDir))
    
    # LD_LIBRARY_PATH for shared libraries needed for LArRecoND executable
    libPath0 = '$_CONDOR_SCRATCH_DIR/LArRecoND/lib'
    libPath1 = '$_CONDOR_SCRATCH_DIR/PandoraSDK/lib'
    libPath2 = '$_CONDOR_SCRATCH_DIR/LArContent/lib'
    libPath3 = '$_CONDOR_SCRATCH_DIR/PandoraMonitoring/lib'
    libPath4 = '$_CONDOR_SCRATCH_DIR/edep-sim/install/lib'
    jobFile.write('export LD_LIBRARY_PATH={0}:{1}:{2}:{3}:{4}:$LD_LIBRARY_PATH\n'.format(libPath0, libPath1,
                                                                                         libPath2, libPath3,
                                                                                         libPath4))
    # Location of Pandora algorithm parameter files
    FWSPath1 = '$_CONDOR_SCRATCH_DIR/LArMachineLearningData'
    FWSPath2 = '$_CONDOR_SCRATCH_DIR/LArRecoND/settings'
    jobFile.write('export FW_SEARCH_PATH={0}:{1}\n'.format(FWSPath1, FWSPath2))

    # Copy inputFile to local batch dir
    jobFile.write('ifdh cp -D {0} $_CONDOR_SCRATCH_DIR\n'.format(iPars.inFileCopy))

    # Copy geometry file to local batch dir if required
    if gPars.geomFile != '':
        jobFile.write('ifdh cp -D {0}/{1} $_CONDOR_SCRATCH_DIR\n'.format(jPars.scratchSample, gPars.geomFile))

    # Print directory contents to check presence of directories and inputFile
    jobFile.write('ls -tral {0}\n\n'.format(jPars.batchDir))

    # LArRecoND run command: required parameters (and print events)
    runCmd = './{0} -i {1} -r {2} -e {3} -N'.format(jPars.recoExe, rPars.xmlFile, rPars.recoOption,
                                                    iPars.inFileCopy.split('/')[-1])
    # Data format, projection, event tree and numbers
    runCmd = '{0} -f {1} -j {2} -k {3} -s {4} -n {5}'.format(runCmd, iPars.dataFormat,
                                                             rPars.projection, iPars.eventTree,
                                                             ePars.startEvt, ePars.NEvents)
    # Geometry parameters
    if gPars.geomFile != '':
        runCmd = '{0} -g {1}'.format(runCmd, gPars.geomFile)
    if gPars.useModular == True:
        runCmd = '{0} -M'.format(runCmd)
    runCmd = '{0} -v {1} -d {2} -t {3}\n'.format(runCmd, gPars.physVol, gPars.sensDet,
                                                 gPars.geomTree)

    jobFile.write('cd {0}\n'.format(jPars.batchDir))
    jobFile.write(runCmd)

    # List all output to make sure everything was generated OK
    jobFile.write('ls -tral {0}\n\n'.format(jPars.batchDir))

    # Copy output ROOT files from batch node to pnfs user scratch area.
    # Note that they will not be copied if they already exist
    if rPars.useHierarchy == True:
        jobFile.write('ifdh cp -D MCHierarchy.root {0}\n'.format(jobDir))
        jobFile.write('ifdh cp -D EventHierarchy.root {0}\n'.format(jobDir))
    else:
        jobFile.write('ifdh cp -D Validation.root {0}\n'.format(jobDir))
        
    # End the job
    jobFile.write('echo Copied ROOT output to {0}\n'.format(jobDir))
    jobFile.write('date\n')
    jobFile.write('exit 0\n')

    jobFile.close()

    # Set job script permissions
    os.chmod(jobScript, 0o0644)


def createJobs(sPars, gPars, rPars):

    # Common job parameters
    jPars = jobPars(sPars.sample, sPars.inputDir, sPars.nEvtJob)

    # Script to run to submit the jobs
    runAllFileName = 'runAllJobs_{0}.sh'.format(sPars.sample)
    runAllFile = open(runAllFileName, 'w')

    # Copy geometry file to scratch area if required
    if gPars.geomFile != '':
        geomFileName = '{0}/{1}'.format(sPars.inputDir, gPars.geomFile)
        geomCopy = '{0}/{1}'.format(jPars.scratchSample, gPars.geomFile)
        if not os.path.exists(geomCopy):
            print('Copying {0} to {1}\n'.format(geomFileName, geomCopy))
            os.system('cp {0} {1}\n'.format(geomFileName, geomCopy))
        else:
            print('Copied geometry file {0} exists\n'.format(geomCopy))

    # Loop over input samples
    for iS in range(sPars.minSample, sPars.nSamples+1):

        print('Sample number {0}'.format(iS))

        # Get sample input file
        inFileLabel = getInFileLabel(sPars.sample, iS)
        inFileName = '{0}/{1}'.format(sPars.inputDir, inFileLabel)

        # Input file copy in equivalent scratch area
        inFileCopy = '{0}/{1}'.format(jPars.scratchSample, inFileName.split('/')[-1])

        # Create input parameters
        iPars = inputPars(inFileName, inFileCopy, sPars.eventTree, sPars.dataFormat)
        if iPars.N == 0:
            print('Zero events. Skipping {0}\n'.format(inFileCopy))
            continue

        print('inFileName = {0}'.format(iPars.inFileName))
        print('inFileCopy = {0}'.format(iPars.inFileCopy))

        print('Number of events in {0} = {1}'.format(iPars.inFileCopy, iPars.N))
        # Set the end event remainder for the final job
        remainder = iPars.N%sPars.nEvtJob
        print('Event remainder = {0}'.format(remainder))

        # Number of jobs required to process sample events
        nJobs = int((iPars.N/sPars.nEvtJob))
        if iPars.N%sPars.nEvtJob != 0:
            nJobs += 1
        print('Number of jobs = {0} for nEvtJobs = {1}'.format(nJobs, sPars.nEvtJob))

        # Script to run jobs for just the given sample
        runSampleName = 'runJobs_{0}_sample{1}.sh'.format(sPars.sample, iS)
        runSampleFile = open(runSampleName, 'w')

        # Loop over jobs
        for iJ in range(nJobs):
            # Start and end events
            startEvt = sPars.nEvtJob*iJ + sPars.firstEvt
            endEvt = startEvt + sPars.nEvtJob - 1
            # Set end event using remainder for final job
            if iJ == nJobs-1 and remainder > 0:
                endEvt = remainder + startEvt

            print('Job {0} startEvt = {1}, endEvt = {2}'.format(iJ, startEvt, endEvt))

            # Create job directory in scratch area
            iJob = iJ + 1
            jobDir = '{0}/sample{1}_job{2}'.format(jPars.scratchSample, iS, iJob)
            print('jobDir = {0}'.format(jobDir))
            if not os.path.exists(jobDir):
                print('Creating {0}'.format(jobDir))
                os.makedirs(jobDir)
                os.chmod(jobDir, 0o744)

            # Create the job script
            jobScript = '{0}/job.sh'.format(jobDir)

            # Set event numbers
            ePars = eventPars(startEvt, endEvt)

            # Create the job script using the various parameter objects
            createJobScript(jobScript, jobDir, ePars, gPars, iPars, jPars, rPars)

            # Set job log file
            logFile = '{0}/submitJob.log'.format(jobDir)

            # Set job submission command
            jobCmd = 'jobsub_submit -N 1 --resource-provides=usage_model=OPPORTUNISTIC ' \
                     '--expected-lifetime={0} --singularity-image={1} ' \
                     '--append_condor_requirements={2} --group=dune --memory={3} -L {4} ' \
                     'file://{5}\n\n'.format(jPars.runtime, jPars.singularity, jPars.cvmfs,
                                             jPars.memory, logFile, jobScript)
            #print('jobCmd = {0}'.format(jobCmd))            
            runSampleFile.write(jobCmd)
            runAllFile.write(jobCmd)

        print('Writing job submission script {0}'.format(runSampleName))
        runSampleFile.close()

    print('Writing job submission script {0}'.format(runAllFileName))
    runAllFile.close()


def runMiniRun4():

    print('runMiniRun4')

    # Setup parameters
    minSample = 1
    nSamples = 10
    nEvtJob = 100
    firstEvt = 0
    sPars = setupPars('MiniRun4', minSample, nSamples, 'MiniRun4', firstEvt,
                      nEvtJob, 'events', 'SPMC')

    # Geometry parameters
    gPars = geomPars('volArgonCubeDetector_PV', 'volLArActive', True,
                     'Merged2x2MINERvA_v3_withRock.root', 'Default')

    # Run parameters
    rPars = recoPars('PandoraSettings_LArRecoND_ThreeD.xml', 'both',
                     'AllHitsSliceNu', True)

    # Create the job scripts
    createJobs(sPars, gPars, rPars)


def runEDepSim():

    print('runEDepSim')

    # Setup parameters
    minSample = 1
    nSamples = 10
    nEvtJob = 100
    firstEvt = 0
    sPars = setupPars('numu_all', minSample, nSamples, 'EDepSimFiles/numu_all',
                      firstEvt, nEvtJob, 'EDepSimEvents', 'EDepSim')

    # Geometry parameters
    gPars = geomPars('volArgonCubeDetector_PV_0', 'ArgonCube', True, '', 'EDepSimGeometry')

    # Run parameters
    rPars = recoPars('PandoraSettings_LArRecoND_ThreeD.xml', 'both', 'AllHitsNu', True)

    # Create the job scripts
    createJobs(sPars, gPars, rPars)


def processArgs(parser):

    parser.add_argument('--option', default='MiniRun4', metavar='Opt',
                        help='Choose run option: MiniRun4, EDepSim')


def run():

    # Process the command line arguments. Use "python createFNALJobs.py --help" to see the full list
    parser = argparse.ArgumentParser(description='List of arguments')
    processArgs(parser)
    args = parser.parse_args()

    if args.option == 'EDepSim':
        runEDepSim()
    else:
        runMiniRun4()


if __name__ == '__main__':

    run()
