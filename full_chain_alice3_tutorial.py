#!/usr/bin/env python3

# Author: Pavel Larionov
# Email:  pavel.larionov@cern.ch
# Event:  ACTS Tutorial, June '25

from acts.examples.reconstruction import (
    addSeeding,
    CkfConfig,
    addCKFTracks,
    addAmbiguityResolution,
    TrackSelectorConfig,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    AmbiguityResolutionConfig,
    SeedingAlgorithm,
    addAmbiguityResolutionML,
    AmbiguityResolutionMLConfig,
    # addSeedFilterML,
    # SeedFilterMLDBScanConfig,
)
from acts.examples.simulation import (
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    addParticleGun,
    addPythia8,
    addGenParticleSelection,
    addFatras,
    addGeant4,
    addDigitization,
    addDigiParticleSelection,
)

from random import *
import os
import argparse
import pathlib

import acts
import acts.examples

import array as arr

import alice3.performance.writers as alice3_writers

#import alice3_iris4_v40 as alice3_geometry
import alice3.alice3_detector_woTOF_Iris as alice3_geometry

from acts import UnitConstants as u
from acts.examples.root import RootParticleReader

parser = argparse.ArgumentParser(
    description="Full chain with the ALICE3 tracker")
parser.add_argument(
    "--output",
    "-o",
    help="Output directory",
    type=pathlib.Path,
    default=pathlib.Path.cwd().parent / "reco_output_gun",
)

parser.add_argument(
    "--events", "-n", help="Number of events", type=int, default=1)

parser.add_argument(
    "--threads", "-nthr", help="Number of threads", type=int, default=-1)

parser.add_argument(
    "--usePythia",
    help="Use Pythia8 instead of particle gun",
    action="store_true",
)

args = parser.parse_args()

outputDir = args.output if not args.usePythia else pathlib.Path.cwd().parent / \
    "reco_output_pythia"

# Partigle gun settings
particle = acts.PdgParticle.ePionPlus
npart = 1

# General settings
bFieldZ = 1
useFieldMap = True
fieldmapName = "1T_7.5m_solenoid_1m_free_bore.txt"
#
etaMin = -4
etaMax = 4
pTmin = 100  # in MeV/c
pTminCkf = 50  # in MeV/c
pTmax = 10  # in GeV/c
nHitsMin = 7
nMeasMin = 7
nMeasCutOff = 1
enableMat = True
initVarInflFactor = 1.0
enableSeedconfirmation = False

# Parameters for seeding: 0 = manually optimized, 1 = ML optimized for Pb-Pb in η = [3, 4]
seedParamOption = 0
deltaRmin = arr.array("f", [1., 1.0592])
deltaRmax = arr.array("f", [60., 61.587])
maxSeedsPerMiddleSp = arr.array("I", [1, 8])
maxSigmaScattering = arr.array("f", [5., 3.163])
radiationLengthPerSeed = arr.array("f", [0.1, 0.0785])
maxImpact = arr.array("f", [3., 22.29])
maxCotTheta = arr.array("f", [27.2899, 32.57108])

acts_dir = pathlib.Path.cwd().parent
#tutorial_dir = pathlib.Path.cwd()
tutorial_dir  = pathlib.Path(os.environ["MAIN_DIR"])

detector = alice3_geometry.buildALICE3Geometry(
    tutorial_dir, enableMat, False)

trackingGeometry = detector.trackingGeometry()
decorators = detector.contextDecorators()

#print ("PF:: Magnetic field map", tutorial_dir._str +"/"+"fieldmaps"+"/"+ fieldmapName)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, bFieldZ * u.T)) if not useFieldMap else acts.MagneticFieldMapRz(
    str(tutorial_dir / "fieldmaps" / fieldmapName))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=args.events, trackFpes=False, numThreads=args.threads, outputDir=str(outputDir))

if not args.usePythia:
    addParticleGun(
        s,
        MomentumConfig(pTmin * u.MeV, pTmax * u.GeV, transverse=True),
        EtaConfig(etaMin, etaMax, uniform=True),
        ParticleConfig(npart, particle,
                       randomizeCharge=True,
                       ),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(
                0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 0.180 * u.ns
            ),
        ),
        multiplicity=10,
        rnd=rnd   
    )
else:
    process_HI_central = ["SoftQCD:inelastic = on", "HeavyIon:bWidth=0.1", "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
                          "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0", "HeavyIon:SigFitNGen = 20",
                          "ParticleDecays:limitTau0 = on", "ParticleDecays:tau0Max = 10"]
    addPythia8(
        s,
        npileup=0,
        # beam=acts.PdgParticle.eLead,
        cmsEnergy=13.6 * acts.UnitConstants.TeV,
        pileupProcess=["SoftQCD:inelastic = on", "ParticleDecays:limitTau0 = on",
                       "ParticleDecays:tau0Max = 10", "Tune:pp = 14", "Random:setSeed = on"],
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 *
                                u.mm, 55.5 * u.mm, 5.0 * u.ns),
            mean=acts.Vector4(0, 0, 0, 0),
        ),
        rnd=rnd,
        outputDirRoot=outputDir,
    )


addGenParticleSelection(
    s,
    ParticleSelectorConfig(
        rho=(0.0, 24 * u.mm),
        absZ=(0.0, 1.0 * u.m),
        eta=(etaMin, etaMax),
        pt=(pTmin * u.MeV, None),
    ),
)

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    enableInteractions=enableMat,
    pMin=pTmin * u.MeV,
    outputDirRoot=outputDir,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=tutorial_dir / "alice3-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)

addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        pt=(pTmin * u.MeV, None),
        eta=(etaMin, etaMax),
        measurements=(nHitsMin, None),
        removeNeutral=True,
    ),
)

addSeeding(
    s,
    trackingGeometry,
    field,
    SeedFinderConfigArg(
        r=(None, 200 * u.mm),
        deltaR=(deltaRmin[seedParamOption] * u.mm,
                deltaRmax[seedParamOption] * u.mm),
        collisionRegion=(-250 * u.mm, 250 * u.mm),
        z=(-1300 * u.mm, 1300 * u.mm),
        maxSeedsPerSpM=maxSeedsPerMiddleSp[seedParamOption],
        sigmaScattering=maxSigmaScattering[seedParamOption],
        radLengthPerSeed=radiationLengthPerSeed[seedParamOption],
        minPt=pTmin * 0.85 * u.MeV,
        impactMax=maxImpact[seedParamOption] * u.mm,
        cotThetaMax=maxCotTheta[seedParamOption],
        rRangeMiddleSP=[
            [5, 60],
            [5, 100],
            [5, 100],
            [10, 130],
            [10, 130],
            [10, 130],
            [10, 130],
            [10, 130],
            [10, 130],
            [10, 130],
            [5, 100],
            [5, 100],
            [5, 60],
        ],
        seedConfirmation=enableSeedconfirmation,
        centralSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-620 * u.mm,
            zMaxSeedConf=620 * u.mm,
            rMaxSeedConf=36 * u.mm,
            nTopForLargeR=1,
            nTopForSmallR=2,
        ),
        forwardSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-1220 * u.mm,
            zMaxSeedConf=1220 * u.mm,
            rMaxSeedConf=36 * u.mm,
            nTopForLargeR=1,
            nTopForSmallR=2,
        ),
        useVariableMiddleSPRange=False,
        # deltaRMiddleMinSPRange=10 * u.mm,
        # deltaRMiddleMaxSPRange=10 * u.mm,
        # deltaRMiddleSPRange=(1 * u.mm, 10 * u.mm),
    ),
    SeedFinderOptionsArg(bFieldInZ=bFieldZ * u.T,
                         beamPos=(0 * u.mm, 0 * u.mm)),
    SeedFilterConfigArg(
        seedConfirmation=enableSeedconfirmation,
        maxSeedsPerSpMConf=5,
        maxQualitySeedsPerSpMConf=5,
    ),
    SpacePointGridConfigArg(
        zBinEdges=[
            -1300.0,
            -1100.0,
            -900.0,
            -700.0,
            -400.0,
            -250.0,
            -50.0,
            50.0,
            250.0,
            400.0,
            700.0,
            900.0,
            1100.0,
            1300.0,
        ],
        impactMax=3. * u.mm,
        phiBinDeflectionCoverage=3,
    ),
    SeedingAlgorithmConfigArg(
        # zBinNeighborsTop=[
        # [0, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 0],
        # ],
        # zBinNeighborsBottom=[
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # ],
        # numPhiNeighbors=1,
    ),
    geoSelectionConfigFile=tutorial_dir /
    "geoSelection-alice3-cfg10.json",
    outputDirRoot=outputDir,
    initialSigmas=[
        1 * u.mm,
        1 * u.mm,
        1 * u.degree,
        1 * u.degree,
        0.1 * u.e / u.GeV,
        1 * u.ns,
    ],
    initialSigmaPtRel=0.1,
    initialVarInflation=[initVarInflFactor] * 6,
    particleHypothesis=acts.ParticleHypothesis.pion,
)

alice3_writers.addCKFTracks(
    s,
    trackingGeometry,
    field,
    TrackSelectorConfig(pt=(pTminCkf * u.MeV, None),
                        nMeasurementsMin=nMeasMin,
                        maxHoles=2,
                        maxOutliers=2,
                        maxSharedHits=2),
    CkfConfig(seedDeduplication=True,
              stayOnSeed=True,
              chi2CutOffMeasurement=15.,
              chi2CutOffOutlier=25.,
              numMeasurementsCutOff=nMeasCutOff),
    outputDirRoot=outputDir,
    writeTrackSummary=True,
    writeTrackStates=False,
)


alice3_writers.addHitRemoverAlgorithm(
    s,
    "measurements",
    "ckf_tracks",
    "filter_measurements")








addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(maximumSharedHits=3, nMeasurementsMin=nMeasMin),
    outputDirRoot=outputDir,
    logLevel=acts.logging.ERROR,
    writeTrackSummary=True,
    writeTrackStates=False,
    writePerformance=True,
    writeCovMat=True,
)

s.run()
