#!/usr/bin/env python3
import os
from datetime import datetime
import shutil
import argparse
import sys
from acts.examples.reconstruction import (
    addSeeding,
    # TruthSeedRanges,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    # ParticleSmearingSigmas,
    addCKFTracks,
    addTruthTrackingGsf,
    #    CKFPerformanceConfig,
    TrackSelectorConfig,
    addKalmanTracks,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    CkfConfig,
    addVertexFitting,
    VertexFinder,
    addSpacePointsMaking,  # May 2025
    addHoughVertexFinding,  # May 2025
)
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addPythia8,
    addGenParticleSelection,
    addSimParticleSelection,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
    addDigiParticleSelection,
)
import pathlib
import acts
import acts.examples
#import alice3
import acts.examples.geant4

# bit ugly but should make it easier to switch geometries
from geom import buildALICE3Geometry

## Configurations
from job_configs import ChainConfig
cfg = ChainConfig.Config()

u = acts.UnitConstants

geo_dir = pathlib.Path(os.path.abspath(os.path.dirname(__file__))) / "geom"


def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(
        description="Command line arguments for full chain setup"
    )

    ##### General parameters
    parser.add_argument(
        "--nThreads", "-nthr", help="Number of threads", type=int, default=-1
    )
    parser.add_argument(
        "--nEv",
        "-n",
        dest="nEvents",
        help="nEvents",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--out_dir_prefix",
        dest="out_dir_prefix",
        help="Dir for output",
        type=str,
        default="test",
    )

    ##### REC: Seeding params
    parser.add_argument(
        "--seedingLayers",
        type=str,
        choices=["VD", "ML3", "MLall", "VDML"],
        required=False,
        default="VD",
        help="Choose seeding layers",
    )
    parser.add_argument(
        "--seedingAlgo",
        type=str,
        choices=["Default", "GridTriplet"],
        required=False,
        default="GridTriplet",
        help="Choose seeding algo",
    )

    parser.add_argument(
        "--minSeedPt",
        dest="minSeedPt",
        help="min pt for seed candidates, GeV/c",
        type=float,
        default=0.08,
    )

    parser.add_argument(
        "--impParForSeeds",
        dest="impParForSeeds",
        help="max impact parameter for track seeds, mm",
        type=float,
        default=1.0,  # mm
    )
    parser.add_argument(
        "--radLengthPerSeed",
        dest="radLengthPerSeed",
        help="Average Radiation Length per seed",
        type=float,
        default=0.05,  # BEST TUNED WITH PARTICLE GUN - May 16. 2025 - if combined with sigmaScattering = 5
    )
    parser.add_argument(
        "--sigmaScattering",
        dest="sigmaScattering",
        help="How many sigmas of scattering to include in seeds",
        type=float,
        default=5.0,  # BEST TUNED WITH PARTICLE GUN - May 16. 2025 - if combined with radLengthPerSeed = 0.05
    )

    parser.add_argument(
        "--maxSeedsPerSpM",
        dest="maxSeedsPerSpM",
        help="max number of seeds per seed middle spacepoint",
        type=int,
        default=2,  # ok for pp, for PbPb should be 3 or more?
    )

    ##### REC: tracking params
    parser.add_argument(
        "--nMeasurementsMin",
        dest="nMeasurementsMin",
        help="nMeasurementsMin",
        type=int,
        default=7,
    )
    parser.add_argument(
        "--ckfMeasPerSurf",
        dest="ckfMeasPerSurf",
        help="ckf maxNumAssocMeasOnSurface",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--ckfChi2Measurement",
        dest="ckfChi2Measurement",
        help="ckfChi2Measurement",
        type=float,
        default=45,
    )
    parser.add_argument(
        "--ckfChi2Outlier",
        dest="ckfChi2Outlier",
        help="ckfChi2Outlier",
        type=float,
        default=100,
    )

    parser.add_argument(
        "--maxSharedHits",
        dest="maxSharedHits",
        help="max number of shared hits. Doesn't have impact if the hit-smearing digi is used.",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--twoWayCKF",
        help="Set twoWayCKF option for CKF",
        action="store_true",  # default is false
        default=True,
    )
    parser.add_argument(
        "--seedDeduplication",
        help="Set seedDeduplication option for CKF",
        action="store_true",  # default is false
        default=True,
    )
    parser.add_argument(
        "--stayOnSeed",
        help="Set stayOnSeed option for CKF",
        action="store_true",  # default is false
        default=False,
    )

    ##### Special flags for tracking in several iterations
    parser.add_argument(
        "--iterationId",
        dest="iterationId",
        help="iterationId",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--isUsedHitsRemoved",
        help="Use hits.root, where the used hits are removed",
        action="store_true",  # default: false
    )
    parser.add_argument("--seed", "--randomSeed",
                        help="seed to use in the random number generator",
                        type=int,
                        default=42)

    parser.add_argument("--digi_configuration", "--digi",
                        help="Configuration of the digitization",
                        default="digi-smearing-config_no_TOFs_iTOF_removed.json")
    #####
    return parser


pars = getArgumentParser().parse_args()


##########################
##### SOME OTHER PARAMS
##########################

IA_collisionRegion_forSeeds = (
    250 if pars.impParForSeeds < 2.0 else 1000
)  # mm; large values - for V0 daughter reconstruction

IA_whichSeedingAlg = (
    SeedingAlgorithm.GridTriplet
)  # Nov 2025: is the default in addSeeding!
if pars.seedingAlgo == "Default":
    IA_whichSeedingAlg = SeedingAlgorithm.Default

IA_maxSeedsPerSpMConf = 1
IA_maxQualitySeedsPerSpMConf = 1
IA_enableMaterial = True

### output directory
IA_outputDirName = (
    "output/"
    + pars.out_dir_prefix
    + "_nEv"
    + str(pars.nEvents)
    + "_PID"
    + str(cfg.particleGun.gunPID)
    + "_nMeasMin"
    + str(pars.nMeasurementsMin)
    + "_ckfChi2Meas"
    + str(pars.ckfChi2Measurement)
    + "_ckfMeasPerSurf"
    + str(pars.ckfMeasPerSurf)
    # + "_B"
    # + str(pars.MF)
    # + "_PU"
    # + str(pars.pileup)
    # + "_iterationId"
    # + str(pars.iterationId)
)
# +"_twoWayCKF_"+str(pars.twoWayCKF) +"_seedDeduplication_"+str(pars.seedDeduplication)  +"_stayOnSeed_"+str(pars.stayOnSeed) \


outputDir = pathlib.Path.cwd() / IA_outputDirName


if not outputDir.exists():
    outputDir.mkdir(mode=0o777, parents=True, exist_ok=True)

detector = buildALICE3Geometry.buildALICE3Geometry(
    geo_dir, IA_enableMaterial, False, acts.logging.INFO
)
trackingGeometry = detector.trackingGeometry()
decorators = detector.contextDecorators()

if cfg.general.fieldMap != "":
    print(">>> !", cfg.general.fieldMap, "field map is used in this sim+rec run.")
    with open(cfg.general.fieldMap) as magFile:  # Checking if it's RZ or XYZ coordinates
        for l in magFile:
            l = l.strip()
            if l.startswith("#"):
                continue
            while "  " in l:
                l = l.replace("  ", " ")
            l = l.split(" ")
            if len(l) == 4:
                field = acts.MagneticFieldMapRz(cfg.general.fieldMap)
            else:
                field = acts.MagneticFieldMapXyz(cfg.general.fieldMap)
            break
else:
    print(">>> !Using Constant B-Field")
    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, cfg.general.MF * u.T))


rnd = acts.examples.RandomNumbers(seed=pars.seed)


# s = acts.examples.Sequencer(events=pars.nEvents, numThreads=-1)
s = acts.examples.Sequencer(events=pars.nEvents, numThreads=pars.nThreads)


# if True: #pars.usedHitsRemoval==0:
if not pars.isUsedHitsRemoved:
    if not cfg.pythia.usePythia:
        addParticleGun(
            s,
            MomentumConfig(
                cfg.particleGun.gunPtRange[0] * u.GeV, cfg.particleGun.gunPtRange[1] * u.GeV, transverse=True
            ),
            EtaConfig(cfg.particleGun.gunEtaRange[0], cfg.particleGun.gunEtaRange[1], uniform=True),
            PhiConfig(cfg.particleGun.gunPhiRange[0], cfg.particleGun.gunPhiRange[1]),
            ParticleConfig(
                cfg.particleGun.gunMult,
                acts.PdgParticle(cfg.particleGun.gunPID),
                randomizeCharge=cfg.particleGun.gunRandCharge,
            ),
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    0.0001 * u.mm, 0.0001 * u.mm, 0.0001 * u.mm, 0.0001 * u.ns
                ),
                mean=acts.Vector4(0, 0, cfg.particleGun.meanVz * u.mm, 0),
            ),
            rnd=rnd,
            logLevel=acts.logging.INFO,
            outputDirRoot=outputDir,
        )

    else:
        if cfg.pythia.system == "pp":
            s = addPythia8(
                s,
                npileup=pars.pileup,
                # for pp:
                beam=acts.PdgParticle.eProton,  # eLead,
                cmsEnergy=13.6 * acts.UnitConstants.TeV,  # 5 * acts.UnitConstants.TeV,
                hardProcess=[
                    "SoftQCD:inelastic = on",
                    "Tune:pp = 14",
                    "ParticleDecays:limitTau0 = on",
                    "ParticleDecays:tau0Max = 10",
                ],  # tune 14 is Monash
                vtxGen=acts.examples.GaussianVertexGenerator(
                    #    stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, pars.sigmaVz * u.mm, 5.0 * u.ns),
                    stddev=acts.Vector4(
                        0.000125 * u.mm,
                        0.000125 * u.mm,
                        cfg.pythia.sigmaVz * u.mm,
                        0.0001 * u.ns,
                    ),  # 5.0 * u.ns),
                    # stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 0.1 * u.mm, 5.0 * u.ns),
                    mean=acts.Vector4(0, 0, cfg.pythia.meanVz * u.mm, 0),
                    # mean=acts.Vector4(0, 0, 100 * u.mm, 0),
                ),
                rnd=rnd,
                logLevel=acts.logging.INFO,
                outputDirRoot=outputDir,
            )
        elif cfg.pythia.system == "PbPbMB":
            s = addPythia8(
                s,
                npileup=cfg.pythia.pileup,
                # for Pb-Pb:
                beam=acts.PdgParticle.eLead,
                cmsEnergy=5.36 * acts.UnitConstants.TeV,
                hardProcess=[
                    "SoftQCD:inelastic = on",
                    "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
                    "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
                    "HeavyIon:SigFitNGen = 20",
                    "ParticleDecays:limitTau0 = on",
                    "ParticleDecays:tau0Max = 10",
                ],
                vtxGen=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(
                        0.000125 * u.mm,
                        0.000125 * u.mm,
                        cfg.pythia.sigmaVz * u.mm,
                        0.0001 * u.ns,
                    ),  # 5.0 * u.ns),
                    mean=acts.Vector4(0, 0, cfg.pythia.meanVz * u.mm, 0),
                ),
                rnd=rnd,
                logLevel=acts.logging.INFO,
                outputDirRoot=outputDir,
            )
        elif cfg.pythia.system == "PbPbCentral":
            s = addPythia8(
                s,
                npileup=pars.pileup,
                # for Pb-Pb:
                beam=acts.PdgParticle.eLead,
                cmsEnergy=5.36 * acts.UnitConstants.TeV,
                hardProcess=[
                    "SoftQCD:inelastic = on",
                    "HeavyIon:bWidth=0.1",
                    "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
                    "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
                    "HeavyIon:SigFitNGen = 20",
                    "ParticleDecays:limitTau0 = on",
                    "ParticleDecays:tau0Max = 10",
                ],
                vtxGen=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(
                        0.000125 * u.mm,
                        0.000125 * u.mm,
                        cfg.pythia.sigmaVz * u.mm,
                        0.0001 * u.ns,
                    ),
                    mean=acts.Vector4(0, 0, cfg.pythia.meanVz * u.mm, 0),
                ),
                rnd=rnd,
                logLevel=acts.logging.INFO,
                outputDirRoot=outputDir,
            )

else:  # use already generated particles
    s.addReader(
        acts.examples.RootParticleReader(
            level=acts.logging.INFO,
            # filePath="./output/testIterations_iterId"+str(pars.iterationId)+"/hits.root",
            # filePath=IA_outputDirName + "/particles.root",
            # filePath="/Users/igor/ALICE3/geometry_2025_12_17_test_Geant4_WORKS_FOR_TRK_FT3/output/TEST_GEANT_Gun_Test_pion_widePt_nEv50000_B2.0_PU1_iterationId0/particles.root",
            # filePath="/Users/igor/ALICE3/geometry_2025_12_17_test_Geant4_WORKS_FOR_TRK_FT3/output/TEST_GEANT_Gun_Test_muon_widePt_GEOM_V3_THIN_LAYERS_nEv200000_PID13_nMeasMin7_ckfChi2Meas45.0_ckfMeasPerSurf1_BEFORE_CLUSTERS/particles.root",
            filePath="/Users/igor/ALICE3/geometry_2025_12_17_test_Geant4_WORKS_FOR_TRK_FT3/output/TEST_GEANT_PYTHIA_V3_nEv5000_PID211_nMeasMin7_ckfChi2Meas45.0_ckfMeasPerSurf1_BEFORE_CLUSTERS/particles.root",
            outputParticles="particles_generated",  # _generated",
            # outputSimHits="simhits",
            # particleCollection="particles",
            # inputDir="output",
            # inputFile="pythia8_particles.root",
        )
    )
    s.addWhiteboardAlias("particles", "particles_generated")

addGenParticleSelection(
    s,
    ParticleSelectorConfig(
        eta=(-cfg.detSim.particleSelectionEta, cfg.detSim.particleSelectionEta),
        pt=(0.001 * u.MeV, None),
        removeNeutral=False,
        # rho=(0.0, 24 * u.mm),
        # absZ=(0.0, 1.0 * u.m),
    ),
)

alice3_gdml = geo_dir / "o2sim_geometry.gdml"

gdml_detector = acts.examples.geant4.GdmlDetector(path=str(alice3_gdml))

if cfg.detSim.simulation == "Fatras":

    addFatras(
        s,
        trackingGeometry,
        field,
        enableInteractions=True,
        rnd=rnd,
        pMin=cfg.detSim.minFatrasPt,  # GeV, May 2025
        outputDirRoot=outputDir,
        logLevel=acts.logging.INFO,
    )
    
elif cfg.detSim.simulation == "Geant4":
    addGeant4(
        s,
        gdml_detector,
        trackingGeometry,
        field,
        # materialMappings=["TRK_SI", "Silicon", "Silicon_elm"],
        materialMappings=["TRK_SILICON", "FT3_SILICON", "FT3_Silicon", "TF3_Silicon"],
        # materialMappings=["TRK_SILICON_SENSITIVE"],
        outputDirRoot=outputDir,
        rnd=rnd,
        logLevel=acts.logging.INFO,
        # logLevel=acts.logging.DEBUG,
        # logLevel=acts.logging.VERBOSE,
        killVolume=trackingGeometry.highestTrackingVolume,
        killAfterTime=40 * u.ns,  # 25 * u.ns,
        recordHitsOfSecondaries=True,  # False,#True, ### IA
        keepParticlesWithoutHits=True,  # True, ### IA
        killSecondaries=False,
    )
    
# addSimParticleSelection(
#     s,
#     ParticleSelectorConfig(
#         # eta=(-4.1,4.1), pt=(1 * u.MeV, None), removeNeutral=True, # May 2025
#         hits=(7, None),  # Nov 2025
#         # rho=(0.0, 24 * u.mm),
#         # absZ=(0.0, 1.0 * u.m),
#     ),
# )


s = addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "../digiConfigurations" / pars.digi_configuration,
    outputDirRoot=outputDir,
    rnd=rnd,
    logLevel=acts.logging.INFO,
    # doMerge=True,   ##!! 
)

# for details on segmented digi: class Channelizer {..} https://github.com/acts-project/acts/blob/v36.3.2/Fatras/include/ActsFatras/Digitization/Channelizer.hpp#L21
# auto channelsRes = m_channelizer.channelize() - https://github.com/acts-project/acts/blob/v36.2.1/Examples/Algorithms/Digitization/src/DigitizationAlgorithm.cpp#L212


# addDigiParticleSelection(
#     s,
#     ParticleSelectorConfig(
#         eta=(-4.1,4.1), pt=(1 * u.MeV, None), removeNeutral=True, # May 2025
#         hits=(9, None),  # Nov 2025
#         # rho=(0.0, 24 * u.mm),
#         # absZ=(0.0, 1.0 * u.m),
#     ),
# )

# s = addHoughVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     inputSpacePoints=addSpacePointsMaking(
#         s,
#         trackingGeometry,
#         geo_dir / "geoSelectionForSeedingInner_BARREL_LayersOnly.json",
#         outputName = "spacePointsForHoughVertexing",
#         logLevel = acts.logging.VERBOSE
#     ),
#     logLevel = acts.logging.VERBOSE
# )


strWhichSeedingLayers = ""
if pars.seedingLayers == "VD":
    strWhichSeedingLayers = "geoSelectionForSeeding_VD.json"
elif pars.seedingLayers == "ML3":
    strWhichSeedingLayers = "geoSelectionForSeeding_3firstML.json"
elif pars.seedingLayers == "MLall":
    strWhichSeedingLayers = "geoSelectionForSeeding_5ML.json"
elif pars.seedingLayers == "VDML":
    strWhichSeedingLayers = "geoSelectionForSeeding_IB_ML.json"
print("strWhichLayers=", strWhichSeedingLayers)


s = addSeeding(
    s,
    trackingGeometry,
    field,
    SeedFinderConfigArg(
        r=(None, 210 * u.mm),  # iTOF is at 190 mm! if we want it for seeding
        # r=(None, 150 * u.mm),
        # r=(None, 30 * u.mm),
        # deltaR=(1 * u.mm, 120 * u.mm),  # deltaR=(1. * u.mm, 60 * u.mm),
        deltaR=(1 * u.mm, 200 * u.mm),  # deltaR=(1. * u.mm, 60 * u.mm),
        collisionRegion=(
            -IA_collisionRegion_forSeeds * u.mm,
            IA_collisionRegion_forSeeds * u.mm,
        ),
        z=(-1000 * u.mm, 1000 * u.mm),
        maxSeedsPerSpM=pars.maxSeedsPerSpM,  # 2 is minimum, >2 is better for Pb-Pb
        sigmaScattering=pars.sigmaScattering,
        # more info: https://github.com/acts-project/acts/blob/main/Core/include/Acts/Seeding/SeedFinderConfig.hpp
        radLengthPerSeed=pars.radLengthPerSeed,
        minPt=pars.minSeedPt * u.GeV,
        impactMax=pars.impParForSeeds
        * u.mm,  # important! IB vs ML seeds (e.g. 1 mm is ok for IB seeds, 5 mm - for ML seeds)
        cotThetaMax=27.2899,
        seedConfirmation=True,
        centralSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-620 * u.mm,
            zMaxSeedConf=620 * u.mm,
            rMaxSeedConf=4.9
            * u.mm,  # 36 * u.mm,  # IA: dramatically affects acceptance at eta ~4. <5 * u.mm  gives best results
            nTopForLargeR=1,  # number of top space points that confirm my seed at larger R, 1 - no confirmation
            nTopForSmallR=2,
        ),
        forwardSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-1220 * u.mm,
            zMaxSeedConf=1220 * u.mm,
            rMaxSeedConf=26 * u.mm,  # 15 * u.mm,  #36 * u.mm,
            nTopForLargeR=1,
            nTopForSmallR=2,
        ),
        # skipPreviousTopSP=True,
        useVariableMiddleSPRange=True,
        deltaRMiddleSPRange=(0.2 * u.mm, 1.0 * u.mm),
    ),
    SeedFinderOptionsArg(bFieldInZ=cfg.seeding.bField * u.T, beamPos=(0 * u.mm, 0 * u.mm)),
    SeedFilterConfigArg(
        seedConfirmation=True if pars.impParForSeeds < 2.0 else False,  # mm
        # If seedConfirmation is true we classify seeds as "high-quality" seeds.
        # Seeds that are not confirmed as "high-quality" are only selected if no
        # other "high-quality" seed has been found for that inner-middle doublet
        # Maximum number of normal seeds (not classified as "high-quality" seeds)
        # in seed confirmation
        maxSeedsPerSpMConf=IA_maxSeedsPerSpMConf,  # 1 - USED_FOR_AUG_2025,#3, # CRUCIAL!!!!!!
        # 1 - USED_FOR_AUG_2025,   # Core/include/Acts/Seeding/SeedFilterConfig.hpp
        maxQualitySeedsPerSpMConf=IA_maxQualitySeedsPerSpMConf,
        # Maximum number of "high-quality" seeds for each inner-middle SP-dublet in
        # seed confirmation. If the limit is reached we check if there is a normal
        # quality seed to be replaced
    ),
    SpacePointGridConfigArg(
        # zBinEdges=[
        # -4000.0,
        # -2500.0,
        # -2000.0,
        # -1320.0,
        # -625.0,
        # -350.0,
        # -250.0,
        # 250.0,
        # 350.0,
        # 625.0,
        # 1320.0,
        # 2000.0,
        # 2500.0,
        # 4000.0,
        # ],
        impactMax=1.0 * u.mm,
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
    # particleHypothesis=acts.ParticleHypothesis.pion, # IA
    # particleHypothesis=acts.ParticleHypothesis.electron, # IA
    geoSelectionConfigFile=geo_dir / "../seedingConfigurations" / strWhichSeedingLayers,
    # seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
    seedingAlgorithm=IA_whichSeedingAlg,  # SeedingAlgorithm.Default, # Nov 2025: default is GridTriplet !!!
    outputDirRoot=outputDir,
    #    initialVarInflation = (50,50,50,50,50,50)  # IA
    #    initialVarInflation = (0.2,0.2,0.2,0.2,0.2,0.2)  #IA
)

# from https://github.com/acts-project/acts/blob/v41.0.0/Examples/Scripts/Python/hashing_seeding.py:
if False:  # if we want to save seeds to root file
    rootSeedsWriter = acts.examples.RootSeedWriter(
        level=acts.logging.VERBOSE,
        inputSeeds="seeds",
        filePath=str(outputDir / "seeds.root"),
        writingMode="big",  # IA
    )
    s.addWriter(rootSeedsWriter)


# s = addHoughVertexFinding(
#     s,
#     outputDirRoot=outputDir,
#     inputSpacePoints="spacepoints",
#     logLevel = acts.logging.VERBOSE
# )

### useful info: https://github.com/acts-project/acts/blob/main/Core/include/Acts/TrackFinding/MeasurementSelector.hpp
s = addCKFTracks(
    s,
    trackingGeometry,
    field,
    TrackSelectorConfig(
        pt=(0.06 * u.GeV, 120 * u.GeV),
        nMeasurementsMin=pars.nMeasurementsMin,
        maxSharedHits=pars.maxSharedHits,
    ),  # IA, was: 500.0 * u.MeV
    ckfConfig=CkfConfig(
        chi2CutOffOutlier=pars.ckfChi2Outlier,
        chi2CutOffMeasurement=pars.ckfChi2Measurement,
        numMeasurementsCutOff=pars.ckfMeasPerSurf,
        seedDeduplication=pars.seedDeduplication,
        stayOnSeed=pars.stayOnSeed,
    ),
    twoWay=pars.twoWayCKF,  # default: True,
    outputDirRoot=outputDir,
    writeTrackSummary=False,
    logLevel=acts.logging.INFO,
)

s = addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(
        maximumSharedHits=pars.maxSharedHits, nMeasurementsMin=pars.nMeasurementsMin
    ),
    outputDirRoot=outputDir,
    logLevel=acts.logging.INFO,
)

s = addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.AMVF,
    outputDirRoot=outputDir,
    seeder=acts.examples.VertexSeedFinder.AdaptiveGridSeeder,
    useTime=False,  # True,
)

s.run()

### save script used

shutil.copyfile(os.path.abspath(__file__), os.path.join(os.getcwd(), IA_outputDirName, os.path.basename(__file__)))
# with open("lastRun.txt", "w") as fh:
#     fh.write(IA_outputDirName + "\n")

### save command line used

cmd = " ".join(sys.argv)
with open(os.path.join(os.getcwd(), IA_outputDirName, "full_chain_used_command_line_args.log"), "a") as f:
    f.write(f"{datetime.now()}\n")
    f.write(f"{cmd}\n")
