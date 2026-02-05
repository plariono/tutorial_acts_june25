from pathlib import Path
from typing import Optional, Union, List
from enum import Enum
from collections import namedtuple

import acts
import acts.examples
import acts.examples.reconstruction as acts_reco

#######################
#### SOME OTHER PARAMS#
#######################

impParForSeeds = 2.0
collisionRegion_forSeeds = 250 # 1000 # mm; large values - for V0 daughter reconstruction
seedingAlg = SeedingAlgorithm.GridTriplet
maxSeedsPerSPMConf = 1
maxQualitySeedsPerSPMConf = 1
enableMaterial = true


##Parameters to be passed
maxSeedsPerSpM   = 2     ## ok for pp, for PbPb higher
sigmaScattering  = 5.0   ## Tuned with particle gun, combine radLengthPerSeed = 0.05
radLengthPerSeed = 0.05
minSeedPt        = 0.08  ## to be changed per iteration
impParForSeeds   = 1.0   #mm
bFieldInZ        = 2.0   #T


DefaultSeedFinderConfig = SeedFinderConfigArg(
        r=(None, 210 * u.mm),  # iTOF is at 190 mm! if we want it for seeding
        # r=(None, 150 * u.mm),
        # r=(None, 30 * u.mm),
        deltaR=(1 * u.mm, 200 * u.mm),  # deltaR=(1. * u.mm, 60 * u.mm),
        collisionRegion=(
            -collisionRegion_forSeeds * u.mm,
            collisionRegion_forSeeds * u.mm,
        ),
        z=(-1000 * u.mm, 1000 * u.mm),
        maxSeedsPerSpM=maxSeedsPerSpM,  # 2 is minimum, >2 is better for Pb-Pb
        sigmaScattering=sigmaScattering,
        radLengthPerSeed=radLengthPerSeed,  
        minPt=minSeedPt * u.GeV,
        impactMax=impParForSeeds * u.mm,  # important! IB vs ML seeds (e.g. 1 mm is ok for IB seeds, 5 mm - for ML seeds)
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
    )

DefaultSeedFinderOptionsArg = SeedFinderOptionsArg(bFieldInZ=bFieldInZ * u.T, beamPos=(0 * u.mm, 0 * u.mm))

DefaultSeedFilterConfigArg = SeedFilterConfigArg(
        seedConfirmation=True if impParForSeeds < 2.0 else False,  # mm
        # If seedConfirmation is true we classify seeds as "high-quality" seeds.
        # Seeds that are not confirmed as "high-quality" are only selected if no
        # other "high-quality" seed has been found for that inner-middle doublet
        # Maximum number of normal seeds (not classified as "high-quality" seeds)
        # in seed confirmation
        maxSeedsPerSpMConf=maxSeedsPerSpMConf,  # 1 - USED_FOR_AUG_2025,#3, # CRUCIAL!!!!!!
        maxQualitySeedsPerSpMConf=maxQualitySeedsPerSpMConf,  # 1 - USED_FOR_AUG_2025,   # Core/include/Acts/Seeding/SeedFilterConfig.hpp
        # Maximum number of "high-quality" seeds for each inner-middle SP-dublet in
        # seed confirmation. If the limit is reached we check if there is a normal
        # quality seed to be replaced
    )


DefaultSpacePointGridConfigArg = SpacePointGridConfigArg(
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
    )
DefaultSeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(
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
    )

s = addSeeding(
    s,
    trackingGeometry,
    field,
    DefaultSeedFinderConfigArg,
    DefaultSeedFinderOptionsArg
    DefaultSeedFilterConfigArg,
    DefaultSpacePointGridConfigArg,
    DefaultSeedingAlgorithmConfigArg,
    # particleHypothesis=acts.ParticleHypothesis.pion, # IA
    # particleHypothesis=acts.ParticleHypothesis.electron, # IA
    geoSelectionConfigFile=geo_dir / strWhichSeedingLayers,
    # seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
    seedingAlgorithm=seedingAlg,  
    outputDirRoot=outputDir,
    #    initialVarInflation = (50,50,50,50,50,50)  # IA
    #    initialVarInflation = (0.2,0.2,0.2,0.2,0.2,0.2)  #IA
)
