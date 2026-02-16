from typing import Literal,Tuple
from dataclasses import dataclass


@dataclass
class General:
    fieldMap: str = "fieldmaps/solenoid_R1625_L7500_B2T_scaled.txt"
    MF: float = 2.0 #T

    #Geo 1
    geo_dir: str = "geometries/geometry_2026_01_30_mockupTiledDisks_for_Geant_tests/geom"
    digi_file: str = "geometries/geometry_2026_01_30_mockupTiledDisks_for_Geant_tests/digiConfigurations/digi-smearing-config_no_TOFs_iTOF_removed.json"

    # Geo 2
    #geo_dir: str = "geometries/geometry_2026_02_12_default_v3_cylindrical_barrel_tiled_disks/geom"
    #digi_file: str = "geometries/geometry_2026_02_12_default_v3_cylindrical_barrel_tiled_disks/digiConfigurations/digi-smearing-config_with_TOFs_noEndcapTOFs_noTimeInTOFs.json"

    enableMaterial : bool = True
    
@dataclass
class ParticleGunConfig:
    gunMult: int = 1
    gunPtRange: Tuple[float,float]  = (0.1,5.0)
    gunEtaRange: Tuple[float,float] = (-2.5,2.5)
    gunPhiRange: Tuple[float,float] = (-3.1415926, 3.1415926)
    gunPID: int = 211
    gunRandCharge: bool = True
    meanVz: float = 0.
    sigmaVz: float = 0.180 ## ns

@dataclass
class PythiaConfig:
    usePythia: bool = False
    system: Literal["pp", "PbPbMB", "PbPbCentral"] = "pp"
    pileup: int = 1
    meanVz: float = 0.
    sigmaVz: float = 50.

@dataclass
class detSimConfig:
    simulation: Literal["Fatras","Geant4"] = "Geant4"
    minFatrasPt: float = 0.001
    particleSelectionEta: float = 4.2

@dataclass
class seedingConfig:
    bField : float = 2.0 #T
    seedingLayers: str = "geoSelectionForSeeding_VD.json"
    seedingAlgo: Literal["GridTriplet","TruthSmeared"] = "GridTriplet"
    minSeedPt: float = 0.08
    impParForSeeds: float = 1.0
    sigmaScattering: float = 5
    radLengthPerSeed: float = 0.05 
    maxSeedsPerSpM: int = 2 
    filter_maxSeedsPerSpMConf: int = 1
    filter_maxQualitySeedsPerSpMConf: int = 1
    
    



@dataclass
class trackingConfig:
    nMeasurementsMin: int = 7
    ckfMeasPerSurf: int = 1
    ckfChi2Measurement: float = 45
    ckfChi2Outlier: float = 100
    seedDeduplication: bool = True
    stayOnSeed: bool = False
    maxSharedHits: int = 2
    twoWayCKF: bool = True
    
    
class Config:
    general = General()
    pythia = PythiaConfig()
    detSim = detSimConfig()
    particleGun = ParticleGunConfig()
    seeding = seedingConfig()
    tracking = trackingConfig()
