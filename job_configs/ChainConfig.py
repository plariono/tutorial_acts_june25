from typing import Literal,Tuple
from dataclasses import dataclass


@dataclass
class General:
    fieldMap: str = "fieldmaps/solenoid_R1625_L7500_B2T_scaled.txt"
    MF: float = 2.0 #T

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
    
class Config:
    general = General()
    pythia = PythiaConfig()
    detSim = detSimConfig()
    particleGun = ParticleGunConfig()
    seeding = seedingConfig()
