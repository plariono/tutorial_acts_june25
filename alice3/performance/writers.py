from pathlib import Path
from typing import Optional, Union, List
from enum import Enum
from collections import namedtuple

import acts
import acts.examples
import acts.examples.reconstruction as acts_reco

#Alice3 specific algorithms
from AliceActsPythonBindings import TrackTruthMatcher
from AliceActsPythonBindings import HitRemoverAlgorithm
from AliceActsPythonBindings import TrackMergerAlgorithm

# Alice3 plotting
import alice3.performance.plotting as alice3_plotting

def addCKFTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    trackSelectorConfig: Optional[
        Union[acts_reco.TrackSelectorConfig, List[acts_reco.TrackSelectorConfig]]
    ] = None,
    ckfConfig: acts_reco.CkfConfig = acts_reco.CkfConfig(),
    twoWay: bool = True,
    reverseSearch: bool = False,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat=False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for ROOT output, None triggers no output
    trackSelectorConfig : TrackSelectorConfig(loc0, loc1, time, eta, absEta, pt, phi, minMeasurements)
        TrackSelector configuration. Each range is specified as a tuple of (min,max).
        Specify as a list(TrackSelectorConfig) for eta-dependent cuts, with binning specified by absEta[1].
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackSelector.hpp
    writeTrackSummary : bool, True
        write tracksummary_ckf.root ntuple?
    writeTrackStates : bool, False
        write trackstates_ckf.root ntuple? This can be quite large.
    writePerformance : bool, True
        write performance_fitting_ckf.root and performance_finding_ckf.root ntuples?
    writeCovMat : bool, False
        write covaraiance matrices to tracksummary_ckf.root ntuple?
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    tslist = (
        []
        if trackSelectorConfig is None
        else (
            [trackSelectorConfig]
            if type(trackSelectorConfig) is acts_reco.TrackSelectorConfig
            else trackSelectorConfig
        )
    )

    if len(tslist) > 1:
        cutSets = []
        for c in tslist:
            defKW = acts_reco.trackSelectorDefaultKWArgs(c)
            defKW.pop("absEtaMax", None)
            cutSets += [acts.TrackSelector.Config(**(defKW))]
    else:
        cutSets = [
            acts.TrackSelector.Config(**(acts_reco.trackSelectorDefaultKWArgs(c))) for c in tslist
        ]

    if len(tslist) == 0:
        trkSelCfg = None
    elif len(tslist) == 1:
        trkSelCfg = cutSets[0]
    else:
        trkSelCfg = acts.TrackSelector.EtaBinnedConfig(
            cutSets=cutSets,
            absEtaEdges=[cutSets[0].absEtaMin] + [c.absEta[1] for c in tslist],
        )

    # Setup the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=customLogLevel(),
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [
                (
                    acts.GeometryIdentifier(),
                    (
                        [],
                        [ckfConfig.chi2CutOffMeasurement],
                        [ckfConfig.chi2CutOffOutlier],
                        [ckfConfig.numMeasurementsCutOff],
                    ),
                )
            ]
        ),
        inputMeasurements="measurements",
        inputInitialTrackParameters="estimatedparameters",
        inputSeeds=(
            "estimatedseeds"
            if ckfConfig.seedDeduplication or ckfConfig.stayOnSeed
            else ""
        ),
        outputTracks="ckf_tracks",
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field, customLogLevel()
        ),
        **acts.examples.defaultKWArgs(
            trackingGeometry=trackingGeometry,
            magneticField=field,
            trackSelectorCfg=trkSelCfg,
            maxSteps=ckfConfig.maxSteps,
            twoWay=twoWay,
            reverseSearch=reverseSearch,
            seedDeduplication=ckfConfig.seedDeduplication,
            stayOnSeed=ckfConfig.stayOnSeed,
            pixelVolumeIds=ckfConfig.pixelVolumes,
            stripVolumeIds=ckfConfig.stripVolumes,
            maxPixelHoles=ckfConfig.maxPixelHoles,
            maxStripHoles=ckfConfig.maxStripHoles,
            trimTracks=ckfConfig.trimTracks,
            constrainToVolumeIds=ckfConfig.constrainToVolumes,
            endOfWorldVolumeIds=ckfConfig.endOfWorldVolumes,
        ),
    )
    s.addAlgorithm(trackFinder)
    s.addWhiteboardAlias("tracks", trackFinder.config.outputTracks)


    truthMatchCfg = TrackTruthMatcher.Config()
    truthMatchCfg.inputTracks=trackFinder.config.outputTracks
    truthMatchCfg.inputParticles="particles_selected"
    truthMatchCfg.inputMeasurementParticlesMap="measurement_particles_map"
    truthMatchCfg.outputTrackParticleMatching="ckf_track_particle_matching"
    truthMatchCfg.outputParticleTrackMatching="ckf_particle_track_matching"
    truthMatchCfg.looperProtection=True
    truthMatchCfg.loop_absEta = 1.5
    truthMatchCfg.loop_maxPt  = 0.2
    truthMatchCfg.loop_maxParticleHits = 11
    
    
    matchAlg = TrackTruthMatcher(
        config = truthMatchCfg,
        level=customLogLevel())
    
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name="ckf",
        tracks=trackFinder.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=writePerformance,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s

def addTrackWriters(
    s: acts.examples.Sequencer,
    name: str,
    tracks: str = "tracks",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeSummary: bool = True,
    writeStates: bool = False,
    writeFitterPerformance: bool = False,
    writeFinderPerformance: bool = False,
    logLevel: Optional[acts.logging.Level] = None,
    writeCovMat=False,
):
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        if writeSummary:
            trackSummaryWriter = acts.examples.root.RootTrackSummaryWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                filePath=str(outputDirRoot / f"tracksummary_{name}.root"),
                treeName="tracksummary",
                writeCovMat=writeCovMat,
            )
            s.addWriter(trackSummaryWriter)

        if writeStates:
            trackStatesWriter = acts.examples.root.RootTrackStatesWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                inputSimHits="simhits",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                filePath=str(outputDirRoot / f"trackstates_{name}.root"),
                treeName="trackstates",
            )
            s.addWriter(trackStatesWriter)

        if writeFitterPerformance:

            trackFitterPerformanceWriter = (
                acts.examples.root.RootTrackFitterPerformanceWriter(
                    level=customLogLevel(),
                    inputTracks=tracks,
                    inputParticles="particles_selected",
                    inputTrackParticleMatching="track_particle_matching",
                    resPlotToolConfig = alice3_plotting.resPlotToolConfig,
                    filePath=str(outputDirRoot / f"performance_fitting_{name}.root"),
                )
            )
            s.addWriter(trackFitterPerformanceWriter)

        if writeFinderPerformance:

            trackFinderPerfWriter = acts.examples.root.RootTrackFinderPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                inputParticleTrackMatching="particle_track_matching",
                inputParticleMeasurementsMap="particle_measurements_map",
                effPlotToolConfig = alice3_plotting.effPlotToolConfig,
                fakePlotToolConfig = alice3_plotting.fakePlotToolConfig,
                duplicationPlotToolConfig = alice3_plotting.duplicationPlotToolConfig,
                trackQualityPlotToolConfig = alice3_plotting.trackQualityPlotToolConfig,
                trackSummaryPlotToolConfig = alice3_plotting.trackSummaryPlotToolConfig,
                filePath=str(outputDirRoot / f"performance_finding_{name}.root"),
            )
            s.addWriter(trackFinderPerfWriter)



def addHitRemoverAlgorithm(
        s : acts.examples.Sequencer,
        inputMeasurements : str,
        inputTracks : str,
        inputMeasurementParticlesMap : str,
        sortByOldIndex : bool,
        used_meas_idxs : str,
        outputMeasurements : str,
        outputMeasurementParticlesMap : str,
        outputParticleMeasurementsMap : str,
        outputIndexingMap : str,
        logLevel : acts.logging.Level = None
        ):

    customLogLevel = acts.examples.defaultLogging(s, logLevel)


    cfg = HitRemoverAlgorithm.Config()
    cfg.inputMeasurements = inputMeasurements
    cfg.inputTracks       = inputTracks # I think we should use the resolved ones?
    cfg.inputMeasurementParticlesMap = inputMeasurementParticlesMap
    cfg.sortByOldIndex = sortByOldIndex
    cfg.usedIndices       = used_meas_idxs
    cfg.outputMeasurements= outputMeasurements
    cfg.outputMeasurementParticlesMap = outputMeasurementParticlesMap
    cfg.outputParticleMeasurementsMap = outputParticleMeasurementsMap
    cfg.outputIndexingMap = outputIndexingMap
    
    HitRemoverAlg = HitRemoverAlgorithm(
        config = cfg,
        level=customLogLevel())
    
    s.addAlgorithm(HitRemoverAlg)
    
    
def addTrackPerformanceWriters(
    sequence: acts.examples.Sequencer,
    outputDirRoot: Union[Path, str],
    tracks: str,
    prototracks: str,
    selectedParticles: str,
    inputParticles: str,
    outputTrackParameters: str,
    logLevel: acts.logging.Level = None,
):


    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    outputDirRoot = Path(outputDirRoot)
    if not outputDirRoot.exists():
        outputDirRoot.mkdir()

    sequence.addWriter(
        acts.examples.root.RootTrackFinderPerformanceWriter(
            level=customLogLevel(),
            inputTracks=tracks,
            inputParticles=selectedParticles,
            inputTrackParticleMatching="seed_particle_matching",
            inputParticleTrackMatching="particle_seed_matching",
            inputParticleMeasurementsMap="particle_measurements_map",
            filePath=str(outputDirRoot / f"performance_seeding.root"),
        )
    )

def addTrackMerger(
        sequence: acts.examples.Sequencer,
        inputTrackCollections : str,
        inputIndexingMaps : str,
        outputTrackCollection : str,
        logLevel : acts.logging.Level = None):

    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    
    cfg = TrackMergerAlgorithm.Config()
    cfg.inputTrackCollections  = inputTrackCollections
    cfg.inputIndexingMaps      = inputIndexingMaps
    cfg.outputTrackCollection  = outputTrackCollection

    TrackMergerAlg = TrackMergerAlgorithm(
    config = cfg,
    level = customLogLevel())
    
        
    sequence.addAlgorithm(TrackMergerAlg)

    
