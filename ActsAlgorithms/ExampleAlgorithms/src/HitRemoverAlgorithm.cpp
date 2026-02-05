#include "HitRemoverAlgorithm.hpp"
#include "Acts/Utilities/Logger.hpp"

using namespace ActsExamples;

namespace AliceActsTrk {

  HitRemoverAlgorithm::HitRemoverAlgorithm(HitRemoverAlgorithm::Config cfg,
                                           Acts::Logging::Level lvl)
    : IAlgorithm("HitRemoverAlgorithm", lvl), m_cfg(std::move(cfg)) {

    if (m_cfg.inputMeasurements.empty()) {
      throw std::invalid_argument("Missing input hit collection");
    }

    if (m_cfg.inputTracks.empty()) {
      throw std::invalid_argument("Missing input tracks");
    }

    if (m_cfg.outputMeasurements.empty()) {
      throw std::invalid_argument("Missing output hit collection");
    }
    
    m_inputMeasurements.initialize(m_cfg.inputMeasurements);
    m_inputTracks.initialize(m_cfg.inputTracks);
    m_usedIndices.initialize(m_cfg.usedIndices);
    m_outputMeasurements.initialize(m_cfg.outputMeasurements);
    
  }


  ProcessCode HitRemoverAlgorithm::execute(const AlgorithmContext& ctx) const {

    // Read input data
    const auto& measurements = m_inputMeasurements(ctx);
    const auto& tracks = m_inputTracks(ctx);
    
    
    ACTS_VERBOSE("Measurements size:: "<< measurements.size());
    ACTS_VERBOSE("Tracks size:: "<< tracks.size());

    // O(1) lookup and shared hits accounted for
    std::unordered_set<size_t> usedIndices;
    MeasurementContainer filteredMeasurements;
    
    computeUsedHits(tracks,
                    measurements,
                    usedIndices,
                    filteredMeasurements);

    ACTS_VERBOSE("Used Measurements size:: "<< usedIndices.size());
    ACTS_VERBOSE("Filtered Measurements size:: "<< filteredMeasurements.size());

    m_usedIndices(ctx, std::move(usedIndices));
    m_outputMeasurements(ctx, std::move(filteredMeasurements));
    
    return ProcessCode::SUCCESS;
  }



  void HitRemoverAlgorithm::computeUsedHits(const ConstTrackContainer& tracks,
                                            const MeasurementContainer& measurements,
                                            std::unordered_set<size_t>& usedIndices,
                                            MeasurementContainer& filteredMeasurements) const {

    
    // Compute used hits from all the reconstruced tracks
    
    for (auto track : tracks) {
      for (auto state : track.trackStatesReversed()) {
        if (!state.typeFlags().isMeasurement()) {
          continue;
        }

        std::size_t hitIndex = state.getUncalibratedSourceLink()
          .template get<IndexSourceLink>()
          .index();

        usedIndices.insert(hitIndex);
        
        
      }
    } // loop on tracks

    filteredMeasurements.reserve(measurements.size() - usedIndices.size());
    for (size_t i = 0; i <  measurements.size(); ++i) {

      if (usedIndices.find(i) == usedIndices.end()) {

        auto meas = measurements.getMeasurement(i);
        filteredMeasurements.copyMeasurement(meas);
      }
    }
  }
} //AliceActsTrk
