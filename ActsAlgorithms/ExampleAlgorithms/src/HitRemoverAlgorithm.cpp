#include "HitRemoverAlgorithm.hpp"

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
    m_outputMeasurements.initialize(m_cfg.outputMeasurements);
    
  }


  ProcessCode HitRemoverAlgorithm::execute(const AlgorithmContext& ctx) const {

    // Read input data
    const auto& measurements = m_inputMeasurements(ctx);
    const auto& tracks = m_inputTracks(ctx);
    
    

    return ProcessCode::SUCCESS;
  }



  void HitRemoverAlgorithm::computeUsedHits(
                                            const TrackContainer& tracks,
                                            const MeasurementContainer& measurements) const {

    // Compute used hits from all the reconstruced tracks


    for (auto track : tracks) {
      for (auto state : track.trackStatesReversed()) {
        if (!state.typeFlags().isMeasurement()) {
          continue;
        }

        std::size_t hitIndex = state.getUncalibratedSourceLink()
          .template get<IndexSourceLink>()
          .index();
        
      }
    } // loop on tracks
    
  }
  

}
