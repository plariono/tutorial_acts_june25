#include "TrackMergerAlgorithm.hpp"
#include "Acts/Utilities/Logger.hpp"

using namespace ActsExamples;

namespace AliceActsTrk {

  TrackMergerAlgorithm::TrackMergerAlgorithm(TrackMergerAlgorithm::Config cfg,
                                              Acts::Logging::Level lvl)
    : IAlgorithm("TrackMergerAlgorithm", lvl), m_cfg(std::move(cfg)) {
  
    
    for (size_t itc =0; itc < m_cfg.inputTrackCollections.size(); itc++) {
      
      auto tc = m_cfg.inputTrackCollections[itc];
      m_inputTrackCollections.push_back(ReadDataHandle<ConstTrackContainer>{this,tc});
      m_inputTrackCollections[itc].initialize(tc);
    }
    
    m_outputTracks.initialize(m_cfg.outputTrackCollection);
    
  }
  
  ProcessCode TrackMergerAlgorithm::execute(const AlgorithmContext& ctx) const {


    // Read input data

    std::vector<ConstTrackContainer> trackCollections;

    for (size_t itc =0; itc < m_cfg.inputTrackCollections.size(); itc++) {
      trackCollections.push_back(m_inputTrackCollections.at(itc)(ctx));
    }

    return ProcessCode::SUCCESS;
  }

  
  
}
