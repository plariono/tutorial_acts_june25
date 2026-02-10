#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "Acts/EventData/TrackContainer.hpp" 
#include "ActsExamples/Framework/DataHandle.hpp"

#include <memory> 
#include <vector>
#include <string>


using namespace ActsExamples;


namespace AliceActsTrk {

  class TrackMergerAlgorithm final : public IAlgorithm {

  public:

    struct Config {
      std::vector<std::string>  inputTrackCollections{};
      std::string outputTrackCollection = "";
    };

    TrackMergerAlgorithm(Config cfg, Acts::Logging::Level lvl);
    ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;

    const Config& config() const {return m_cfg;}

  private:
    Config m_cfg;
    std::vector<ReadDataHandle<ConstTrackContainer>> m_inputTrackCollections;
    WriteDataHandle<TrackContainer> m_outputTracks{this, "outputTracks"};
    
  }; //TrackMergerAlgorithm
}
