#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

#include <memory>
#include <string>
#include <vector>

using namespace ActsExamples;

namespace AliceActsTrk {

  class HitRemoverAlgorithm final : public IAlgorithm {
    
  public:

    struct Config {
      std::string inputMeasurements = "";
      std::string inputTracks = "";
      std::string outputMeasurements = "";
      
    };

    HitRemoverAlgorithm(Config cfg, Acts::Logging::Level lvl);
    ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;

    const Config& config() const { return m_cfg; }

  private:
    Config m_cfg;

    void computeUsedHits(const TrackContainer& tracks,
                         const MeasurementContainer& measurements) const;


    ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
      "InputMeasurements"};
    ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
        
    WriteDataHandle<MeasurementContainer> m_outputMeasurements{this, "OutputMeasurments"};
    
  }; //HitRemoverAlgorithm
    
} //namespace AliceActsTrk


