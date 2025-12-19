#include "ExampleAlgorithm.hpp"

using namespace ActsExamples;

namespace AliceActsTrk {

  ExampleAlgorithm::ExampleAlgorithm(ExampleAlgorithm::Config cfg,
                                     Acts::Logging::Level lvl)
    : IAlgorithm("ExampleAlgorithm", lvl), m_cfg(std::move(cfg)) {
    
    std::cout<<"PF:: Constructed ExampleAlgorithm"<<std::endl;
  }
  
  
  ProcessCode ExampleAlgorithm::execute(const AlgorithmContext& ctx) const {
    
    std::cout<<"PF:: Executing ExampleAlgorithm" <<std::endl;
    
    return ProcessCode::SUCCESS;
    
  }
  
}

