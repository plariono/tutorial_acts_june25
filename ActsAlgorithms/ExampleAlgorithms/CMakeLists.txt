add_library(ExampleAlgorithms SHARED)

# ----------------------------------------
# Sources
# ----------------------------------------
target_sources(ExampleAlgorithms
    PRIVATE
    src/ExampleAlgorithm.cpp
    src/TrackTruthMatcher.cpp
        # add more .cpp files here
)

# ----------------------------------------
# Include directories
# ----------------------------------------
target_include_directories(ExampleAlgorithms
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#        $<INSTALL_INTERFACE:include>
)

# ----------------------------------------
# Link against Acts
# ----------------------------------------
target_link_libraries(ExampleAlgorithms
    PUBLIC
    Acts::Core
    /home/pibutti/sw/acts/install/lib64/libActsExamplesFramework.so
  )

