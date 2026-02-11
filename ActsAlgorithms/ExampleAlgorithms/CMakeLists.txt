add_library(ExampleAlgorithms SHARED)

# ----------------------------------------
# Sources
# ----------------------------------------
target_sources(ExampleAlgorithms
    PRIVATE
    src/ExampleAlgorithm.cpp
    src/TrackTruthMatcher.cpp
    src/HitRemoverAlgorithm.cpp
    src/TrackMergerAlgorithm.cpp
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

#link_directories(${CMAKE_PREFIX_PATH}/lib64)

# ----------------------------------------
# Link against Acts
# ----------------------------------------
target_link_libraries(ExampleAlgorithms
    PUBLIC
    Acts::Core
    ${CMAKE_PREFIX_PATH}/lib64/libActsExamplesFramework.so
  )
