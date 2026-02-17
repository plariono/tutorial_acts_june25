# Add the examples python binding
pybind11_add_module(AliceActsPythonBindings
  src/AliceActsAlgorithmsModule.cpp
)

target_compile_features(AliceActsPythonBindings PUBLIC cxx_std_20)


target_include_directories(AliceActsPythonBindings
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../ExampleAlgorithms/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../Validation/include
    )

#link_directories(${CMAKE_PREFIX_PATH}/lib64)
#message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
    
target_link_libraries(AliceActsPythonBindings
  PUBLIC
  ExampleAlgorithms
  Validation
  ${CMAKE_PREFIX_PATH}/lib64/libActsExamplesFramework.so
  ${CMAKE_PREFIX_PATH}/lib64/libActsExamplesIoRoot.so
  Acts::Core
)
