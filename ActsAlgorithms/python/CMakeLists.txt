# Add the examples python binding
pybind11_add_module(AliceActsPythonBindings
  src/AliceActsAlgorithmsModule.cpp
)

target_compile_features(AliceActsPythonBindings PUBLIC cxx_std_20)


target_include_directories(AliceActsPythonBindings
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../ExampleAlgorithms/include
    )

target_link_libraries(AliceActsPythonBindings
  PUBLIC
  ExampleAlgorithms
  /home/pibutti/sw/acts/install/lib64/libActsExamplesFramework.so
  Acts::Core
)
