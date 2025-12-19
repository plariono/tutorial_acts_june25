#pragma once


#include "ExampleAlgorithm.hpp"
#include <algorithm>

#include "ActsPython/Utilities/Macros.hpp"
#include <pybind11/pybind11.h>


namespace AliceActsPython {

// Inline ensures the symbol is available in the .so
inline void addExampleAlgorithm(pybind11::module& mex) {
    ACTS_PYTHON_DECLARE_ALGORITHM(
        AliceActsTrk::ExampleAlgorithm, mex, "ExampleAlgorithm", testString
    );
}

}  // namespace AliceActsPython
