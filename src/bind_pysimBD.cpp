//
//  bind_pysimBD.cpp
//  BDmodel
//
//  Created by Francesco PINOTTI on 20/10/2025.
//

#include "pysimBD.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MODULE(pysimBD, m) {
    m.doc() = "python binding for c++ code simulating SEIR dynamics in a market"; // optional module docstring
    
    //==== Individual-based model functions
    
    m.def("simulate_BD_tree", &simulate_BD, "Returns diversity metrics single type simulations",
          py::arg("seed"),
          py::arg("max_cases"),
          py::arg("max_samples"),
          py::arg("R0"),
          py::arg("dI"),
          py::arg("rho") ) ;
    
}
