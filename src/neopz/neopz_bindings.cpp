#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

namespace py = pybind11;
using namespace py::literals;

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzgmesh.h"

PYBIND11_MODULE(neopz, m) {
    m.doc() = R"pbdoc(
        -------------------------
        Python bindings for NeoPZ
        -------------------------
    )pbdoc";

    py::class_<TPZManVector<double>>(m, "PZVecDouble")
       .def(py::init())
       .def(py::init<int64_t>())
       .def(py::init<int64_t, double>())
       .def("Resize", [](TPZManVector<double>& vec, const int64_t& newsize) { return vec.Resize(newsize); })
       .def("Size", [](const TPZManVector<double>& vec) { return vec.size(); })
       .def("__getitem__",
           [](const TPZManVector<double>& vec, int64_t position) {
               if (position >= vec.size() || position < 0) throw py::index_error();
               return vec[position];
           },
           py::is_operator()
           )
       .def("__setitem__",
           [](TPZManVector<double>& vec, int64_t position, double value) {
               if (position >= vec.size() || position < 0) throw py::index_error();
               vec[position] = value;
           },
           py::is_operator()
           )
       ;
   py::class_<TPZFMatrix<double>>(m, "PZMatrix")
       .def(py::init())
       .def(py::init<int64_t, int64_t>())
       .def(py::init<int64_t, int64_t, double>())
       ;

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
