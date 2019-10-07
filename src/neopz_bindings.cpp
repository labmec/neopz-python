//#include <pybind11/pybind11.h>
//#include <pybind11/operators.h>
//
//namespace py = pybind11;
//using namespace py::literals;
//
//#include "pzmanvector.h"
//#include "pzfmatrix.h"
//#include "tpzpoint.h"
//#include "tpzline.h"
//#include "tpztriangle.h"
//#include "tpzquadrilateral.h"
//
//PYBIND11_MODULE(neopz, m) {
//    m.doc() = R"pbdoc(
//        -------------------------
//        Python bindings for NeoPZ
//        -------------------------
//    )pbdoc";
//
//    // TPZManVector<double> bindings
//    py::class_<TPZManVector<double>>(m, "PZVecDouble")
//        .def(py::init())
//        .def(py::init<int64_t>())
//        .def(py::init<int64_t, double>())
//        .def("Resize", [](TPZManVector<double>& vec, const int64_t& newsize) { return vec.Resize(newsize); })
//        .def("Size", [](const TPZManVector<double>& vec) { return vec.size(); })
//        .def("__getitem__",
//            [](const TPZManVector<double>& vec, int64_t position) {
//                if (position >= vec.size() || position < 0) throw py::index_error();
//                return vec[position];
//            },
//            py::is_operator()
//        )
//        .def("__setitem__",
//            [](TPZManVector<double>& vec, int64_t position, double value) {
//                if (position >= vec.size() || position < 0) throw py::index_error();
//                vec[position] = value;
//            },
//            py::is_operator()
//        )
//    ;
//    // TPZFMatrix<double> bindings
//    py::class_<TPZFMatrix<double>>(m, "PZMatrix")
//        .def(py::init())
//        .def(py::init<int64_t, int64_t>())
//        .def(py::init<int64_t, int64_t, double>())
//        .def("GetVal", [](TPZFMatrix<double>& matrix, const int64_t& row, const int64_t& col) {
//            if (row >= matrix.Rows() || row < 0) throw py::index_error();
//            if (col >= matrix.Cols() || col < 0) throw py::index_error();
//            return matrix.GetVal(row, col);
//        })
//    ;
//
//    // TPZTriangle bindings
//    py::class_<pztopology::TPZTriangle>(m, "PZTriangle")
//            .def(py::init())
//            .def_static("SideDimension", &pztopology::TPZTriangle::SideDimension)
//    ;
//
//    // TPZQuadrilateral bindings
//    py::class_<pztopology::TPZQuadrilateral>(m, "PZQuad")
//        .def(py::init())
//        .def_static("SideDimension", &pztopology::TPZQuadrilateral::SideDimension)
//    ;
//
//#ifdef VERSION_INFO
//    m.attr("__version__") = VERSION_INFO;
//#else
//    m.attr("__version__") = "dev";
//#endif
//}
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

namespace py = pybind11;
using namespace py::literals;

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pztrnsform.h"

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzpyramid.h"
#include "tpzprism.h"
#include "tpzcube.h"

PYBIND11_MODULE(neopz, m) {
    m.doc() = R"pbdoc(
        -------------------------
        Python bindings for NeoPZ
        -------------------------
    )pbdoc";

    // TPZManVector<double> bindings
    py::class_<TPZManVector<double>>(m, "TPZVecDouble")
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

    // TPZFMatrix<double> bindings
    py::class_<TPZFMatrix<double>>(m, "TPZMatrix")
        .def(py::init())
        .def(py::init<int64_t, int64_t>())
        .def(py::init<int64_t, int64_t, double>())
        .def("GetVal", [](TPZFMatrix<double>& matrix, const int64_t& row, const int64_t& col) {
            if (row >= matrix.Rows() || row < 0) throw py::index_error();
            if (col >= matrix.Cols() || col < 0) throw py::index_error();
            return matrix.GetVal(row, col);
        })
    ;

    // TPZTransform bindings
    py::class_<TPZTransform<double>>(m, "TPZTransform")
        .def(py::init())
        .def(py::init<int>())
        .def(py::init<int, int>())
        .def("Mult", py::overload_cast<>(&TPZTransform<double>::Mult))
        .def("Sum", py::overload_cast<>(&TPZTransform<double>::Sum))
        .def("Mult", py::overload_cast<>(&TPZTransform<double>::Mult, py::const_))
        .def("Sum", py::overload_cast<>(&TPZTransform<double>::Sum, py::const_))
        .def("SetMatrix", &TPZTransform<double>::SetMatrix)
        .def("Multiply", &TPZTransform<double>::Multiply)
        .def("Apply", &TPZTransform<double>::Apply)
    ;

    // TPZPoint bindings
    py::class_<pztopology::TPZPoint>(m, "TPZPoint")
        .def(py::init())
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZPoint::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZPoint::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZPoint::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZPoint::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZPoint::SideNodeLocId)
        .def_static("CenterPoint", &pztopology::TPZPoint::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZPoint& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZPoint::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZPoint::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZPoint::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZPoint::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZPoint::CreateSideIntegrationRule)
    ;

    // TPZLine bindings
    py::class_<pztopology::TPZLine>(m, "TPZLine")
        .def(py::init())
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZLine::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZLine::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZLine::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZLine::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZLine::SideNodeLocId)
        .def_static("CenterPoint", &pztopology::TPZLine::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZLine& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZLine::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZLine::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZLine::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZLine::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZLine::CreateSideIntegrationRule)
    ;

    // TPZTriangle bindings
    py::class_<pztopology::TPZTriangle>(m, "TPZTriangle")
        .def(py::init())
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZTriangle::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZTriangle::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZTriangle::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZTriangle::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZTriangle::SideNodeLocId)
        .def_static("CenterPoint", &pztopology::TPZTriangle::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZTriangle& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZTriangle::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZTriangle::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZTriangle::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZTriangle::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZTriangle::CreateSideIntegrationRule)
    ;

    // TPZQuadrilateral bindings
    py::class_<pztopology::TPZQuadrilateral>(m, "TPZQuadrilateral")
        .def(py::init())
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZQuadrilateral::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZQuadrilateral::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZQuadrilateral::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZQuadrilateral::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZQuadrilateral::SideNodeLocId)
        .def_static("CenterPoint", &pztopology::TPZQuadrilateral::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZQuadrilateral& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZQuadrilateral::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZQuadrilateral::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZQuadrilateral::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZQuadrilateral::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZQuadrilateral::CreateSideIntegrationRule)
    ;

    // TPZTetrahedron bindings
    py::class_<pztopology::TPZTetrahedron>(m, "TPZTetrahedron")
        .def(py::init())
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZTetrahedron::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZTetrahedron::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZTetrahedron::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZTetrahedron::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZTetrahedron::SideNodeLocId)
        .def_static("CenterPoint", &pztopology::TPZTetrahedron::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZTetrahedron& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZTetrahedron::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZTetrahedron::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZTetrahedron::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZTetrahedron::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZTetrahedron::CreateSideIntegrationRule)
    ;

    // TPZPyramid bindings
    py::class_<pztopology::TPZPyramid>(m, "TPZPyramid")
        .def(py::init())
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZPyramid::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZPyramid::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZPyramid::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZPyramid::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZPyramid::SideNodeLocId)
        .def_static("CenterPoint", &pztopology::TPZPyramid::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZPyramid& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZPyramid::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZPyramid::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZPyramid::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZPyramid::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZPyramid::CreateSideIntegrationRule)
    ;

    // TPZPrism bindings
    py::class_<pztopology::TPZPrism>(m, "TPZPrism")
        .def(py::init())
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZPrism::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZPrism::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZPrism::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZPrism::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZPrism::SideNodeLocId)
        .def_static("CenterPoint", &pztopology::TPZPrism::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZPrism& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZPrism::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZPrism::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZPrism::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZPrism::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZPrism::CreateSideIntegrationRule)
    ;

    // TPZCube bindings
    py::class_<pztopology::TPZCube>(m, "TPZCube")
        .def(py::init())
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZCube::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZCube::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZCube::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZCube::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZCube::SideNodeLocId)
        .def_static("CenterPoint", &pztopology::TPZCube::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZCube& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZCube::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZCube::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZCube::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZCube::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZCube::CreateSideIntegrationRule)
    ;

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
