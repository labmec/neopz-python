#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

// cpp STL
#include <map>                    // for map
#include <set>                    // for set

namespace py = pybind11;
using namespace py::literals;

// NeoPZ container classes
#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pztrnsform.h"

// NeoPZ topology classes
#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzpyramid.h"
#include "tpzprism.h"
#include "tpzcube.h"

// Mesh
#include "pzgmesh.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeoelbc.h"
#include "pzgeoelside.h"

#include "pzcmesh.h"

// Matrix
#include "pzmatrix.h"

//StrMatrix
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZStructMatrixBase.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

// Material
#include "TPZMaterial.h"
#include "pzpoisson3d.h"
#include "TPZMatElasticity2D.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzpostprocanalysis.h"
#include "TPZMatLaplacian.h"

#include "pzsolve.h"
#include "tpzautopointer.h"

// Elastoplasticity 
#include "TPZBndCondWithMem.h" 
#include "TPZBndCondWithMem_impl.h"
#include "TPZMatElastoPlastic2D.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElasticResponse.h"
#include "TPZElastoPlasticMem.h"
#include "TPZPlasticCriterion.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"
#include "TPZBndCondWithMem.h"

// For SBFEM simulations
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"

PYBIND11_MAKE_OPAQUE(std::map<int,int>)

PYBIND11_MODULE(neopz, m) {
    m.doc() = R"pbdoc(
        -------------------------
        Python bindings for NeoPZ
        -------------------------
    )pbdoc";

    // TPZVec<double> bindings
    py::class_<TPZVec <double>>(m, "TPZVecDouble")
        .def(py::init())
        .def(py::init<int64_t>())
        .def(py::init<int64_t, double>())
        .def("Resize", [](TPZVec<double>& vec, const int64_t& newsize) { return vec.Resize(newsize); })
        .def("Size", [](const TPZVec<double>& vec) { return vec.size(); })
        .def("__getitem__",
             [](const TPZVec<double>& vec, int64_t position) {
                 if (position >= vec.size() || position < 0) throw py::index_error();
                 return vec[position];
             },
             py::is_operator()
        )
        .def("__setitem__",
             [](TPZVec<double>& vec, int64_t position, double value) {
                 if (position >= vec.size() || position < 0) throw py::index_error();
                 vec[position] = value;
             },
             py::is_operator()
        )
        .def("__repr__",
             [](const TPZVec<double>& vec) {
                 std::string r("TPZVecDouble [");
                 r += std::to_string(vec[0]);
                 for (int i = 1; i < vec.NElements(); i++) {
                     r += ", ";
                     r += std::to_string(vec[i]);
                 }
                 r += "]";
                 return r;
             }
        )
    ;

    // TPZManVector<double> bindings
    py::class_<TPZManVector<double>, TPZVec<double>>(m, "TPZManVecDouble")
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
        .def("__repr__",
             [](const TPZManVector<double>& vec) {
                 std::string r("TPZManVecDouble [");
                 r += std::to_string(vec[0]);
                 for (int i = 1; i < vec.NElements(); i++) {
                     r += ", ";
                     r += std::to_string(vec[i]);
                 }
                 r += "]";
                 return r;
             }
        )
    ;

    // TPZVec<int64_t> bindings
    py::class_<TPZVec <int>>(m, "TPZVecInt")
        .def(py::init())
        .def(py::init<int64_t>())
        .def(py::init<int64_t, int64_t>())
        .def("Resize", [](TPZVec<int64_t>& vec, const int64_t& newsize) { return vec.Resize(newsize); })
        .def("Size", [](const TPZVec<int64_t>& vec) { return vec.size(); })
        .def("__getitem__",
             [](const TPZVec<int64_t>& vec, int64_t position) {
                 if (position >= vec.size() || position < 0) throw py::index_error();
                 return vec[position];
             },
             py::is_operator()
        )
        .def("__setitem__",
             [](TPZVec<int64_t>& vec, int64_t position, int64_t value) {
                 if (position >= vec.size() || position < 0) throw py::index_error();
                 vec[position] = value;
             },
             py::is_operator()
        )
        .def("__repr__",
             [](const TPZVec<int64_t>& vec) {
                 std::string r("TPZVecInt [");
                 r += std::to_string(vec[0]);
                 for (int i = 1; i < vec.NElements(); i++) {
                     r += ", ";
                     r += std::to_string(vec[i]);
                 }
                 r += "]";
                 return r;
             }
        )
    ;

    // TPZManVector<int> bindings
    py::class_<TPZManVector<int>, TPZVec<int>>(m, "TPZManVecInt")
        .def(py::init())
        .def(py::init<int64_t>())
        .def(py::init<int64_t, int64_t>())
        .def("Resize", [](TPZManVector<int64_t>& vec, const int64_t& newsize) { return vec.Resize(newsize); })
        .def("Size", [](const TPZManVector<int64_t>& vec) { return vec.size(); })
        .def("__getitem__",
             [](const TPZManVector<int64_t>& vec, int64_t position) {
                 if (position >= vec.size() || position < 0) throw py::index_error();
                 return vec[position];
             },
             py::is_operator()
        )
        .def("__setitem__",
             [](TPZManVector<int64_t>& vec, int64_t position, double value) {
                 if (position >= vec.size() || position < 0) throw py::index_error();
                 vec[position] = value;
             },
             py::is_operator()
        )
        .def("__repr__",
             [](const TPZManVector<int64_t>& vec) {
                 std::string r("TPZManVecInt [");
                 r += std::to_string(vec[0]);
                 for (int i = 1; i < vec.NElements(); i++) {
                     r += ", ";
                     r += std::to_string(vec[i]);
                 }
                 r += "]";
                 return r;
             }
        )
    ;
    
    // TPZFMatrix<double> bindings
    py::class_<TPZFMatrix<double>>(m, "TPZFMatrix")
        .def(py::init())
        .def(py::init<int64_t, int64_t>())
        .def(py::init<int64_t, int64_t, double>())
        .def("GetVal", [](TPZFMatrix<double>& matrix, const int64_t& row, const int64_t& col) {
            if (row >= matrix.Rows() || row < 0) throw py::index_error();
            if (col >= matrix.Cols() || col < 0) throw py::index_error();
            return matrix.GetVal(row, col);
        })
        .def("SetItem",
             [](TPZFMatrix<double>& mat, int64_t rows, int64_t cols, double value) {
                 if (rows >= mat.Rows() || rows < 0) throw py::index_error();
                 if (cols >= mat.Cols() || cols < 0) throw py::index_error();
                 mat(rows,cols) = value;
             }//,             py::is_operator()
        )
        .def("__repr__",
             [](TPZFMatrix<double>& matrix) {
                 std::string r("TPZFMatrix ");
                 r += "'(";
                 r += std::to_string(matrix.Rows());
                 r += " x ";
                 r += std::to_string(matrix.Cols());
                 r += ")' = [\n";
                 for (int64_t row = 0; row < matrix.Rows(); row++) {
                     r += "\t";
                     for (int64_t col = 0; col < matrix.Cols(); col++) {
                         r += std::to_string(matrix.GetVal(row, col));
                         r += "  ";
                     }
                     r += "\n";
                 }
                 r += "]";
                 return r;
             })
        .def("Resize", & TPZFMatrix<double>::Resize)
        .def("Zero", &TPZFMatrix<double>::Zero)
    ;
    
    // TPZStack<int> bindings
    py::class_<TPZStack<int>>(m, "TPZStackInt")
        .def(py::init())
        .def(py::init<int, int>())
        .def("Push", [](TPZStack<int>& stack, const int object) {
            stack.Push(object);
        })
        .def("Pop", [](TPZStack<int>& stack) {
            return stack.Pop();
        })
        .def("Peek", [](TPZStack<int>& stack) {
            return stack.Peek();
        })
        .def("__repr__",
            [](const TPZStack<int>& stack) {
                std::string r("TPZStackInt [");
                for (int i = 0; i < stack.NElements(); i++) {
                    r += std::to_string(stack[i]);
                    if (i != stack.NElements() - 1) {
                        r += ", ";
                    }
                }
                r += "]";
                return r;
            }
        )
    ;
    
    // TPZAdmChunkVectorNodes bindings
    py::class_<TPZAdmChunkVector<TPZGeoNode>>(m, "TPZAdmChunkVectorNodes")
        .def(py::init<int>())
        .def("NElements", &TPZAdmChunkVector<TPZGeoNode>::NElements)
        .def("Resize", &TPZAdmChunkVector<TPZGeoNode>::Resize)
        .def("__getitem__",
             [](const TPZAdmChunkVector<TPZGeoNode>& vec, int64_t position) {
                 if (position >= vec.NElements() || position < 0) throw py::index_error();
                 return vec[position];
             },
             py::is_operator()
        )
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
        .def("__repr__",
             [](const TPZTransform<double>& trans) {
                //  std::stringstream repr;
                //  trans.Print(repr);
                //  return repr.str();
                 std::string r("TPZTransform\n");
                 r += "  Mult ";
                 r += "'(";
                 r += std::to_string(trans.Mult().Rows());
                 r += " x ";
                 r += std::to_string(trans.Mult().Cols());
                 r += ")' = [\n";
                 for (int64_t row = 0; row < trans.Mult().Rows(); row++) {
                     r += "\t";
                     for (int64_t col = 0; col < trans.Mult().Cols(); col++) {
                         r += std::to_string(trans.Mult().GetVal(row, col));
                         r += "  ";
                     }
                     r += "\n";
                 }
                 r += "]\n";
                 r += "  Sum";
                 r += "'(";
                 r += std::to_string(trans.Sum().Rows());
                 r += " x ";
                 r += std::to_string(trans.Sum().Cols());
                 r += ")' = [\n";
                 for (int64_t row = 0; row < trans.Sum().Rows(); row++) {
                     r += "\t";
                     for (int64_t col = 0; col < trans.Sum().Cols(); col++) {
                         r += std::to_string(trans.Sum().GetVal(row, col));
                         r += "  ";
                     }
                     r += "\n";
                 }
                 r += "]";
                 return r;
             }
        )
    ;

    // TPZPoint bindings
    py::class_<pztopology::TPZPoint>(m, "TPZPoint")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZPoint::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZPoint::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZPoint::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZPoint::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZPoint::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZPoint::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZPoint::NumSides))
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
        .def_static("SideDimension", &pztopology::TPZLine::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZLine::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZLine::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZLine::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZLine::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZLine::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZLine::NumSides))
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
        .def_static("SideDimension", &pztopology::TPZTriangle::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZTriangle::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZTriangle::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZTriangle::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZTriangle::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZTriangle::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZTriangle::NumSides))
        .def_static("CenterPoint", &pztopology::TPZTriangle::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZTriangle& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZTriangle::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZTriangle::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZTriangle::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZTriangle::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZTriangle::CreateSideIntegrationRule)
    ;

    // TPZQuadrilateral bindings
    py::class_<pztopology::TPZQuadrilateral>(m, "TPZQuadrilateral")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZQuadrilateral::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZQuadrilateral::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZQuadrilateral::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZQuadrilateral::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZQuadrilateral::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZQuadrilateral::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZQuadrilateral::NumSides))
        .def_static("CenterPoint", &pztopology::TPZQuadrilateral::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZQuadrilateral& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZQuadrilateral::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZQuadrilateral::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZQuadrilateral::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZQuadrilateral::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZQuadrilateral::CreateSideIntegrationRule)
    ;

    // TPZTetrahedron bindings
    py::class_<pztopology::TPZTetrahedron>(m, "TPZTetrahedron")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZTetrahedron::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZTetrahedron::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZTetrahedron::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZTetrahedron::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZTetrahedron::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZTetrahedron::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZTetrahedron::NumSides))
        .def_static("CenterPoint", &pztopology::TPZTetrahedron::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZTetrahedron& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZTetrahedron::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZTetrahedron::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZTetrahedron::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZTetrahedron::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZTetrahedron::CreateSideIntegrationRule)
    ;

    // TPZPyramid bindings
    py::class_<pztopology::TPZPyramid>(m, "TPZPyramid")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZPyramid::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZPyramid::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZPyramid::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZPyramid::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZPyramid::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZPyramid::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZPyramid::NumSides))
        .def_static("CenterPoint", &pztopology::TPZPyramid::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZPyramid& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZPyramid::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZPyramid::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZPyramid::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZPyramid::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZPyramid::CreateSideIntegrationRule)
    ;

    // TPZPrism bindings
    py::class_<pztopology::TPZPrism>(m, "TPZPrism")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZPrism::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZPrism::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZPrism::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZPrism::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZPrism::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZPrism::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZPrism::NumSides))
        .def_static("CenterPoint", &pztopology::TPZPrism::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZPrism& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZPrism::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZPrism::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZPrism::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZPrism::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZPrism::CreateSideIntegrationRule)
    ;

    // TPZCube bindings
    py::class_<pztopology::TPZCube>(m, "TPZCube")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZCube::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZCube::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZCube::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZCube::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZCube::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZCube::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZCube::NumSides))
        .def_static("CenterPoint", &pztopology::TPZCube::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZCube& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZCube::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZCube::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZCube::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZCube::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZCube::CreateSideIntegrationRule)
    ;

    // TPZGeoMesh bindings
    py::class_<TPZGeoMesh>(m, "TPZGeoMesh")
        .def(py::init())
        .def("Print", [](TPZGeoMesh &GeoMesh){ return GeoMesh.Print();})
        .def("BuildConnectivity", &TPZGeoMesh::BuildConnectivity)
        .def("NElements", &TPZGeoMesh::NElements)
        .def("Dimension", &TPZGeoMesh::Dimension)
        .def("NNodes", &TPZGeoMesh::NNodes)
        .def("NodeVec", py::overload_cast<>(&TPZGeoMesh::NodeVec))
        .def("NElements", &TPZGeoMesh::NElements)
        .def("Element", &TPZGeoMesh::Element)
    ;

    // TPZGeoNode bindings
    py::class_<TPZGeoNode>(m, "TPZGeoNode")
        .def(py::init())
        .def("GetCoordinates", &TPZGeoNode::GetCoordinates)
    ;

    // TPZGeoEl
    py::class_<TPZGeoEl, std::unique_ptr<TPZGeoEl, py::nodelete>  >(m, "TPZGeoEl")
        .def("NSides", &TPZGeoEl::NSides)
        .def("NSideNodes", &TPZGeoEl::NSideNodes)
        .def("SideNodeIndex", &TPZGeoEl::SideNodeIndex)
        .def("SideDimension", &TPZGeoEl::SideDimension)
        .def("Dimension", &TPZGeoEl::Dimension)
        .def("Neighbour", &TPZGeoEl::Neighbour)
    ;

    // TPZGeoElBC
    py::class_<TPZGeoElBC >(m, "TPZGeoElBC")
        .def(py::init<TPZGeoEl*, int, int>())
    ;

    // TPZGeoElSide
    py::class_<TPZGeoElSide, std::unique_ptr<TPZGeoElSide, py::nodelete> >(m, "TPZGeoElSide")
        .def(py::init())
        .def(py::init<TPZGeoEl*, int>())
        .def("Neighbour", &TPZGeoElSide::Neighbour)
    ;

    // TPZGMshReader
    py::class_<TPZGmshReader>(m, "TPZGmshReader")
        .def(py::init())
        .def("GeometricGmshMesh3", &TPZGmshReader::GeometricGmshMesh3, "Reads geometric mesh from GMsh (3.x) .msh file.") 
        .def("GeometricGmshMesh4", &TPZGmshReader::GeometricGmshMesh4, "Reads geometric mesh from GMsh (4.x) .msh file.") 
    ;

    //TPZMaterial
    py::class_<TPZMaterial, std::unique_ptr<TPZMaterial, py::nodelete>>(m, "TPZMaterial")
        .def("CreateBC", &TPZMaterial::CreateBC)
        .def("SetId", &TPZMaterial::SetId)
        .def("Id", &TPZMaterial::Id)
    ;

    py::class_<TPZBndCond, TPZMaterial , std::unique_ptr<TPZBndCond, py::nodelete>>(m, "TPZBndCond")
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init< TPZMaterial * ,int ,int  , TPZFMatrix<STATE> & ,TPZFMatrix<STATE> & >())
    ;
    
    py::class_<TPZMatLaplacian, TPZMaterial, std::unique_ptr<TPZMatLaplacian, py::nodelete>>(m, "TPZMatLaplacian")
        .def(py::init<int,int>())
    ;

    py::class_<TPZMatPoisson3d, TPZMaterial, std::unique_ptr<TPZMatPoisson3d, py::nodelete>>(m, "TPZMatPoisson3d")
        .def(py::init<int, int>())
    ;

    py::class_<TPZMatElasticity2D, TPZMaterial >(m, "TPZMatElasticity2D")
        .def(py::init<int>())
    ;
    
    py::class_<TPZElasticity3D, TPZMaterial, std::unique_ptr<TPZElasticity3D, py::nodelete>>(m, "TPZElasticity3D")
        .def(py::init<int>())
        .def("SetMaterialDataHook", & TPZElasticity3D::SetMaterialDataHook)
    ;

    //
    py::class_<TPZCompMesh , std::unique_ptr<TPZCompMesh, py::nodelete>>(m, "TPZCompMesh")
        .def(py::init())
        .def(py::init<TPZGeoMesh *>())
        .def("AutoBuild", [](TPZCompMesh &compmesh){ return compmesh.AutoBuild();})
        .def("SetDimModel", &TPZCompMesh::SetDimModel )
        .def("InsertMaterialObject", [](TPZCompMesh &compmesh, TPZMaterial *mat){ return compmesh.InsertMaterialObject(mat);} )
        .def("SetAllCreateFunctionsContinuous", &TPZCompMesh::SetAllCreateFunctionsContinuous )
        .def("NMaterials", &TPZCompMesh::NMaterials )
        .def("NElements", &TPZCompMesh::NElements)
        .def("Print", [](TPZCompMesh &compmesh){ return compmesh.Print();})
        .def("SetDefaultOrder",&TPZCompMesh::SetDefaultOrder)
        .def("FindMaterial", &TPZCompMesh::FindMaterial)
        .def("NEquations", &TPZCompMesh::NEquations)
        .def("__repr__",
             [](TPZCompMesh & comp) {
                 std::ofstream printstream;
                 comp.Print(printstream);
                 return printstream;
             }
        )
    ;
    

    py::class_<TPZMatrixSolver<STATE> >(m, "TPZMatrixSolver")
    ;
    
    py::class_<TPZStepSolver<STATE>, TPZMatrixSolver<STATE> >(m, "TPZStepSolver")
        .def(py::init())
        .def("SetDirect", &TPZStepSolver<STATE>::SetDirect)
    ;
    
    py::class_<TPZMatrix<double> > mat(m, "TPZMatrix")
    ;
    
    py::enum_<DecomposeType>(m, "DecomposeType")
        .value("ECholesky", DecomposeType::ECholesky)
        .export_values()
    ;

    py::class_<TPZStructMatrix >(m, "TPZStructMatrix")
    ;

    py::class_<TPZFStructMatrix,TPZStructMatrix >(m, "TPZFStructMatrix")

        .def(py::init<TPZCompMesh *>())
    ;

    py::class_<TPZSymetricSpStructMatrix, TPZStructMatrix >(m, "TPZSymetricSpStructMatrix")
        .def(py::init<TPZCompMesh *>())
        .def("Create", [](TPZSymetricSpStructMatrix & spmatrix) {
        })
        .def("SetupMatrixData", [](TPZSymetricSpStructMatrix & spmatrix,TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex) {
            return spmatrix.SetupMatrixData(elgraph, elgraphindex);
    })
    
    ;
//
    py::class_<TPZVec<std::string> >(m, "TPZVecString")
        .def(py::init())
        .def(py::init<int64_t>())
        .def(py::init<int64_t, std::string>())
        .def("__getitem__",
             [](const TPZVec<std::string>& vec, int64_t position) {
                 if (position >= vec.size() || position < 0) throw py::index_error();
                 return vec[position];
             },
             py::is_operator()
             )
        .def("__setitem__",
             [](TPZVec<std::string>& vec, int64_t position, std::string value) {
                 if (position >= vec.size() || position < 0) throw py::index_error();
                 vec[position] = value;
             },
             py::is_operator()
             )
    ;
    //
    
    py::class_<TPZAnalysis >(m, "TPZAnalysis")
        .def(py::init())
        .def(py::init<TPZCompMesh *, bool>())
        .def("SetStructuralMatrix",py::overload_cast<TPZStructMatrix &>(&TPZAnalysis::SetStructuralMatrix))
        .def("SetSolver", &TPZAnalysis::SetSolver)
        .def("Assemble", &TPZAnalysis::Assemble)
        .def("Solve", &TPZAnalysis::Solve)
        .def("Solution", &TPZAnalysis::Solution)
        .def("DefineGraphMesh", py::overload_cast<int,const TPZVec<std::string> &,const TPZVec<std::string>&, const std::string &  >(&TPZAnalysis::DefineGraphMesh))
        .def("DefineGraphMesh", py::overload_cast<int,const TPZVec<std::string> &,const TPZVec<std::string>&, const TPZVec<std::string>&, const std::string &  >(&TPZAnalysis::DefineGraphMesh))
        .def("PostProcess", py::overload_cast<int, int>(&TPZAnalysis::PostProcess))
        .def("NormRhs", &TPZAnalysis::NormRhs)
        .def("Mesh", &TPZAnalysis::Mesh)
        .def("SetExact", &TPZAnalysis::SetExact)
        .def("SetExact2", [](TPZAnalysis &self){
            std::function<void(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)>
                Sol_exact = [](const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
                    dsol.Resize(3,3);
                    sol.resize(3);
                    dsol.Zero();
                    
                    REAL xv = x[0];
                    REAL yv = x[1];
                //  REAL zv = x[2];
                    
                    REAL v_x =  cos(2*M_PI*yv)*sin(2*M_PI*xv);
                    REAL v_y =  -(cos(2*M_PI*xv)*sin(2*M_PI*yv));
                    
                    sol[0]=v_x;
                    sol[1]=v_y;
                    
                    // vx direction
                    dsol(0,0)= 2*M_PI*cos(2*M_PI*xv)*cos(2*M_PI*yv);
                    dsol(0,1)= 2*M_PI*sin(2*M_PI*xv)*sin(2*M_PI*yv);

                    
                    // vy direction
                    dsol(1,0)= -2*M_PI*sin(2*M_PI*xv)*sin(2*M_PI*yv);
                    dsol(1,1)= -2*M_PI*cos(2*M_PI*xv)*cos(2*M_PI*yv);
                    
                };

                self.SetExact(Sol_exact);
                self.Mesh()->ElementSolution().Redim(self.Mesh()->NElements(),9);
            return 0;
        })
        // .def("SetForcingFunction", &TPZAnalysis::SetForcingFunction)
    ;

    py::class_<TPZPostProcAnalysis, TPZAnalysis> (m,"TPZPostProcAnalysis")
        .def(py::init())
        .def("SetCompMesh", &TPZPostProcAnalysis::SetCompMesh)
        .def("SetPostProcessVariables", &TPZPostProcAnalysis::SetPostProcessVariables)
        .def("SetStructuralMatrix",py::overload_cast<TPZStructMatrix &>(&TPZAnalysis::SetStructuralMatrix))
        .def("TransferSolution", &TPZPostProcAnalysis::TransferSolution)
        // .def("DefineGraphMesh", py::overload_cast<int,const TPZVec<std::string> &,const TPZVec<std::string>&, const TPZVec<std::string>&, const std::string &  >(&TPZPostProcAnalysis::DefineGraphMesh))
        // .def("PostProcess", py::overload_cast<int, int>(&TPZPostProcAnalysis::PostProcess))
    ;

    py::class_<TPZSBFemVolume, std::unique_ptr<TPZSBFemVolume, py::nodelete> >(m, "TPZSBFemVolume")
        .def(py::init())
        .def("ReadUNSWSBGeoFile", &TPZSBFemVolume::ReadUNSWSBGeoFile)
    ;

    // TPZGeoMesh bindings
    py::class_<TPZVTKGeoMesh, std::unique_ptr<TPZVTKGeoMesh, py::nodelete>>(m, "TPZVTKGeoMesh")
        .def(py::init())
        .def_static("PrintGMeshVTK",  py::overload_cast<TPZGeoMesh*, const char *, int>(&TPZVTKGeoMesh::PrintGMeshVTK))
    ;

    py::bind_map<std::map<int, int>>(m, "MapIntInt");

    // TPZBuildSBFem bindings
    py::class_<TPZBuildSBFem, std::unique_ptr<TPZBuildSBFem, py::nodelete>>(m, "TPZBuildSBFem")
        .def(py::init<TPZGeoMesh*, int, std::map<int,int> &>())
        .def("SetPartitions", &TPZBuildSBFem::SetPartitions)
        .def("BuildComputationalMeshFromSkeleton", &TPZBuildSBFem::BuildComputationalMeshFromSkeleton)
    ;
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
