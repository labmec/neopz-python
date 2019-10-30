#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

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

// Geometric mesh
#include "TPZGmshReader.h"
#include "pzgmesh.h"

//
#include "TPZSSpStructMatrix.h"
#include "pzcmesh.h"


//
#include "TPZMaterial.h"
#include "pzpoisson3d.h"
#include "TPZMatElasticity2D.h"
#include "pzbndcond.h"
#include <map>                    // for map
#include <set>                    // for set
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzmatrix.h"


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
        .def("__repr__",
             [](const TPZManVector<double>& vec) {
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

    // // chunk vector
    //  py::class_<TPZChunkVector<double>>(m, "TPZVecDouble")
    //     .def(py::init())
    //     .def(py::init<int64_t>())
    //     .def(py::init<int64_t, double>())
    //     .def("Resize", [](TPZManVector<double>& vec, const int64_t& newsize) { return vec.Resize(newsize); })
    //     .def("Size", [](const TPZManVector<double>& vec) { return vec.size(); })
    //     .def("__getitem__",
    //          [](const TPZManVector<double>& vec, int64_t position) {
    //              if (position >= vec.size() || position < 0) throw py::index_error();
    //              return vec[position];
    //          },
    //          py::is_operator()
    //     )
    //     .def("__setitem__",
    //          [](TPZManVector<double>& vec, int64_t position, double value) {
    //              if (position >= vec.size() || position < 0) throw py::index_error();
    //              vec[position] = value;
    //          },
    //          py::is_operator()
    //     )
    //     .def("__repr__",
    //          [](const TPZManVector<double>& vec) {
    //              std::string r("TPZVecDouble [");
    //              r += std::to_string(vec[0]);
    //              for (int i = 1; i < vec.NElements(); i++) {
    //                  r += ", ";
    //                  r += std::to_string(vec[i]);
    //              }
    //              r += "]";
    //              return r;
    //          }
    //     )
    // ;
    
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
             }
        )
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
        .def_static("IsInParametricDomain", &pztopology::TPZTriangle::IsInParametricDomain)
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
        .def_static("IsInParametricDomain", &pztopology::TPZQuadrilateral::IsInParametricDomain)
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
        .def_static("IsInParametricDomain", &pztopology::TPZTetrahedron::IsInParametricDomain)
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
        .def_static("IsInParametricDomain", &pztopology::TPZPyramid::IsInParametricDomain)
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
        .def_static("IsInParametricDomain", &pztopology::TPZPrism::IsInParametricDomain)
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
        .def_static("IsInParametricDomain", &pztopology::TPZCube::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZCube::CreateSideIntegrationRule)
    ;

    // TPZGeoMesh bindings
    py::class_<TPZGeoMesh>(m, "TPZGeoMesh")
        .def(py::init())
        .def("Print", [](TPZGeoMesh &GeoMesh){ return GeoMesh.Print();})
        .def("BuildConnectivity", &TPZGeoMesh::BuildConnectivity)
        .def("NElements", &TPZGeoMesh::NElements)
    ;

    // TPZGMshReader
    py::class_<TPZGmshReader>(m, "TPZGmshReader")
        .def(py::init())
        .def("GeometricGmshMesh3", &TPZGmshReader::GeometricGmshMesh3, "Reads geometric mesh from GMsh (3.x) .msh file.") 
        .def("GeometricGmshMesh4", &TPZGmshReader::GeometricGmshMesh4, "Reads geometric mesh from GMsh (4.x) .msh file.") 
    ;
    //TPZMaterial
    py::class_<TPZMaterial, std::unique_ptr<TPZMaterial, py::nodelete>>(m, "TPZMaterial")
//    TPZMaterial.def("CreateBC", [](TPZMaterial *mat, int id, int type, TPZFMatrix<double> &val1, TPZFMatrix<double> &val2  ){ return mat->CreateBC(mat, id, type, val1, val2);})
    .def("CreateBC", &TPZMaterial::CreateBC)
    .def("SetId", &TPZMaterial::SetId)
    
    
    ;
    
//     TPZBndCond *TPZMaterial::CreateBC(TPZMaterial * reference, int id, int typ, TPZFMatrix<STATE> &val1, TPZFMatrix<STATE> &val2)
    
    
    py::class_<TPZBndCond, TPZMaterial , std::unique_ptr<TPZBndCond, py::nodelete>>(m, "TPZBndCond")
    .def(py::init<>())
    .def(py::init<int>())
//    .def(py::init< TPZMaterial *,int ,int , TPZFMatrix<STATE> &,TPZFMatrix<STATE> &>())
     .def(py::init< TPZMaterial * ,int ,int  , TPZFMatrix<STATE> & ,TPZFMatrix<STATE> & >())
    ;
    
    //TPZMaterial
//    py::class_<TPZMatPoisson3d, TPZMaterial >(m, "TPZMatPoisson3d")
    py::class_<TPZMatPoisson3d, TPZMaterial , std::unique_ptr<TPZMatPoisson3d, py::nodelete>>(m, "TPZMatPoisson3d")
        .def(py::init<int, int>())
    
    ;

    py::class_<TPZMatElasticity2D, TPZMaterial >(m, "TPZMatElasticity2D")
            .def(py::init<int>())

            ;

  
    
    //
    py::class_<TPZCompMesh , std::unique_ptr<TPZCompMesh, py::nodelete>>(m, "TPZCompMesh")
//    py::class_<TPZCompMesh >(m, "TPZCompMesh")
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
    .export_values();
    
    py::class_<TPZStructMatrix >(m, "TPZStructMatrix")
    
    
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
    .def("DefineGraphMesh", py::overload_cast<int,const TPZVec<std::string> &,const TPZVec<std::string>&, const std::string &  >(&TPZAnalysis::DefineGraphMesh))
    .def("PostProcess", py::overload_cast<int, int>(&TPZAnalysis::PostProcess))
    
//    void TPZAnalysis::DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile)
    
    
    ;
  // .def("Mult", py::overload_cast<>(&TPZTransform<double>::Mult))
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
