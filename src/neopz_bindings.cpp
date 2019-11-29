#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace py::literals;

// STL headers
#include <map>
#include <set>

// NeoPZ linear algebra and container classes
#include "pzmatrix.h"
#include "pzstack.h"
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

// Computational mesh
#include "TPZSSpStructMatrix.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZHybridizeHDiv.h"

// Material classes
#include "TPZMaterial.h"
#include "pzpoisson3d.h"
#include "TPZMatElasticity2D.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"

// Python binding utilities
#include "TypedPythonIterator.h"
#include "TemplatedContainersBindings.h"
#include "TemplatedTopologyBindings.h"

PYBIND11_MODULE(neopz, m) {
    m.doc() = R"pbdoc(

    )pbdoc";

    // TPZVec<T> bindings
    declareTPZVec<int>(m, "Int");
    declareTPZVec<REAL>(m, "Real");
    declareTPZVec<std::string>(m, "String");

    // TPZManVector<T> bindings
    declareTPZManVector<int>(m, "Int");
    declareTPZManVector<REAL>(m, "Real");
    declareTPZManVector<std::string>(m, "String");

    // TPZStack<T> bindings
    declareTPZStack<int>(m, "Int");
    declareTPZStack<REAL>(m, "Real");
    declareTPZStack<std::string>(m, "String");

    // TPZFMatrix<T> bindings
    declareTPZFMatrix<int>(m, "Int");
    declareTPZFMatrix<REAL>(m, "Real");

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
                 [](TPZTransform<double>& trans) {
                     std::ostringstream stream;
                     trans.PrintInputForm(stream);
                     return stream.str();
                 }
            );

    // Topology bindings
    declareTopology<pztopology::TPZPoint>(m, "Point");
    declareTopology<pztopology::TPZLine>(m, "Line");
    declareTopology<pztopology::TPZTriangle>(m, "Triangle");
    declareTopology<pztopology::TPZQuadrilateral>(m, "Quadrilateral");
    declareTopology<pztopology::TPZTetrahedron>(m, "Tetrahedron");
    declareTopology<pztopology::TPZPyramid>(m, "Pyramid");
    declareTopology<pztopology::TPZPrism>(m, "Prism");
    declareTopology<pztopology::TPZCube>(m, "Cube");

    // TPZGeoMesh bindings
    py::class_<TPZGeoMesh>(m, "TPZGeoMesh")
            .def(py::init())
            .def("Print", [](TPZGeoMesh& GeoMesh) { return GeoMesh.Print(); })
            .def("BuildConnectivity", &TPZGeoMesh::BuildConnectivity)
            .def("NElements", &TPZGeoMesh::NElements);

    // TPZGMshReader
    py::class_<TPZGmshReader>(m, "TPZGmshReader")
            .def(py::init())
            .def("GeometricGmshMesh3", &TPZGmshReader::GeometricGmshMesh3,
                 "Reads geometric mesh from GMsh (3.x) .msh file.")
            .def("GeometricGmshMesh4", &TPZGmshReader::GeometricGmshMesh4,
                 "Reads geometric mesh from GMsh (4.x) .msh file.");

    // TPZMaterial
    py::class_<TPZMaterial, std::unique_ptr<TPZMaterial, py::nodelete>>(m, "TPZMaterial")
            .def("CreateBC", &TPZMaterial::CreateBC)
            .def("SetId", &TPZMaterial::SetId);

    // TPZBndCond
    py::class_<TPZBndCond, TPZMaterial, std::unique_ptr<TPZBndCond, py::nodelete>>(m, "TPZBndCond")
            .def(py::init<>())
            .def(py::init<int>())
            .def(py::init<TPZMaterial*, int, int, TPZFMatrix<STATE>&, TPZFMatrix<STATE>&>());

    // TPZMatPoisson3d
    py::class_<TPZMatPoisson3d, TPZMaterial, std::unique_ptr<TPZMatPoisson3d, py::nodelete>>(m, "TPZMatPoisson3d")
            .def(py::init<int, int>());

    // TPZMatElasticity2D
    py::class_<TPZMatElasticity2D, TPZMaterial>(m, "TPZMatElasticity2D")
            .def(py::init<int>());

    // TPZCompMesh
    py::class_<TPZCompMesh, std::unique_ptr<TPZCompMesh, py::nodelete>>(m, "TPZCompMesh")
            .def(py::init())
            .def(py::init<TPZGeoMesh*>())
            .def("AutoBuild", [](TPZCompMesh& compmesh) { return compmesh.AutoBuild(); })
            .def("SetDimModel", &TPZCompMesh::SetDimModel)
            .def("InsertMaterialObject",
                 [](TPZCompMesh& compmesh, TPZMaterial* mat) { return compmesh.InsertMaterialObject(mat); })
            .def("SetAllCreateFunctionsContinuous", &TPZCompMesh::SetAllCreateFunctionsContinuous)
            .def("NMaterials", &TPZCompMesh::NMaterials)
            .def("NElements", &TPZCompMesh::NElements)
            .def("Print", [](TPZCompMesh& compmesh) { return compmesh.Print(); })
            .def("SetDefaultOrder", &TPZCompMesh::SetDefaultOrder);


    py::class_<TPZMatrixSolver<STATE> >(m, "TPZMatrixSolver");

    py::class_<TPZStepSolver<STATE>, TPZMatrixSolver<STATE> >(m, "TPZStepSolver")
            .def(py::init())
            .def("SetDirect", &TPZStepSolver<STATE>::SetDirect);

    py::class_<TPZMatrix<double> > mat(m, "TPZMatrix");

    py::enum_<DecomposeType>(m, "DecomposeType")
            .value("ECholesky", DecomposeType::ECholesky)
            .export_values();

    py::class_<TPZStructMatrix>(m, "TPZStructMatrix");

    py::class_<TPZSymetricSpStructMatrix, TPZStructMatrix>(m, "TPZSymetricSpStructMatrix")

            .def(py::init<TPZCompMesh*>())
            .def("Create", [](TPZSymetricSpStructMatrix& spmatrix) {
            })
            .def("SetupMatrixData",
                 [](TPZSymetricSpStructMatrix& spmatrix, TPZStack<int64_t>& elgraph, TPZVec<int64_t>& elgraphindex) {
                     return spmatrix.SetupMatrixData(elgraph, elgraphindex);
                 });

    // TPZAnalysis
    py::class_<TPZAnalysis>(m, "TPZAnalysis")
            .def(py::init())
            .def(py::init<TPZCompMesh*, bool>())
            .def("SetStructuralMatrix", py::overload_cast<TPZStructMatrix&>(&TPZAnalysis::SetStructuralMatrix))
            .def("SetSolver", &TPZAnalysis::SetSolver)
            .def("Assemble", &TPZAnalysis::Assemble)
            .def("Solve", &TPZAnalysis::Solve)
            .def("DefineGraphMesh",
                 py::overload_cast<int, const TPZVec<std::string>&, const TPZVec<std::string>&, const std::string&>(
                         &TPZAnalysis::DefineGraphMesh))
            .def("PostProcess", py::overload_cast<int, int>(&TPZAnalysis::PostProcess))
            ;

    py::class_<TPZMultiphysicsCompMesh>(m, "TPZMultiphysicsCompMesh")
            .def(py::init<TPZGeoMesh *>())
            .def("BuildMultiphysicsSpace", py::overload_cast<TPZVec<int> & , TPZVec<TPZCompMesh * > & >(&TPZMultiphysicsCompMesh::BuildMultiphysicsSpace))
            .def("SetNMeshes", &TPZMultiphysicsCompMesh::SetNMeshes)
            .def("SetAllCreateFunctionsMultiphysicElem", &TPZMultiphysicsCompMesh::SetAllCreateFunctionsMultiphysicElem)
            .def("LoadSolutionFromMeshes", &TPZMultiphysicsCompMesh::LoadSolutionFromMeshes)
            .def("InitializeBlock", &TPZMultiphysicsCompMesh::InitializeBlock)
            ;

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
