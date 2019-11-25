//
// Created by gustavo on 25/11/2019.
//

#ifndef TEMPLATEDTOPOLOGYBINDINGS_H
#define TEMPLATEDTOPOLOGYBINDINGS_H

#include <pybind11/pybind11.h>
#include <type_traits>

namespace py = pybind11;
using namespace py::literals;

#include "tpzpoint.h"
#include "tpzpoint.h"
#include "tpztriangle.h"

// Templated bindings for different NeoPZ topology classes

template<typename Class>
void declareTopology(py::module& m, const std::string& suffix) {

    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    PyClass cls(m, ("TPZ" + suffix).c_str());

    cls.def(py::init());
    cls.def_static("SideDimension", &Class::SideDimension);
    cls.def_static("LowerDimensionSides",
                   py::overload_cast<int, TPZStack<int>&>(&Class::LowerDimensionSides));
    cls.def_static("LowerDimensionSides",
                   py::overload_cast<int, TPZStack<int>&, int>(&Class::LowerDimensionSides));
    cls.def_static("HigherDimensionSides", &Class::HigherDimensionSides);
    cls.def_static("NSideNodes", &Class::NSideNodes);
    cls.def_static("SideNodeLocId", &Class::SideNodeLocId);
    cls.def_static("NumSides", py::overload_cast<>(&Class::NumSides));
    cls.def_static("CenterPoint", &Class::CenterPoint);
    cls.def_static("RefElVolume", [](Class& topology) { return topology.RefElVolume(); });
    cls.def_static("SideToSideTransform", &Class::SideToSideTransform);
    cls.def_static("TransformSideToElement", &Class::TransformSideToElement);
    cls.def_static("TransformElementToSide", &Class::TransformElementToSide);
    cls.def_static("IsInParametricDomain",
                   (bool (*)(const TPZVec<REAL>&, REAL)) &Class::IsInParametricDomain);
    cls.def_static("CreateSideIntegrationRule", &Class::CreateSideIntegrationRule);
}

#endif //TEMPLATEDTOPOLOGYBINDINGS_H
