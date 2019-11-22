//
// Created by gustavo on 22/11/2019.
//

#ifndef TEMPLATEDCONTAINERSBINDINGS_H
#define TEMPLATEDCONTAINERSBINDINGS_H

#include <pybind11/pybind11.h>
#include <type_traits>

namespace py = pybind11;
using namespace py::literals;

#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

// TPZVec<T> bindings
template<typename T>
void declareTPZVec(py::module& m, const std::string& suffix) {
    using Class = TPZVec<T>;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    PyClass cls(m, ("TPZVec" + suffix).c_str());

    cls.def(py::init());
    cls.def(py::init<int>());
    cls.def(py::init<int64_t, T>());
    cls.def("Resize", [](Class& vec, const int64_t& newsize) { return vec.Resize(newsize); });
    cls.def("Size", [](const Class& vec) { return vec.size(); });
    cls.def("__getitem__",
            [](const Class& vec, int64_t position) {
                if (position >= vec.size() || position < 0) throw py::index_error();
                return vec[position];
            },
            py::is_operator()
    );
    cls.def("__setitem__",
            [](Class& vec, int64_t position, T value) {
                if (position >= vec.size() || position < 0) throw py::index_error();
                vec[position] = value;
            },
            py::is_operator()
    );
    cls.def("__repr__",
            [](Class& cls) {
                std::ostringstream stream;
                cls.Print(stream);
                return stream.str();
            }
    );
}

// TPZManVector<T> bindings
template<typename T>
void declareTPZManVector(py::module& m, const std::string& suffix) {
    using Class = TPZManVector<T>;
    using PyClass = py::class_<Class, std::shared_ptr<Class>, TPZVec<T> >;

    PyClass cls(m, ("TPZManVector" + suffix).c_str());

    cls.def(py::init());
    cls.def(py::init<int>());
    cls.def(py::init<int64_t, T>());
    cls.def("Resize", [](Class& vec, const int64_t& newsize) { return vec.Resize(newsize); });
    cls.def("Size", [](const Class& vec) { return vec.size(); });
    cls.def("__getitem__",
            [](const Class& vec, int64_t position) {
                if (position >= vec.size() || position < 0) throw py::index_error();
                return vec[position];
            },
            py::is_operator()
    );
    cls.def("__setitem__",
            [](Class& vec, int64_t position, T value) {
                if (position >= vec.size() || position < 0) throw py::index_error();
                vec[position] = value;
            },
            py::is_operator()
    );
    cls.def("__repr__",
            [](Class& cls) {
                std::ostringstream stream;
                cls.Print(stream);
                return stream.str();
            }
    );
}

// TPZStack<T> bindings
template<typename T>
void declareTPZStack(py::module& mod, std::string const& suffix) {
    using Class = TPZStack<T>;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    PyClass cls(mod, ("TPZStack" + suffix).c_str());

    cls.def(py::init());
    cls.def(py::init<int, T>());
    cls.def("Push", [](Class& stack, T object) {
        stack.Push(object);
    });
    cls.def("Pop", [](Class& stack) {
        return stack.Pop();
    });
    cls.def("Peek", [](Class& stack) {
        return stack.Peek();
    });
    cls.def("__repr__",
            [](Class& cls) {
                std::ostringstream stream;
                cls.Print(stream);
                return stream.str();
            }
    );
}

// TPZFMatrix<T> bindings
template<typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
void declareTPZFMatrix(py::module& m, std::string const& suffix) {
    using Class = TPZFMatrix<T>;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    PyClass cls(m, ("TPZFMatrix" + suffix).c_str());

    cls.def(py::init());
    cls.def(py::init<int64_t, int64_t>());
    cls.def(py::init<int64_t, int64_t, double>());
    cls.def("GetVal", [](TPZFMatrix<double>& matrix, const int64_t& row, const int64_t& col) {
        if (row >= matrix.Rows() || row < 0) throw py::index_error();
        if (col >= matrix.Cols() || col < 0) throw py::index_error();
        return matrix.GetVal(row, col);
    });
    cls.def("__repr__",
            [](TPZFMatrix<double>& matrix) {
                std::ostringstream stream;
                matrix.Print(stream);
                return stream.str();
            }
    );
}

#endif // TEMPLATEDCONTAINERSBINDINGS_H
