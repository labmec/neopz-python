//
// Created by gustavo on 21/11/2019.
//

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

namespace py = pybind11;
using namespace py::literals;

#ifndef TYPEDPYTHONITERATOR_H
#define TYPEDPYTHONITERATOR_H

template <typename T>
class TypedPythonIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = T*;
    using reference = T&;

    explicit TypedPythonIterator(py::iterator& py_iter) :
            py_iter_(py_iter)
    {}

    explicit TypedPythonIterator(py::iterator&& py_iter) :
            py_iter_(py_iter)
    {}

    value_type operator*() {
        return (*py_iter_).template cast<value_type>();
    }

    TypedPythonIterator operator++(int) {
        auto copy = *this;
        ++py_iter_;
        return copy;
    }

    TypedPythonIterator& operator++() {
        ++py_iter_;
        return *this;
    }

    bool operator!=(TypedPythonIterator &rhs) {
        return py_iter_ != rhs.py_iter_;
    }

    bool operator==(TypedPythonIterator &rhs) {
        return py_iter_ == rhs.py_iter_;
    }

private:
    py::iterator py_iter_;
};

#endif // TYPEDPYTHONITERATOR_H
