# neopz-python
This project implements Python bindings for [NeoPZ](https://github.com/labmec/neopz) finite element library using [pybind11](https://github.com/pybind/pybind11).
Make sure you have both libraries (and CMake!) installed in your environment.

## Requirements

### 1. Pybind11
Download Pybind11-v2.4.3 release:
```
$ wget https://github.com/pybind/pybind11/archive/v2.4.3.zip
```
Extract the source files and create pybind11-2.4.3-build directory:
```
$ unzip v2.4.3.zip
$ mkdir pybind11-2.4.3-build
```
Enter in pybind11-2.4.3 directory:
```
$ cd pybind11-2.4.3-build
```
Generate the Makefile with CMake and install Pybind11:
```
$ cmake -DPYBIND11_TEST=off ../pybind11-2.4.3
$ sudo make install
```

### 2. NeoPZ

The tutorial to install NeoPZ can be found at http://www.labmec.org.br/wiki/howto/roteiro_pzlinux_eng. However, an additional step has to be executed to create PZPy module. After cloning the source files in item 3.3 of NeoPZ install tutorial, execute:

```
$ cd neopz
$ git checkout neopz-for-python
$ cd ..
```
then, follow the next steps of the tutorial.


## Installing the module

Download the latest source tree of neopz-python:
```
$ git clone https://github.com/labmec/neopz-python
```
Enter in neopz-python directory and generate de module:
```
$ cd neopz-python
$ python3 setup.py develop
```
(You'll may need `sudo` for this depending on where `python3` is installed in your machine.)
