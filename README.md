# segy-stack

A library with C++ and Python API to read 3D Post-stack seismic data in SEG-Y format.

## Installing

``` shell
sudo apt-get install build-essential protobuf-compiler 
git clone --recurse-submodules sso://user/wrlewis/segy-stack
cd segy-stack && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
