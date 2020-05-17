# segy-stack

A library with C++ and Python API to read 3D Post-stack seismic data in SEG-Y format.

## Installing

``` shell
sudo apt-get install build-essential protobuf-compiler libprotobuf-dev libgflags-dev
git clone --recurse-submodules https://github.com/google/segy-stack.git
cd segy-stack && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## License 

Apache 2.0

This is not an official Google product.
