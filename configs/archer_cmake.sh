#!/bin/bash

cmake . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/work/e281/e281/jomotani/soltransport/BOUT-configs/archer/BOUT-dev/build -DCMAKE_CXX_COMPILER=CC
cmake --build build -j 32
