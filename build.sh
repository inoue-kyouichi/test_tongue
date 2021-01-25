#!/bin/sh
mkdir build
cd build
cmake -D compiler=intel \
      -D CMAKE_INSTALL_PREFIX=/mnt/c/code/.git/femSolidSolver \
      -D TP_DIR=/home/matsumuramasahiro/lib/TextParser \
      -D EIGEN_DIR=/mnt/c/code/eigen \
      -D enable_GLOG=ON \
      -D GLOG_DIR=/mnt/c/glog-0.3.4 \
      ..

make && make install