#!/bin/bash

#./install-dependencies.sh

threads=1
BRANCH_GEMS3K=trunk
BuildType=Release
#BuildType=Debug
InstallPrefix=../../Resources
#InstallPrefix=/home/sveta/tmp
workfolder=${PWD}


mkdir -p $InstallPrefix
mkdir -p build
cd build

# Auto-detect platform and use appropriate cmake command
case "$(uname -s)" in
    Darwin)
        echo "Detected macOS - using default cmake generator"
        cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix
        ;;
    MINGW*|MSYS*)
        echo "Detected Windows (MSYS2/MinGW) - using MinGW Makefiles generator"
        cmake -G "MinGW Makefiles" .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix
        ;;
    Linux)
        echo "Detected Linux - using default cmake generator"
        cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix
        ;;
    *)
        echo "Unknown platform: $(uname -s) - trying default cmake generator"
        cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix
        ;;
esac

make -j $threads
make install

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
   sudo ldconfig
fi
