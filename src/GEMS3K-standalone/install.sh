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

# For Windows using MinGW use this line
cmake -G "MinGW Makefiles" .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix

# For Mac OS using the clang set use this line
#cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix

# For Mac OS using the GNU compiler set installed via Homebrew, use the next line
#cmake .. -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++ -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX14.5.sdk -DCMAKE_INSTALL_PREFIX=$InstallPrefix 
make -j $threads 
make install

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
   sudo ldconfig
fi
