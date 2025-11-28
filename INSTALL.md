# INSTALL

## INTRODUCTION

This file describes how to install THAMES-Hydration v5.2.  The instructions in this file
are for the most common use cases, and cover the command line tools.

For further information, or in case of problems, please contact the author,
Jeff Bullard (jwbullard@tamu.edu).

The API documentation of THAMES-Hydration is available in PDF and HTML formats.  It
is not bundled with the software but can be created during the installation
of the software.

For more information about THAMES-Hydration, see the API documentation overview, or
the accompanying user guide (coming soon).  More information about THAMES-Hydration and
its applications can also be found in the following references:

   * Bullard, J.W., Lothenbach, B., Stutzman, P.E., Snyder, K.A., Coupling thermodynamics and digital image models to simulate hydration and microstructure development of portland cement pastes, _Journal of Materials Research_ 26, (2011) 609-622.

   * Feng, P., Miao, C., Bullard, J.W., A model of phase stability, microstructure and properties during leaching of portland cement binders, _Cement and Concrete Composites_ 49, (2014) 9-19.

   * Li, X., Grasley, Z.C., Garboczi, E.J., Bullard, J.W., Modeling the apparent intrinsic viscoelastic relaxation of hydrating cement paste, _Cement and Concrete Composites_ 55, (2014) 322-330.

   * Feng, P., Garboczi, E.J., Miao, C., Bullard, J.W., Microstructural origins of cement paste degradation by external sulfate attack, _Construction and Building Materials_ 96, (2015) 391-403.

   * Li, X., Grasley, Z.C., Bullard, J.W., Garboczi, E.J., Computing the time evolution of the apparent viscoelastic/viscoplastic Poisson's ratio of hydrating cement paste, _Cement and Concrete Composites_, 56 (2015) 121-133.


## PREREQUISITES

* One of the following compiler collections:
    * Clang 17.0 or newer
    * GNU Compiler Collection (gcc/g++ 14.0 or newer)

* CMake (>= 3.30), the build system used by THAMES-Hydration
 * Required for building THAMES-Hydration

* Doxygen (>= 1.8.13), the API documentation software
 * Required for creating the API documentation

* (Optional) LaTeX 2e, the document preparation system
 * Required only for creating the PDF version of the API documentation

* GEM-Selektor (optional)
 * Needed for producing input data files for thermodynamic calculations.
 * Source code available at https://bitbucket.org/gems4/gems3gui.git
 * Binaries available at https://gems.web.psi.ch

## Building on Mac OS

Modern Mac OS uses the Apple Clang compiler by default, and it seems to work for THAMES-Hydration without further modifications.

### Build GEMS3K library

* cd /PathToTHAMES-Hydration/THAMES-Hydration/src/GEMS3K-standalone
* Open the file `install.sh` with your favorite text editor.


#### Clang Compiler from Apple
No changes are needed to `install.sh` so you can
close it.

#### GNU Compiler via Homebrew

* Comment out the line that begins with `cmake ..` by placing a `#` in the first column.
* Uncomment the line directly below that one; it should look something like this:

```
cmake .. -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++ -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX14.5.sdk -DCMAKE_INSTALL_PREFIX=$InstallPrefix 
```

* Save the `install.sh` file and close it
* Run the command `./install.sh`

### Build THAMES-Hydration
Next, build and install THAMES-Hydration. The recommended way to configure THAMES-Hydration is to do an out-of-source
build, which means that the original files and directories are left untouched.
Doing this makes the re-compiling and cleaning of the installation files
much simpler.

* cd /PathToTHAMES-Hydration/THAMES-Hydration/build
* Run one of the following commands, depending on your C++ compiler

#### Apple's Clang Compiler

```
cmake ..
```

#### GNU Compiler via Homebrew

```
cmake -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++ -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX14.5.sdk ..
```

* `make`
* `make install`

This will install the "thames" executable in the /PathToTHAMES-Hydration/bin directory, and the
static libraries in the PathToTHAMES-Hydration/lib directory.


## Building on Unix or Linux

### Build GEMS3K library

1. cd /PathToTHAMES-Hydration/THAMES-Hydration/src/GEMS3K-standalone
2. ./install.sh

### Build THAMES-Hydration
Next, build and install THAMES-Hydration. The recommended way to configure THAMES-Hydration is to do an out-of-source
build, which means that the original files and directories are left untouched.
Doing this makes the re-compiling and cleaning of the installation files
much simpler.

1. cd /PathToTHAMES-Hydration/THAMES-Hydration/build
2. `cmake ..`
3. `make`
4. `make install`
5. `make doc`

This will install the "thames" executable in the /PathToTHAMES-Hydration/bin directory, and the
static libraries in the PathToTHAMES-Hydration/lib directory.

## Building on Windows

Windows does not come prebuilt with any kind
of system for compiling C/C++ code, but it they have made it easier to get a full Linux environment set up, which is part of these instructions.  Therefore

### Install WSL

- Switch to admin user
- Type `powershell` in the search bar.
- In the powershell window, type `wsl --install`
- Reboot computer, select userid and password

### Install GEMS3K Library

- Start up Ubuntu by typing `ubuntu` in the search bar.
- sudo apt update
- sudo apt install git
- sudo apt install cmake
- sudo apt install build-essential
- git clone https://github.com/jwbullard/THAMES-Hydration.git
- Wait for the THAMES-Hydration folder to fully download
- In the Windows file explorer, navigate to Linux/home/THAMES-Hydration/src/GEMS3K-standalone
- Open the file install.h with Notepad
- Comment out the line that begins with `cmake ..` by placing a `#` in the first column.
- Add a line directly below it that reads

```
cmake .. -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../../Resources 
```
- Save the `install.sh` file and close the file explorer
- Back in the powershell, `cd THAMES-Hydration/src/GEMS3K-standalone`
- ./install.sh
- cd ../../build
- Run this command:
```
cmake .. -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../bin 
```
- make
- make install

This procedure should install the `thames` and `viz` executables in the Linux/home/THAMES-Hydration/bin folder

## UNINSTALLING

To uninstall everything except the original source files, simply delete
recursively everything in the build/ directory.
