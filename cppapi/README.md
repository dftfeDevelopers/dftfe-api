DFT-FE's cpp API is already in-built into the library in the form of dftfeWrapper class. The doxygen documentation of this class is provided [here](https://dftfedevelopers.github.io/dftfe/classdftfe_1_1dftfe_wrapper.html). 

* The wrapper currently only sets up GGA PBE ground-state DFT calculations using ONCV pseudopotentials (http://www.pseudo-dojo.org for example). All types of boundary conditions (non-periodic, periodic and semi-periodic) are supported by the wrapper. Further collinear spin-polarization DFT calculations are also supported.

* The constructor as well as reinit functions of the wrapper takes the following objects: mpi communicator, list of atomic cartesian coordinates (origin at corner), list of atomic numbers, cell vectors, cell periodic boundary conditions, finite-element mesh size (related to plane-wave kinetic energy cutoff typically used in plane-wave based DFT codes), k-point Monkhorst-Pack grid related inputs, spin polarization toggle (default is off), GPU usage toggle (default if off) and other DFT related input parameters. Please refer to doxygen documentation and the demo example source files for further details on datastructure of the input parameters and other regarding how to use the wrapper.  

* The following two static function calls are mandatory to use the API: *globalHandlesInitialize* just after MPI_Init and *globalHandlesFinalize* just before MPI_Finalize


cpp API demo example
==========================================

* Install DFT-FE current development branch (publicGithubDevelop) from [DFT-FE github repo](https://github.com/dftfeDevelopers/dftfe). The installation instructions for DFT-FE and its dependencies are provied in the development version manual [here](https://github.com/dftfeDevelopers/dftfe/blob/manual/manual-develop.pdf). To use the dftfe library as an API we suggest performing a prefix based installation using CMake. For example, please refer to an example setup script [here](https://github.com/dftfeDevelopers/dftfe/blob/publicGithubDevelop/helpers/NERSCPerlmutterGPU/setupUserPerlmutterPrefixInstall.sh). Please note that two separate libraries will be installed, one for real datatype (Gamma point) and the other for the complex datatype (multiple k-points).

* Choose one of the demo examples by changing the source file in the CMakelists.txt. The description of the demo examples are given inside their respective .cpp files. The only header from the dftfe library which is required to be included in the .cpp file is *dftfeWrapper.h*.

* Create build directory and compile the example using Cmake by passing the dftfe installation path for the complex datatype as that would work all cases. For example on NERSC Perlmutter machine we used:
```
$ mkdir build 
$ cd build
$ cmake -DDFTFE_INSTALL_PATH=/global/common/software/m3916/softwareDFTFE/dftfe/installComplex -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc ../
```

* set *DFTFE_PSP_PATH* environment variable using export. The pseudpotential directory must contain ONCV files in the format: AtomicSymbol.upf

* Run the compiled executable as a MPI job
