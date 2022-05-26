DFT-FE's cpp API is already in-built into the library in the form of dftfeWrapper class. The doxygen documentation of this class is provided [here](https://dftfedevelopers.github.io/dftfe/classdftfe_1_1dftfe_wrapper.html). 

* The wrapper currently only sets up GGA PBE ground-state DFT calculations using ONCV pseudopotentials (http://www.pseudo-dojo.org for example) on either CPU or GPUs. All types of boundary conditions (non-periodic, periodic and semi-periodic) are supported by the wrapper. Further collinear spin-polarization DFT calculations are also supported. Currently, the wrapper uses Anderson mixing strategy for preconditioning of the SCF fixed point iteration with a default mixing parameter of 0.2, which can be modified via the wrapper. More functionality and more efficient black-box mixing strategies will be added to the API soon.

* We use atomic units throughout the wrapper.

* The constructor as well as reinit functions of the wrapper takes the following objects: mpi communicator, list of atomic cartesian coordinates (origin at corner), list of atomic numbers, cell vectors, cell periodic boundary conditions, finite-element mesh size (related to plane-wave kinetic energy cutoff typically used in plane-wave based DFT codes), k-point Monkhorst-Pack grid related inputs, spin polarization toggle (default is off), GPU usage toggle (default if off), mixing parameter and other DFT related input parameters. Please refer to doxygen documentation and the demo example source files for further details on datastructure of the input parameters and other regarding how to use the wrapper. 

* Geometry update calls to update the atomic coordinates (`updateAtomPositions`) and/or the cell vectors (`deformCell`) are also provided which perform a fast update of the related dftfe data structures without deleting the object and performing initialization from scratch. Further, this allows the electronic fields to be reused as optimal initial guesses for subsequent ground-state solve calls.

* dftfeWrapper object can be deleted using the `clear` member function. 

* The following two static member function calls of dftfeWrapper class are mandatory to use the API: `globalHandlesInitialize` just after MPI_Init and `globalHandlesFinalize` just before MPI_Finalize and after deleting all dftfeWrapper objects. 

* **CAUTION**: Due to the nature of the electrostatics formulation implemented in DFT-FE we strongly recommend to tile periodic cell lengths so that they are more than 10 atomic units


cpp API demo example
==========================================

* Install DFT-FE current development branch (publicGithubDevelop) from [DFT-FE github repo](https://github.com/dftfeDevelopers/dftfe). The installation instructions for DFT-FE and its dependencies are provied in the development version manual [here](https://github.com/dftfeDevelopers/dftfe/blob/manual/manual-develop.pdf). To use the dftfe library as an API we suggest performing a prefix based installation using CMake. For example, please refer to an example setup script [here](https://github.com/dftfeDevelopers/dftfe/blob/publicGithubDevelop/helpers/NERSCPerlmutterGPU/setupUserPerlmutterPrefixInstall.sh). Please note that two separate libraries will be installed, one for real datatype (Gamma point) and the other for the complex datatype (multiple k-points).

* Choose one of the demo examples by changing the source file in the CMakelists.txt. The description of the demo examples are given inside their respective .cpp files. The only header from the dftfe library which is required to be included in the .cpp file is `dftfeWrapper.h`.

* Create build directory and compile the example using CMake by passing the dftfe installation path for the complex datatype as that would work all cases. For example on NERSC Perlmutter machine we used:
```
$ module load PrgEnv-gnu
$ module load cudatoolkit
$ module unload cray-libsci/21.08.1.2
$ module load cmake/3.22.0
$ module load nccl/2.11.4
$ mkdir build 
$ cd build
$ cmake -DDFTFE_INSTALL_PATH=/global/common/software/m3916/softwareDFTFE/dftfe/installComplex -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc ../
$ make
```

* set `DFTFE_PSP_PATH` environment variable using export. The pseudpotential directory must contain ONCV files in the format: AtomicSymbol.upf

* Run the compiled `dftfeapi` executable as a MPI job
