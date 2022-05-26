DFT-FE's cpp API is already in-built into the library in the form of dftfeWrapper class. The doxygen documentation of this class is provided [here](https://dftfedevelopers.github.io/dftfe/classdftfe_1_1dftfe_wrapper.html). The wrapper currently only sets up GGA PBE pseudopotential DFT calculations using ONCV pseudopotentials. All types of boundary conditions (non-periodic, periodic and semi-periodic) are supported by the wrapper.



cpp API demo example
==========================================

* Install DFT-FE current development branch (publicGithubDevelop) from [DFT-FE github repo](https://github.com/dftfeDevelopers/dftfe). The installation instructions for DFT-FE and its dependencies are provied in the development version manual [here](https://github.com/dftfeDevelopers/dftfe/blob/manual/manual-develop.pdf). To use the dftfe library as an API we suggest performing a prefix based installation using CMake. For example, please refer to an example setup script [here](https://github.com/dftfeDevelopers/dftfe/blob/publicGithubDevelop/helpers/NERSCPerlmutterGPU/setupUserPerlmutterPrefixInstall.sh). Please note that two separate libraries will be installed, one for real datatype (Gamma point) and the other for the complex datatype (multiple k-points).

* Choose one of the demo examples by changing the source file in the CMakelists.txt. The description of the demo examples are given inside their respective .cpp files. The only header from the dftfe library which is required to be included in the .cpp file is ``dftfeWrapper.h''.

* Create build directory and compile the example using Cmake by passing the dftfe installation path for the complex datatype as that would work all cases. For example on NERSC Perlmutter machine we used:
```
cmake -DDFTFE_INSTALL_PATH=/global/common/software/m3916/softwareDFTFE/dftfe/installComplex -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc
```

* set ``DFTFE_PSP_PATH'' environment variable using export. The pseudpotential directory must contain ONCV files in the format: AtomicSymbol.upf

* Run the compiled executable as a MPI job
