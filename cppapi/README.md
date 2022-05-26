DFT-FE's cpp API is already in-built into the library in the form of dftfeWrapper class. The doxygen documentation of this class is provided [here](https://dftfedevelopers.github.io/dftfe/classdftfe_1_1dftfe_wrapper.html)



Instructions to use cpp API
==========================================

* Install DFT-FE current development branch (publicGithubDevelop) from [DFT-FE github repo](https://github.com/dftfeDevelopers/dftfe). The installation instructions for DFT-FE and its dependencies are provied in the development version manual [here](https://github.com/dftfeDevelopers/dftfe/blob/manual/manual-develop.pdf). To use the dftfe library as an API we suggest performing a prefix based installation using CMake. For example, please refer to an example setup script [here](https://github.com/dftfeDevelopers/dftfe/blob/publicGithubDevelop/helpers/NERSCPerlmutterGPU/setupUserPerlmutterPrefixInstall.sh).

* set DFTFE_PSP_PATH environment variable

* cmake -DDFTFE_INSTALL_PATH=/global/common/software/m3916/softwareDFTFE/dftfe/installReal -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc  
