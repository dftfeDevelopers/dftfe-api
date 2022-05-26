// ---------------------------------------------------------------------
//
// Copyright (c) 2017-2022 The Regents of the University of Michigan and DFT-FE
// authors.
//
// This file is part of the DFT-FE code.
//
// The DFT-FE code is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the DFT-FE distribution.
//
// ---------------------------------------------------------------------
//
// @author Sambit Das
//

//
// dftfe header
//
#include <dftfeWrapper.h>

//mpi header
#include <mpi.h>

//
// C++ headers
//
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <sys/stat.h>

//
//This demo example runs two material systems with MPI on two separate MPI
//communicators created with even and odd MPI ranks respectively:
//
//Material system 1 on even MPI rank communicator:
//Periodic BCC AlNi 2x2x2 16 atom super cell using DFT-FE with ONCV pseudopotentials,
//GGA PBE exchange-correlation, 500 K Fermi-Dirac smearing temperature, and 2x2x2
//MP shifted k-point grid. Spin un-polarized DFT calculations are used.
//GPU capability can be toggled on or off depending on the architecture
//
//Material system 2 on odd MPI rank communicator:
//Periodic BCC Fe 2x2x2 16 atom super cell using DFT-FE with ONCV pseudopotentials,
//GGA PBE exchange-correlation, 500 K Fermi-Dirac smearing temperature, and 2x2x2
//MP shifted k-point grid. Spin polarized DFT calculations are used.
//GPU capability can be toggled on or off depending on the architecture



int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  // Check that even number of MPI processes are used
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if(comm_size%2!= 0)
  {
      printf("This application is meant to be run with even number of processes, not %d.\n", comm_size);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // Get my rank in the global communicator
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // Determine the colour and key based on whether my rank is even.
  int colour;
  int key;
  if(my_rank % 2 == 0)
  {
      colour = 0;
      key = my_rank;
  }
  else
  {
      colour = 1;
      key = comm_size - my_rank;
  }

  // Split the global communicator
  MPI_Comm new_comm;
  MPI_Comm_split(MPI_COMM_WORLD, colour, key, &new_comm);

  // Get my rank in the new communicator
  int my_new_comm_rank;
  MPI_Comm_rank(new_comm, &my_new_comm_rank);


  dftfe::dftfeWrapper::globalHandlesInitialize();


  //Material system 1
  if (new_comm!=MPI_COMM_NULL && colour==0)
  {
    //in Bohr units
    std::vector<std::vector<double>> cell(3, std::vector<double>(3, 0.0));
    cell[0][0] = 10.9146877116;
    cell[1][1] = 10.9146877116;
    cell[2][2] = 10.9146877116;

    //in Bohr units
    std::vector<std::vector<double>> atomicPositionsCart={{0.000000000000,0.000000000000,0.000000000000},
                                      {0.250000000000,0.250000000000,0.250000000000},
                                      {0.000000000000,0.000000000000,0.500000000000},
                                      {0.250000000000,0.250000000000,0.750000000000},
                                      {0.000000000000,0.500000000000,0.000000000000},
                                      {0.250000000000,0.750000000000,0.250000000000},
                                      {0.000000000000,0.500000000000,0.500000000000},
                                      {0.250000000000,0.750000000000,0.750000000000},
                                      {0.500000000000,0.000000000000,0.000000000000},
                                      {0.750000000000,0.250000000000,0.250000000000},
                                      {0.500000000000,0.000000000000,0.500000000000},
                                      {0.750000000000,0.250000000000,0.750000000000},
                                      {0.500000000000,0.500000000000,0.000000000000},
                                      {0.750000000000,0.750000000000,0.250000000000},
                                      {0.500000000000,0.500000000000,0.500000000000},
                                      {0.750000000000,0.750000000000,0.750000000000}};

    for (unsigned int i=0;i<atomicPositionsCart.size();++i)
      for (unsigned int j=0;j<3;++j)
        atomicPositionsCart[i][j]*=10.9146877116;

    std::vector<unsigned int> atomicNumbers={13,28,13,28,13,28,13,28,13,28,13,28,13,28,13,28};

    //constructs dftfe wrapper object
    dftfe::dftfeWrapper dftfeWrappedObject(new_comm,
                                  true, //GPU toggle 
                                  atomicPositionsCart,
                                  atomicNumbers,
                                  cell,
                                  std::vector< bool >{true, true, true}, //pbc
                                  std::vector< unsigned int >{2, 2, 2},//MP grid
                                  std::vector< bool >{true, true, true});//MP grid shift

    //performs ground-state DFT calculation and computes ground-state free energy.
    //ionic forces (first boolean flag) and cell stress (second boolean flag)
    //computation are set to false
    const double energy = dftfeWrappedObject.computeDFTFreeEnergy(false, false);

    if (my_new_comm_rank==0)
      std::cout << "DFT free energy for system 1: " << energy << std::endl;

    //clear call is not required as the object is automatically destroyed when
    //it goes out of scope
    //dftfeWrappedObject.clear();
  }

  //Material system 2
  if (new_comm!=MPI_COMM_NULL && colour==1)
  {
    //in Bohr units
    std::vector<std::vector<double>> cell(3, std::vector<double>(3, 0.0));
    cell[0][0] = 10.73572298;
    cell[1][1] = 10.73572298;
    cell[2][2] = 10.73572298;

    //in Bohr units
    std::vector<std::vector<double>> atomicPositionsCart={{0.000000000000,0.000000000000,0.000000000000},
                                      {0.250000000000,0.250000000000,0.250000000000},
                                      {0.000000000000,0.000000000000,0.500000000000},
                                      {0.250000000000,0.250000000000,0.750000000000},
                                      {0.000000000000,0.500000000000,0.000000000000},
                                      {0.250000000000,0.750000000000,0.250000000000},
                                      {0.000000000000,0.500000000000,0.500000000000},
                                      {0.250000000000,0.750000000000,0.750000000000},
                                      {0.500000000000,0.000000000000,0.000000000000},
                                      {0.750000000000,0.250000000000,0.250000000000},
                                      {0.500000000000,0.000000000000,0.500000000000},
                                      {0.750000000000,0.250000000000,0.750000000000},
                                      {0.500000000000,0.500000000000,0.000000000000},
                                      {0.750000000000,0.750000000000,0.250000000000},
                                      {0.500000000000,0.500000000000,0.500000000000},
                                      {0.750000000000,0.750000000000,0.750000000000}};

    for (unsigned int i=0;i<atomicPositionsCart.size();++i)
      for (unsigned int j=0;j<3;++j)
        atomicPositionsCart[i][j]*=10.73572298;

    std::vector<unsigned int> atomicNumbers={26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26};

    //constructs dftfe wrapper object
    dftfe::dftfeWrapper dftfeWrappedObject(new_comm,
                                  true, //GPU toggle
                                  atomicPositionsCart,
                                  atomicNumbers,
                                  cell,
                                  std::vector< bool >{true, true, true}, //pbc
                                  std::vector< unsigned int >{2, 2, 2},//MP grid
                                  std::vector< bool >{true, true, true},//MP grid shift
                                  true,//spin polarization toggle
                                  0.1);//starting magnetization

    //performs ground-state DFT calculation and computes ground-state free energy.
    //ionic forces (first boolean flag) and cell stress (second boolean flag)
    //computation are set to false
    const double energy = dftfeWrappedObject.computeDFTFreeEnergy(false, false);


    if (my_new_comm_rank==0)
      std::cout << "DFT free energy for system 2: " << energy << std::endl;

    //clear call not required as the object is automatically destroyed when it
    //goes out of scope
    //dftfeWrappedObject.clear();
  }



  dftfe::dftfeWrapper::globalHandlesFinalize();

  MPI_Comm_free(&new_comm);

  MPI_Finalize();
  return 0;
}
