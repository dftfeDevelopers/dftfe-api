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

//
// C++ headers
//
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <sys/stat.h>



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



  if (new_comm!=MPI_COMM_NULL && colour==0)
  {
    std::vector<std::vector<double>> cell1(3, std::vector<double>(3, 0.0));
    cell1[0][0] = 20.0;
    cell1[1][1] = 20.0;
    cell1[2][2] = 20.0;

    std::vector<std::vector<double>> atomicPositionsCart1(
      2, std::vector<double>(3, 0.0));
    atomicPositionsCart1[0][0] = 8.0;
    atomicPositionsCart1[0][1] = 10.0;
    atomicPositionsCart1[0][2] = 10.0;
    atomicPositionsCart1[1][0] = 11.0;
    atomicPositionsCart1[1][1] = 10.0;
    atomicPositionsCart1[1][2] = 10.0;

    std::vector<unsigned int> atomicNumbers1(atomicPositionsCart1.size());
    atomicNumbers1[0] = 8;
    atomicNumbers1[1] = 6;

    std::vector<bool> pbc1(3, false);

    dftfe::dftfeWrapper dftfeWrappedObject1(
      new_comm, true, atomicPositionsCart1, atomicNumbers1, cell1, pbc1);

    const double energy = dftfeWrappedObject1.computeDFTFreeEnergy(true, false);

    if (my_new_comm_rank==0)
      std::cout << "DFT free energy for system 1: " << energy << std::endl;

    dftfeWrappedObject1.clear();
  }

  if (new_comm!=MPI_COMM_NULL && colour==1)
  {
    std::vector<std::vector<double>> cell2(3, std::vector<double>(3, 0.0));
    cell2[0][0] = 20.0;
    cell2[1][1] = 20.0;
    cell2[2][2] = 20.0;

    std::vector<std::vector<double>> atomicPositionsCart2(
      2, std::vector<double>(3, 0.0));
    atomicPositionsCart2[0][0] = 8.2;
    atomicPositionsCart2[0][1] = 10.0;
    atomicPositionsCart2[0][2] = 10.0;
    atomicPositionsCart2[1][0] = 10.8;
    atomicPositionsCart2[1][1] = 10.0;
    atomicPositionsCart2[1][2] = 10.0;

    std::vector<unsigned int> atomicNumbers2(atomicPositionsCart2.size());
    atomicNumbers2[0] = 8;
    atomicNumbers2[1] = 6;

    std::vector<bool> pbc2(3, false);

    dftfe::dftfeWrapper dftfeWrappedObject2(
      new_comm, true, atomicPositionsCart2, atomicNumbers2, cell2, pbc2);

    const double energy = dftfeWrappedObject2.computeDFTFreeEnergy(true, false);

    if (my_new_comm_rank==0)
      std::cout << "DFT free energy for system 2: " << energy << std::endl;

    dftfeWrappedObject2.clear();
  }



  dftfe::dftfeWrapper::globalHandlesFinalize();

  MPI_Comm_free(&new_comm);

  MPI_Finalize();
  return 0;
}
