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
//This demo example runs a ground-state DFT on a single material system with MPI:
//
//Periodic BCC AlNi 2x2x2 16 atom super cell using DFT-FE with ONCV pseudopotentials,
//GGA PBE exchange-correlation, 500 K Fermi-Dirac smearing temperature, and 2x2x2
//MP shifted k-point grid. Spin un-polarized DFT calculations are used.
//GPU capability can be toggled on or off depending on the architecture
//
//deformCell and updateAtomPositions functionality are also tested
//


int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  //Initialize global handles for dftfe's dependencies
  dftfe::dftfeWrapper::globalHandlesInitialize();

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
  dftfe::dftfeWrapper dftfeWrapped(MPI_COMM_WORLD, 
		                  true, //GPU toggle 
				  atomicPositionsCart,
				  atomicNumbers,
				  cell,
				  std::vector< bool >{true, true, true}, //pbc
				  std::vector< unsigned int >{2, 2, 2},//MP grid
				  std::vector< bool >{true, true, true});//MP grid shift

  //performs ground-state DFT calculation and computes ground-state free energy.
  //ionic forces (first boolean flag) and cell stress (second boolean flag)
  //computation are set to true
  double energy = dftfeWrapped.computeDFTFreeEnergy(true, true);

  if (world_rank==0)
     std::cout << "DFT free energy (in Ha) : " << energy << std::endl;

  std::vector<std::vector<double>> forces=dftfeWrapped.getForcesAtoms();
  if (world_rank==0)
  {
    std::cout<<"---------Forces (in Ha/Bohr)----------"<<std::endl;
    for (unsigned int i=0; i<forces.size();i++)
      std::cout<<"Atomid: "<<i<<", x: "<<forces[i][0]<<", y: "<<forces[i][1]<<", z: "<<forces[i][2]<<std::endl;
  }


  std::vector<std::vector<double>> stress=dftfeWrapped.getCellStress();
  if (world_rank==0)
  {       
    std::cout<<"---------Stress (in Ha/Bohr^3)----------"<<std::endl; 
    for (unsigned int i=0; i<3;i++)
     for (unsigned int j=0; j<3;j++)
	std::cout<<"Stress component["<<i<<","<<j<<"]: "<<stress[i][j]<<std::endl;
  }

  //
  //Deform cell by applying affine volumetric strain
  //
  
  dftfeWrapped.deformCell(std::vector<std::vector<double>>{{1.002,0.0,0.0},{0.0,1.002,0.0},{0.0,0.0,1.002}});

  energy = dftfeWrapped.computeDFTFreeEnergy(true, true);

  if (world_rank==0)
  {
     std::cout<<std::endl;
     std::cout << "DFT free energy (in Ha) after cell deform : " << energy << std::endl;

  }

  forces=dftfeWrapped.getForcesAtoms();
  if (world_rank==0)
  {
    std::cout<<"---------Forces (in Ha/Bohr) after cell deform----------"<<std::endl;
    for (unsigned int i=0; i<forces.size();i++)
      std::cout<<"Atomid: "<<i<<", x: "<<forces[i][0]<<", y: "<<forces[i][1]<<", z: "<<forces[i][2]<<std::endl;
  }


  stress=dftfeWrapped.getCellStress();
  if (world_rank==0)
  {
    std::cout<<"---------Stress (in Ha/Bohr^3) after cell deform----------"<<std::endl;
    for (unsigned int i=0; i<3;i++)
     for (unsigned int j=0; j<3;j++)
        std::cout<<"Stress component["<<i<<","<<j<<"]: "<<stress[i][j]<<std::endl;
  }

  //
  //Update x position of first atom by 0.05 Bohr
  //

  std::vector<std::vector<double>> displacementsCart={{0.05,0.000000000000,0.000000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.00000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000},
                                      {0.000000000,0.0000000000,0.0000000000}};

  dftfeWrapped.updateAtomPositions(displacementsCart);

  energy = dftfeWrapped.computeDFTFreeEnergy(true, true);

  if (world_rank==0)
  {
     std::cout<<std::endl;
     std::cout << "DFT free energy (in Ha) after atom position update : " << energy << std::endl;
  }

  forces=dftfeWrapped.getForcesAtoms();
  if (world_rank==0)
  {
    std::cout<<"---------Forces (in Ha/Bohr) after atom position update----------"<<std::endl;
    for (unsigned int i=0; i<forces.size();i++)
      std::cout<<"Atomid: "<<i<<", x: "<<forces[i][0]<<", y: "<<forces[i][1]<<", z: "<<forces[i][2]<<std::endl;
  }


  stress=dftfeWrapped.getCellStress();
  if (world_rank==0)
  {
    std::cout<<"---------Stress (in Ha/Bohr^3) after atom position update----------"<<std::endl;
    for (unsigned int i=0; i<3;i++)
     for (unsigned int j=0; j<3;j++)
        std::cout<<"Stress component["<<i<<","<<j<<"]: "<<stress[i][j]<<std::endl;
  }

  //deletes dftfe wrapper object
  dftfeWrapped.clear();

  //Finalize global handles for dftfe's dependencies
  dftfe::dftfeWrapper::globalHandlesFinalize();

  MPI_Finalize();
  return 0;
}
