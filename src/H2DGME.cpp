/*
***********************************************************************

H2DGME.cpp	(Main program)
This file is the main program of H2DGME.

***********************************************************************

Dec 01, 2023
Copyright 2023

Zuwei Huang
hzw1498218560@tongji.edu.cn
School of Ocean and Earth Science, Tongji University
Integrated Geophysics Group


Chongjin Zhao
zcjadc@126.com
School of Ocean and Earth Science, Tongji University
Integrated Geophysics Group


version 2.0.0

* ***********************************************************************
*/
#include "../include/CommonHeaders.h"
#include "../include/Gme_inverse2D.h"
#include "../include/Gme_fwd2D.h"
#include "../include/Gme_mesh2D.h"
#include "../include/Gme_obs2D.h"
#include "../include/Gme_fileread.h"

int main(int argc, char** argv)
{
	//mpi initialize//
	int myid, np;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	//
	if (myid == 0)
	{
		// Copyright and license notice:
		std::cout << std::endl;
		std::cout << "***************************************************************************" << std::endl << std::endl;
		std::cout << "                   H2DGME                                 " << std::endl << std::endl;
		std::cout << "  Copyright 2024, Zuwei Huang         " << std::endl << std::endl;

		std::cout << "  H2DGME is free software: you can redistribute it and/or modify it under the " << std::endl <<
			"  terms of the GNU Lesser General Public License as published by the Free " << std::endl <<
			"  Software Foundation, version 3 of the License. " << std::endl << std::endl;

		std::cout << "  H2DGME is distributed in the hope that it will be useful, but WITHOUT ANY " << std::endl <<
			"  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS " << std::endl <<
			"  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for " << std::endl <<
			"  more details. " << std::endl << std::endl;
		std::cout << "***************************************************************************" << std::endl << std::endl;
	}
	// The user should specify the configuration file as input to H2DGME:
	//************************** Read Startup File *************************************
	string startup = "Startup_GME";
	int GM, fwd_only;
	string Fwd_Setting_File, Inverse_Setting_File, Mesh_Setting_File, Observed_data_File;
	GmeMesh2D model, mapr; GmeObs2D data_obs; GmeObs2D data_cal; GmeFwd2D fwd; GmeInv2D Inv;
	GmestartupRead(startup, fwd_only, GM, Fwd_Setting_File, Inverse_Setting_File, Mesh_Setting_File, Observed_data_File);
	Inv.myid = myid; fwd.myid = myid;
	Inv.np = np;	fwd.np = np;
	///MPI_INIT_IN_FWD_AND_INV
		///
	if (myid == 0)cerr << "Starting constructing velocity model!" << endl;
	model.Init_model(Mesh_Setting_File, GM);
	if (myid == 0)std::cerr << "Constructing velocity model data sucessfully!" << endl;
	if (myid == 0)std::cerr << endl;
	mapr = model;
	
	if (myid == 0)std::cerr << "Starting reading observed data!" << endl;
	data_obs.Read_obsData(Observed_data_File, model);
	data_cal = data_obs;
	if (myid == 0)std::cerr << "Reading observed data sucessfully!" << endl;
	if (myid == 0)std::cerr << endl;
	///////
	std::cerr << endl;
	std::cerr << "NX: " << model.Num_xCen << " " << "NY:" << model.Num_yCen << endl;
	std::cerr << "XMIN XMAX: " << model.xmin << " " << model.xmax << " " << "YMIN YMAX: " << model.ymin << " " << model.ymax << endl;
	std::cerr << "Total observation station: " << data_obs.Num_site << endl;
	std::cerr << endl;
	///////
	if (fwd_only == 1)
	{
		fwd.Init_fwd(Fwd_Setting_File, model);
		fwd.Gme_get_Jacobi(model, data_obs, GM);
		fwd.Fwd(model, data_cal);
		model.Outputmesh("Result/model.dat");
		if (GM == 0)data_cal.Output(fwd.caldata_outroot + "fwd_grav_data.dat");
		else data_cal.Output(fwd.caldata_outroot + "fwd_mag_data.dat");
	}
	else
	{
		fwd.Init_fwd(Fwd_Setting_File, model);
		Inv.Init_Module(Inverse_Setting_File);
		Inv.ConjugateGradient(fwd, model, mapr, data_obs, data_cal);
	}

	MPI_Finalize();
}
