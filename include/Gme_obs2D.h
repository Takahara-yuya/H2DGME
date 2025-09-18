#pragma once
#ifndef GME_OBS2D_H_
#define GME_OBS2D_H_
#include "CommonHeaders.h"
#include "Gme_mesh2D.h"

class GmeObs2D
{
public:

	GmeObs2D() {
	}

	~GmeObs2D() {
	}

	//Basic Mesh
	int Num_site;//site number
	vector<double> obsx, obsy, data, data_err;

	void Read_obsData(string filename, GmeMesh2D model);

	void Output(string filename);

protected:

};

#endif
