#pragma once
#ifndef GME_MESH2D_H_
#define GME_MESH2D_H_

#include "CommonHeaders.h"
#include "Interpolation.h"
#include "tool.h"

class GmeMesh2D
{
public:
	//Basic Mesh
	GmeMesh2D() {
	}

	~GmeMesh2D() {
	}
	bool Marine;

	int GM;

	int Num_xCen, Num_yCen, Num_Cen;
	int Num_xNode, Num_yNode, Num_Node;

	double xmin, xmax, ymin, ymax;

	vector<double>  xCen, yCen, xNode, yNode, xTopo, dX, dY, upbound;
	vector<vector<double>>  fix;

	VectorXd Val1D;

	void Init_model(string filename, int GM1);

	void Outputmesh(string filename);

protected:

	double den_start, d_den;

	void _Make1Dmodel(int type, int GM);

	void _Read_Gridfile(string filename);

	void _Read_Initfile(string filename);

	void _Read_Topofile(string filename);
};

#endif