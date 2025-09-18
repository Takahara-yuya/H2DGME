/*
***********************************************************************

Gme_mesh2D.cpp	(Constructing Gme mesh)
This file is part of H2DGME.

***********************************************************************

Nov 22, 2023
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


***********************************************************************
*/

#include "../include/Gme_mesh2D.h"

#define den_water 1.05

void GmeMesh2D::Init_model(string filename, int GM1)
{
	string tmp, uniform, gridfile, initfile, topofile, model_type, data_type;
	double dx, dy;
	ifstream fin(filename);
	GM = GM1;
	if (!fin)
	{
		cerr << "Error reading meshsetting file!" << endl;
		exit(1);
	}
	fin >> tmp >> tmp >> data_type;
	fin >> tmp >> tmp >> Num_xCen >> Num_yCen;
	fin >> tmp >> tmp >> uniform;
	fin >> tmp >> tmp >> dx >> dy;
	fin >> tmp >> tmp >> gridfile;
	fin >> tmp >> tmp >> topofile;
	fin >> tmp >> tmp >> xmin >> ymin;
	fin >> tmp >> tmp >> model_type;
	fin >> tmp >> tmp >> initfile;
	fin >> tmp >> tmp >> den_start >> d_den;
	fin.close();
	//
	dx = dx * 1000; dy = dy * 1000;
	xmin = xmin * 1000; ymin = ymin * 1000.;
	//read finished! construct model
	//string change:
	for (char& c : model_type) {
		c = std::tolower(c);
	}
	for (char& c : uniform) {
		c = std::tolower(c);
	}
	for (char& c : data_type) {
		c = std::tolower(c);
	}
	//data_type
	if (data_type == "land")Marine = false;
	else Marine = true;
	//allocated parameter
	Num_xNode = Num_xCen + 1; Num_yNode = Num_yCen + 1;
	Num_Cen = Num_xCen * Num_yCen; Num_Node = Num_xNode * Num_yNode;
	dX.resize(Num_xCen); dY.resize(Num_yCen);
	xCen.resize(Num_xCen); yCen.resize(Num_yCen);
	xTopo.resize(Num_xCen);
	upbound.resize(Num_xCen);
	xNode.resize(Num_xNode); yNode.resize(Num_yNode);
	//
	Val1D.resize(Num_Cen);
	//add topo
	fix.resize(Num_xCen);
	for (int i = 0; i < Num_xCen; ++i) {
		fix[i].resize(Num_yCen);
	}
	for (int i = 0; i < Num_xCen; ++i)
	{
		for (int j = 0; j < Num_yCen; ++j)
			fix[i][j] = 0;
	}
	//decide dx,dy,xCen,yCen
	if (uniform == "yes")
	{
		for (int i = 0; i < Num_xCen; i++)dX[i] = dx;
		for (int i = 0; i < Num_yCen; i++)dY[i] = dx;
		for (int i = 0; i < Num_xNode; i++)xNode[i] = xmin + dx * i;
		for (int i = 0; i < Num_yNode; i++)yNode[i] = ymin + dy * i;
		for (int i = 0; i < Num_xCen; i++)xCen[i] = 0.5 * (xNode[i] + xNode[i + 1]);
		for (int i = 0; i < Num_yCen; i++)yCen[i] = 0.5 * (yNode[i] + yNode[i + 1]);
	}
	else
	{
		_Read_Gridfile(gridfile);
	}
	//model make
	if (model_type == "gradient1")
	{	//read topo
		_Read_Topofile(topofile);
		_Make1Dmodel(1, GM1);
	}
	else if (model_type == "gradient0")
	{
		_Read_Topofile(topofile);
		_Make1Dmodel(0, GM1);
	}
	else
	{
		_Read_Topofile(topofile);
		_Read_Initfile(initfile);
	}
}

void GmeMesh2D::Outputmesh(string filename)
{
	ofstream fout(filename);
	if (!fout)
	{
		cerr << "Error opening output file!" << endl;
		exit(1);
	}
	for (int ix = 0; ix < Num_xCen; ix++)
		for (int iy = 0; iy < Num_yCen; iy++)
		{
			fout << xCen[ix] / 1000. << " " << yCen[iy] / 1000. << " " << Val1D(iy + ix * Num_yCen) << endl;
		}
	fout.close();
}

void GmeMesh2D::_Make1Dmodel(int type, int GM)
{
	if (GM == 0)
	{
		for (int ix = 0; ix < Num_xCen; ix++)
		{
			for (int iy = 0; iy < upbound[ix]; iy++)
			{
				if (!Marine)Val1D(iy + ix * Num_yCen) = 0 * 1000;
				else Val1D(iy + ix * Num_yCen) = den_water * 1000;
			}
			for (int iy = upbound[ix]; iy < Num_yCen; iy++)
			{
				Val1D(iy + ix * Num_yCen) = den_start + d_den * (iy - upbound[ix] * double(type));
				Val1D(iy + ix * Num_yCen) = Val1D(iy + ix * Num_yCen) * 1000.;
			}
		}
	}
	else
	{
		Val1D.setZero();
	}
}

void GmeMesh2D::_Read_Gridfile(string filename)
{
	string tmp;
}

void GmeMesh2D::_Read_Initfile(string filename)
{
	double tmp1, tmp2;
	ifstream fin(filename);
	if (!fin)
	{
		cerr << "Error reading Init file!" << endl;
		exit(1);
	}
	for (int ix = 0; ix < Num_xCen; ix++)
		for (int iy = 0; iy < Num_yCen; iy++)
		{
			fin >> tmp1 >> tmp2 >> Val1D(iy + ix * Num_yCen);
		}
	fin.close();
}

void GmeMesh2D::_Read_Topofile(string filename)
{
	double tmp1, tmp2;
	vector<double> topo_tmp, x_tmp;
	ifstream fin(filename);
	if (!fin)
	{
		cerr << "Error reading Topo file!" << endl;
		exit(1);
	}
	while (!fin.eof())
	{
		fin >> tmp1 >> tmp2;
		//if (!Land && tmp2 < 0)
		//{
		////	cerr << "Marine topo should larger than zero!" << endl;
		////	exit(1);
		//	tmp2 = 0;
		//}
		x_tmp.push_back(tmp1 * 1000.); topo_tmp.push_back(tmp2 * 1000.);
	}
	fin.close();
	//interplote
	Interpolation inter;
	inter.linearinter1D_init(x_tmp, topo_tmp);
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		xTopo[ix] = inter.linearinter1D(xCen[ix]);
	}
	//adjust yMin
	auto minElement = std::min_element(xTopo.begin(), xTopo.end());
	tmp1 = *minElement;
	if (tmp1 < 0)
	{
		ymin = tmp1;
		for (int iy = 0; iy < Num_yNode; iy++)
		{
			yNode[iy] += ymin;
			if (iy < Num_yCen)
				yCen[iy] += ymin;
		}
	}
	else
		ymin = yNode[0];
	//
	xmin = xNode[0]; xmax = xNode[Num_xNode - 1];
	ymin = yNode[0]; ymax = yNode[Num_yNode - 1];
	//
	for (int i = 0; i < Num_xCen; i++)
	{
		for (int j = 0; j < Num_yCen; j++)
		{
			if (yCen[j] <= xTopo[i])
			{
				if (!Marine)
					fix[i][j] = 1;
				else
					fix[i][j] = 2;
			}
			else fix[i][j] = 0;
		}
		for (int j = 0; j < Num_yCen - 1; j++)
		{
			if (fix[i][0] == 0)
			{
				upbound[i] = 0;
				break;
			}
			if (fix[i][j] != 0 && fix[i][j + 1] == 0)
			{
				upbound[i] = j + 1;
				break;
			}
		}
	}
}

