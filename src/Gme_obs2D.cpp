/*
***********************************************************************

Gme_obs2D.cpp	(Constructing Gme obs)
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

#include "../include/Gme_obs2D.h"

void GmeObs2D::Read_obsData(string filename, GmeMesh2D model)
{
	string tmp;
	ifstream fin(filename);
	if (!fin)
	{
		cerr << "Error reding srs file!" << endl;
		exit(0);
	}
	while (!fin.eof())
	{
		double xx, yy, zz, zz2;
		fin >> xx >> yy >> zz >> zz2;
		if (fin.fail())break;
		obsx.push_back(xx); obsy.push_back(yy); data.push_back(zz); data_err.push_back(zz2);
	}
	Num_site = obsx.size();
	for (int is = 0; is < Num_site; is++)
	{
		obsx[is] *= 1000.; obsy[is] *= 1000.;
	}
}

void GmeObs2D::Output(string filename)
{
	ofstream fout(filename);
	if (!fout)
	{
		cerr << "Error writing data file" << endl;
		exit(0);
	}
	for (int is = 0; is < Num_site; is++)
	{
		fout << obsx[is] / 1000. << " " << obsy[is] / 1000. << " " << data[is] << " " << data_err[is] << endl;
	}
	fout.close();
}
