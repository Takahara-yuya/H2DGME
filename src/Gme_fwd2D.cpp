/*
***********************************************************************

Gme_fwd2D.cpp	(Forward modeling)
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

#include "../include/Gme_fwd2D.h"

void GmeFwd2D::Init_fwd(string filename, GmeMesh2D model)
{
	string tmp, tmp2, tmp3;
	ifstream fin(filename);
	fin >> tmp >> tmp >> caldata_outroot;
	fin >> tmp >> tmp >> tmp >> caldata_outroot_tmp;
	fin >> tmp >> tmp >> expand;
	fin >> tmp >> tmp >> height;
	height *= 1000.;
	fin.close();
}

void GmeFwd2D::Gme_get_Jacobi(GmeMesh2D model, GmeObs2D data_obs, int GM)
{
	Jacobi_inv.resize(data_obs.Num_site, model.Num_Cen);
	Jacobi_fwd.resize(data_obs.Num_site, (model.Num_xCen + 2 * expand) * model.Num_yCen);
	double zk, dx, dy, rnk, g, mu0, pi;
	g = 6.672e-3;
	pi = 3.1415926535;
	mu0 = 4 * pi * 1e-7 * 1e9;
	for (int is = 0; is < data_obs.Num_site; is++)
	{
		for (int ix = 0; ix < model.Num_xCen; ix++)
		{
			for (int iy = 0; iy < model.Num_yCen; iy++)
			{
				int np = iy + ix * model.Num_yCen;
				zk = abs(data_obs.obsy[is] - model.yCen[iy]) + height;
				dx = model.dX[ix]; dy = model.dY[iy];
				rnk = sqrt(pow(data_obs.obsx[is] - model.xCen[ix], 2) + pow(zk, 2));
				if (GM == 0)	Jacobi_inv.coeffRef(is, np) = 2 * g * dx * dy * zk / (rnk * rnk);
				else Jacobi_inv.coeffRef(is, np) = mu0 / (2 * pi) * (2.0 * zk * zk / (rnk * rnk) - 1.0) / (rnk * rnk) * dx * dy;
				if (iy == 0)Jacobi_inv.coeffRef(is, np) = 0;
			}
		}
		for (int ix = 0; ix < model.Num_xCen + 2 * expand; ix++)
		{
			for (int iy = 0; iy < model.Num_yCen; iy++)
			{
				int np = iy + ix * model.Num_yCen;
				zk = abs(data_obs.obsy[is] - model.yCen[iy]) + height;
				if (ix < expand)
				{
					dx = model.dX[0];
					rnk = sqrt(pow(data_obs.obsx[is] - (model.xCen[0] - dx * (expand - ix)), 2)
						+ pow(zk, 2));
				}
				else if (ix >= model.Num_xCen + expand)
				{
					dx = model.dX[model.Num_xCen - 1];
					rnk = sqrt(pow(data_obs.obsx[is] - (model.xCen[model.Num_xCen - 1] + dx * (ix - model.Num_xCen - expand + 1.)), 2)
						+ pow(zk, 2));
				}
				else
				{
					dx = model.dX[ix - expand];
					rnk = sqrt(pow(data_obs.obsx[is] - model.xCen[ix - expand], 2) + pow(zk, 2));
				}
				dy = model.dY[iy];
				if (GM == 0)	Jacobi_fwd.coeffRef(is, np) = 2 * g * dx * dy * zk / (rnk * rnk);
				else Jacobi_fwd.coeffRef(is, np) = mu0 / (2 * pi) * (2.0 * zk * zk / (rnk * rnk) - 1.0) / (rnk * rnk) * dx * dy;
			}
		}
	}
}

void GmeFwd2D::Fwd(GmeMesh2D model, GmeObs2D& data_cal)
{
	//expand
	VectorXd m;
	int np2;
	m.resize((model.Num_xCen + 2 * expand) * model.Num_yCen);
	for (int ix = 0; ix < model.Num_xCen + 2 * expand; ix++)
	{
		for (int iy = 0; iy < model.Num_yCen; iy++)
		{
			int np = iy + ix * model.Num_yCen;
			if (ix < expand)
			{
				np2 = iy;
				m(np) = model.Val1D(np2);
			}
			else if (ix >= model.Num_xCen + expand)
			{
				np2 = iy + (model.Num_xCen - 1) * model.Num_yCen;
				m(np) = model.Val1D(np2);
			}
			else
			{
				np2 = iy + (ix - expand) * model.Num_yCen;
				m(np) = model.Val1D(np2);
			}
		}
	}

	VectorXd data;
	data = Jacobi_fwd * m;
	for (int is = 0; is < data_cal.Num_site; is++)
		data_cal.data[is] = data(is);
}






