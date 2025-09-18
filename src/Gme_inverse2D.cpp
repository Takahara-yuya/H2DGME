/*
***********************************************************************

Gme_inverse2D.cpp
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

#include "../include/Gme_inverse2D.h"

void GmeInv2D::Init_Module(string Filename)
{
	ifstream fin(Filename);
	if (!fin) {
		cerr << "Error Opening " << Filename << " !" << endl;
		exit(1);
	}
	///
	string tmp1, tmp2, tmp3;
	fin >> tmp1 >> tmp2 >> tmp3 >> OutputFile_Root;
	fin >> tmp1 >> tmp2 >> tmp3 >> OutputTmp_Root;
	fin >> tmp1 >> tmp2 >> tmp3 >> LogFile_Root;
	fin >> tmp1 >> tmp2 >> alpha;
	fin >> tmp1 >> tmp2 >> damp;
	fin >> tmp1 >> tmp2 >> dt;
	fin >> tmp1 >> tmp2 >> stabilizer;
	fin >> tmp1 >> tmp2 >> ee;
	fin >> tmp1 >> tmp2 >> v2h;
	fin >> tmp1 >> tmp2 >> target_misfit;
	fin >> tmp1 >> tmp2 >> Max_iteration;
	fin.close();
	cerr << "--------------------------------Inversion parameters---------------------------------" << endl;
	cerr << "Max iteration:" << Max_iteration << endl;
	cerr << "Starting regularization trad-off: " << alpha << endl;
	cerr << "Data mistfit reduction threshold: " << dt << "%" << endl;
	if (stabilizer == 0)cerr << "Stabilizer: Minimum model matrix" << endl;
	if (stabilizer == 1)cerr << "Stabilizer: First-order smooth matrix" << endl;
	if (stabilizer > 1)cerr << "Focusing parameter: " << ee << endl;
	cerr << "Regularization trad-off attenuation rate: " << damp << endl;
	cerr << "-------------------------------------------------------------------------------------" << endl;
}

void GmeInv2D::ConjugateGradient(GmeFwd2D& fwd, GmeMesh2D& model, GmeMesh2D& mapr,
	GmeObs2D& data_obs, GmeObs2D& data_cal)
{
	//
	_Data_Weighting_Matrix_Cons(data_obs);
	//Jacobi
	fwd.Gme_get_Jacobi(model, data_obs, model.GM);
	//
	int nit = 0;
	VectorXd Wm, Rn, bn, g, g_last, d;
	MatrixXd Wm1;
	Rn.resize(data_obs.Num_site);
	double Data_Misfit, Data_Misfit_last, Model_Misfit, k;
	Data_Misfit_last = 999999.; Model_Misfit = 0;
	GmeMesh2D model_A;
	model_A = model;
	Wm1 = (fwd.Jacobi_inv.transpose() * fwd.Jacobi_inv);
	Wm = Wm1.diagonal().array().sqrt();
	for (int i = 0; i < Wm.size(); i++)Wm(i) = pow(Wm(i) + 0.6, 0.8);
	//Aw
	for (int i = 0; i < model.Num_Cen; i++)model.Val1D(i) *= Wm(i);
	for (int j = 0; j < data_obs.Num_site; j++)
		for (int i = 0; i < model.Num_Cen; i++)
		{
			fwd.Jacobi_inv.coeffRef(j, i) /= Wm(i);
		}
	//Aw
	while (1)
	{
		nit++;
		if (nit > Max_iteration)break;
		//fwd
		fwd.Fwd(model_A, data_cal);
		_Misfit_Calculate(Rn, Data_Misfit, data_obs, data_cal);
		if ((Data_Misfit_last - Data_Misfit) < dt * 0.01 * Data_Misfit_last && nit > 1)alpha = alpha * damp;
		Data_Misfit_last = Data_Misfit;
		_sm_bn(Model_Misfit, bn, Wm, model, mapr);
		if (myid == 0) {
			_Print_Information(nit, Data_Misfit, Model_Misfit, model_A, data_cal, fwd);
		}
		g = fwd.Jacobi_inv.transpose() * Wd.transpose() * Wd * Rn + alpha * bn;
		//exit(0);
		//conjugate
		if (nit == 1)d = g;
		else
		{
			double tmp;
			tmp = g.transpose() * g;
			tmp = tmp / (g_last.transpose() * g_last);
			d = g + tmp * g_last;
		}
		g_last = g;
		//
		_steplength_calulate(k, Data_Misfit, g, d, Wm, model, fwd, data_obs, data_cal);
		//
		model.Val1D = model.Val1D - k * d;
		for (int i = 0; i < model.Num_Cen; i++)
		{
			if (model.Val1D(i) < 0)
				model.Val1D(i) = 0;
		}
		for (int i = 0; i < model.Num_Cen; i++)model_A.Val1D(i) = model.Val1D(i) / Wm(i);
		//
	}
	//
}

void GmeInv2D::_Misfit_Calculate(VectorXd& Rn, double& Data_Misfit, GmeObs2D data_obs, GmeObs2D data_cal)
{
	int exclude = 0;
	for (int i = 0; i < data_obs.Num_site; i++)
	{
		if (data_obs.data_err[i] > 9000.)exclude += 1;
		Rn(i) = -data_obs.data[i] + data_cal.data[i];
	}
	for (int i = 0; i < data_obs.Num_site; i++)
	{
		Data_Misfit += pow(Rn(i) * Wd.coeff(i, i), 2);
	}
	Data_Misfit = sqrt(Data_Misfit / (data_obs.Num_site - exclude));
}

void GmeInv2D::_Data_Weighting_Matrix_Cons(GmeObs2D data_obs)
{
	cerr << "Constructing Data Weighting Matrix...... " << endl;
	vector<Triplet<double> > triplets;
	Wd.resize(data_obs.Num_site, data_obs.Num_site);
	for (int is = 0; is < data_obs.Num_site; is++)
	{
		triplets.push_back(Triplet<double>(is, is, 1. / data_obs.data_err[is]));
	}
	Wd.setZero(); // 
	for (const auto& triplet : triplets) {
		Wd.insert(triplet.row(), triplet.col()) = triplet.value();
	}
	Wd.makeCompressed(); // 
}

void GmeInv2D::_sm_bn(double& Model_Misfit, VectorXd& bn, VectorXd Wm, GmeMesh2D model, GmeMesh2D mapr)
{
	vector<Triplet<double> > triplets, triplets2, triplets3;
	if (stabilizer == 0)
	{
		bn = model.Val1D - mapr.Val1D;
	}
	else if (stabilizer == 1)
	{
		SparseMatrix<double> Wex, Wey;
		//vertical smooth
		for (int i = 0; i < model.Num_xCen; i++)
		{
			for (int j = 0; j < model.Num_yCen - 1; j++)
			{
				int np = j + i * model.Num_yCen;
				if (model.fix[i][j] == 0 && model.fix[i][j + 1] == 0)
				{
					triplets.push_back(Triplet<double>(np, np, 1. * v2h));//EVENLY WEIGHT
					triplets.push_back(Triplet<double>(np, np + 1, -1. * v2h));//EVENLY WEIGHT
				}
			}
		}
		Wey.resize(model.Num_Cen, model.Num_Cen);
		Wey.setFromTriplets(triplets.begin(), triplets.end());
		Wey.makeCompressed();
		//horizontal smooth
		double v, h, ssm;
		for (int nl = 0; nl < model.Num_yCen; nl++)//L layer
		{
			v = abs(model.yNode[nl + 1] - model.yNode[nl]);
			for (int np = 0; np < model.Num_xCen - 1; np++)
			{
				int tmp = nl + np * model.Num_yCen;
				if (model.fix[np][nl] == 0 && model.fix[np + 1][nl] == 0)
				{
					triplets2.push_back(Triplet<double>(np, np, 1.));//EVENLY WEIGHT
					triplets2.push_back(Triplet<double>(np, np + 1, -1.));//EVENLY WEIGHT				
				}
			}
		}
		Wex.resize(model.Num_Cen, model.Num_Cen);
		Wex.setFromTriplets(triplets2.begin(), triplets2.end());
		Wex.makeCompressed();
		//calculate bn
		bn = (Wex.transpose() * Wex + Wey.transpose() * Wey) * (model.Val1D - mapr.Val1D);
		Model_Misfit = alpha * (model.Val1D - mapr.Val1D).transpose() * Wex.transpose() * Wex * (model.Val1D - mapr.Val1D);
		Model_Misfit = +alpha * (model.Val1D - mapr.Val1D).transpose() * Wey.transpose() * Wey * (model.Val1D - mapr.Val1D);
	}
	else
	{

	}
}

void GmeInv2D::_steplength_calulate(double& k, double Data_Misfit, VectorXd g, VectorXd d, VectorXd Wm, GmeMesh2D model, GmeFwd2D fwd,
	GmeObs2D data_obs, GmeObs2D data_cal)
{
	GmeMesh2D model_test;
	GmeMesh2D model_A;
	//
	model_test = model;
	model_A = model;
	//start
	int count = 0;
	VectorXd Wmln, WdFmln, Rn;
	Rn.resize(data_obs.Num_site);
	double Data_Misfit_test;
	Wmln.resize(g.size());
	for (int i = 0; i < d.size(); i++)Wmln(i) = Wm(i) * d(i);
	WdFmln = Wd * fwd.Jacobi_inv * d;
	//
	double tmp1, tmp2;
	tmp1 = d.transpose() * g;
	tmp2 = WdFmln.transpose() * WdFmln;
	tmp2 += alpha * Wmln.transpose() * Wmln;
	k = tmp1 / tmp2;
	//while (1)
	//{
	//	count++;
	//	if (count > 10)break;
	//	model_test.Val1D = model.Val1D - k * d;
	//	for (int i = 0; i < model.Num_Cen; i++)model_A.Val1D(i) = model_test.Val1D(i) / Wm(i);
	//	fwd.Fwd(model_A, data_cal);
	//	_Misfit_Calculate(Rn, Data_Misfit_test, data_obs, data_cal);
	//	if (Data_Misfit_test < Data_Misfit)break;
	//	else k = k * 0.5;
	//}
}

void GmeInv2D::_Print_Information(double nit, double Data_Misfit, double Model_Misfit, GmeMesh2D model, GmeObs2D data_cal, GmeFwd2D fwd)
{
	cerr << endl;
	cerr << "********************************INVERSION********************************" << endl;
	cerr << "NUMBER OF ITERATION: " << nit << endl;
	cerr << "Data Misfit: " << Data_Misfit << " " << "Model Misfit: " << Model_Misfit << endl;
	cerr << "Reguralization trad-off: " << alpha << " " << "Damping factor: " << damp << endl;
	cerr << "***************************************************************************" << endl;
	model.Outputmesh(OutputFile_Root + "Inv_model.dat");
	model.Outputmesh(OutputTmp_Root + "inv_itr" + to_string(int(nit)) + ".dat");
	data_cal.Output(fwd.caldata_outroot + "data_cal.dat");
	data_cal.Output(fwd.caldata_outroot_tmp + "data-cal-itr" + to_string(int(nit)) + ".dat");
}







