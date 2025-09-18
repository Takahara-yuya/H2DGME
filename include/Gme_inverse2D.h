#pragma once
#ifndef GME_INVERSE2D_H_
#define GME_INVERSE2D_H_
#include "CommonHeaders.h"
#include "Gme_fwd2D.h"
#include "Gme_mesh2D.h"
#include "Gme_obs2D.h"

class GmeInv2D
{
public:
	//! Constructor
	GmeInv2D() {

	};
	//! Destructor
	~GmeInv2D() {

	};
	int myid, np;
	//
	double ee = 0.1;

	string OutputFile_Root, LogFile_Root, OutputTmp_Root;

	int stabilizer, Max_iteration, max_nls, pre_type, cg_type, cg_break;

	double target_misfit;

	double alpha, damp, dt, v2h, gamma, smdw, t1, t2, static_shift, step_gn;

	VectorXd Rn;//data misfit
	SparseMatrix<double> Wd, Wm, We, Wex, Wey;//data weighting ;//model weighting matrix
	//terminate condition
	int Bad_reduction;

	void Init_Module(string Filename);

	void ConjugateGradient(GmeFwd2D& fwd, GmeMesh2D& model, GmeMesh2D& mapr,
		GmeObs2D& data_obs, GmeObs2D& data_cal);

private:
	void _sm_bn(double& Model_Misfit, VectorXd& bn, VectorXd Wm, GmeMesh2D model, GmeMesh2D mapr);

	void _Misfit_Calculate(VectorXd& Rn, double& Data_Misfit, GmeObs2D data_obs, GmeObs2D data_cal);

	void _Data_Weighting_Matrix_Cons(GmeObs2D data_obs);

	void _steplength_calulate(double& k, double Data_Misfit, VectorXd g, VectorXd d, VectorXd Wm, GmeMesh2D model, GmeFwd2D fwd,
		GmeObs2D data_obs, GmeObs2D data_cal);

	void _Print_Information(double nit, double Data_Misfit, double Model_Misfit, GmeMesh2D model, GmeObs2D data_cal, GmeFwd2D fwd);

protected:



};
#endif
