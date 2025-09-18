#pragma once
#ifndef GME_FWD2D_H_
#define GME_FWD2D_H_

#include "CommonHeaders.h"
#include "Gme_mesh2D.h"
#include "Gme_obs2D.h"
#include "tool.h"
#include "Filter.h"

class GmeFwd2D
{
public:
	GmeFwd2D() {
	}

	~GmeFwd2D() {
	}
	//mpi parameters//
	int np, myid;
	//

	string caldata_outroot, caldata_outroot_tmp;

	int expand;

	double height;//up

	MatrixXd Jacobi_fwd;
	MatrixXd Jacobi_inv;

	void Init_fwd(string filename, GmeMesh2D model);

	void Gme_get_Jacobi(GmeMesh2D model, GmeObs2D data_obs, int GM);

	void Fwd(GmeMesh2D model, GmeObs2D& data_cal);

protected:

};

#endif
