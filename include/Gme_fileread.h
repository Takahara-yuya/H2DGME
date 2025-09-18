#pragma once
#ifndef TOMO_FILEREAD_H_
#define TOMO_FILEREAD_H_

#include "CommonHeaders.h"
#include "Gme_fwd2D.h"
#include "Gme_mesh2D.h"
#include "Gme_obs2D.h"

void GmestartupRead(string filename, int& fwd_only,int &GM, string& Fwd_Setting_File, string& Inverse_Setting_File,
	string& Mesh_Setting_File, string& Observed_data_File);

#endif