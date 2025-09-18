/*
***********************************************************************

Gme_fileread.cpp
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

#include "../include/Gme_fileread.h"

void GmestartupRead(string filename, int& fwd_only, int& GM, string& Fwd_Setting_File, string& Inverse_Setting_File,
	string& Mesh_Setting_File, string& Observed_data_File)
{
	int tt;
	ifstream fin(filename);
	if (!fin)
	{
		cerr << "Starup file read failed!" << endl;
		exit(0);
	}
	string input_line, tmp1, tmp2, tmp3;
	fin >> tmp1 >> input_line;
	for (char& c : input_line) {
		c = std::toupper(c);
	}
	if (input_line == "G" || input_line == "GRAV" || input_line == "GRA" || input_line == "GRAVITY")GM = 0;
	else
	{
		GM = 1;
	}
	/////
	fin >> tmp1 >> tmp2 >> tmp3 >> input_line;
	for (char& c : input_line) {
		c = std::toupper(c);
	}
	if (input_line == "FWD" || input_line == "FORWARD" || input_line == "FWDONLY")fwd_only = 1;
	else
	{
		fwd_only = 2;
	}
	fin >> tmp1 >> tmp2 >> tmp3 >> Fwd_Setting_File;
	ifstream fin2(Fwd_Setting_File);
	if (!fin2)
	{
		cerr << "Error reading Fwd_Setting_File" << endl;
		exit(0);
	}
	fin2.close();
	fin >> tmp1 >> tmp2 >> tmp3 >> Inverse_Setting_File;
	ifstream fin3(Inverse_Setting_File);
	if (!fin3)
	{
		cerr << "Error reading Inverse_Setting_File" << endl;
		exit(0);
	}
	fin3.close();
	fin >> tmp1 >> tmp2 >> tmp3 >> Mesh_Setting_File;
	ifstream fin4(Mesh_Setting_File);
	if (!fin4)
	{
		cerr << "Error reading Mesh_Setting_File" << endl;
		exit(0);
	}
	fin4.close();
	fin >> tmp1 >> tmp2 >> tmp3 >> Observed_data_File;
	ifstream fin5(Observed_data_File);
	if (!fin5)
	{
		cerr << "Error reading Observed_data_File" << endl;
		exit(0);
	}
	fin5.close();
	fin.close();
}