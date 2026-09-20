#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib> // For std::system

#include "BreastCancer_3D.h"

using std::vector;
using std::string;

void BreastCancer_3D::SaveConfig_Adipocyte(string filename)
{
    FILE*   otptfl;
    if ((otptfl = fopen(filename.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s.\n", filename.c_str());
		exit(-1);
	}

	fprintf(otptfl, "%.32e\n", Lx);

	for (int n = 0; n < Na; n++)
		fprintf(otptfl, "%.32e\n", Kb_all[n]);
	for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", xa[n]);
    for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", ya[n]);
	for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", za[n]);

	fclose(otptfl);
}

void BreastCancer_3D::SaveConfig_Adipocyte_Comp(string filename)
{
    FILE*   otptfl;
    if ((otptfl = fopen(filename.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s.\n", filename.c_str());
		exit(-1);
	}

	fprintf(otptfl, "%.32e\n", Lx);

	for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", xa[n]);
    for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", ya[n]);
	for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", za[n]);

	fclose(otptfl);
}

void BreastCancer_3D::SaveConfig_Cancer(string filename)
{
    FILE*   otptfl;
    if ((otptfl = fopen(filename.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s.\n", filename.c_str());
		exit(-1);
	}

	fprintf(otptfl, "%.32e\n", Ly);

	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", Rc[n]);
	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", xc[n]);
    for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", yc[n]);
	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", zc[n]);

	fclose(otptfl);
}

void BreastCancer_3D::SaveConfig_All(string filename)
{
    FILE*   otptfl;
    if ((otptfl = fopen(filename.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s.\n", filename.c_str());
		exit(-1);
	}

	fprintf(otptfl, "%.32e\n", Lx);
	fprintf(otptfl, "%.32e\n", Ly);
	fprintf(otptfl, "%.32e\n", Lz);

	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", Rc[n]);

	for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", xa[n]);
    for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", ya[n]);
	for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", za[n]);
	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", xc[n]);
    for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", yc[n]);
	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", zc[n]);

	fclose(otptfl);
}

void BreastCancer_3D::SaveInv_Growth(string posname,
										string velname,
										vector<double> &Vx_a,
										vector<double> &Vy_a,
										vector<double> &Vz_a,
										vector<double> &Vx_c,
										vector<double> &Vy_c,
										vector<double> &Vz_c)
{
	FILE*	otptfl;
	int 	n;

	if ((otptfl = fopen(posname.c_str(), "a")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}

	fprintf(otptfl, "%.16e\n", Lx);
	fprintf(otptfl, "%.16e\n", Ly);
	fprintf(otptfl, "%.16e\n", Lz);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", xa[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", ya[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", za[n]);
	fprintf(otptfl, "%d\n", Nc);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Rc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", xc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", yc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", zc[n]);

	fclose(otptfl);

	if ((otptfl = fopen(velname.c_str(), "a")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}

	fprintf(otptfl, "%.16e\n", Vw);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", Vx_a[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", Vy_a[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", Vz_a[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vx_c[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vy_c[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vz_c[n]);

	fclose(otptfl);
}

void BreastCancer_3D::SaveInv_Growth_LipidLoss(string posname,
												string velname,
												vector<double> &rsc_all,
												vector<double> &Vx_a,
												vector<double> &Vy_a,
												vector<double> &Vz_a,
												vector<double> &Vx_c,
												vector<double> &Vy_c,
												vector<double> &Vz_c)
{
	FILE*	otptfl;
	int 	n;

	if ((otptfl = fopen(posname.c_str(), "a")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}

	fprintf(otptfl, "%.16e\n", Lx);
	fprintf(otptfl, "%.16e\n", Ly);
	fprintf(otptfl, "%.16e\n", Lz);
	for (n = 0; n < Na; n++)
		fprintf(otptfl, "%.16e\n", rsc_all[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", xa[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", ya[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", za[n]);
	fprintf(otptfl, "%d\n", Nc);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Rc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", xc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", yc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", zc[n]);

	fclose(otptfl);

	if ((otptfl = fopen(velname.c_str(), "a")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}

	fprintf(otptfl, "%.16e\n", Vw);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", Vx_a[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", Vy_a[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", Vz_a[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vx_c[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vy_c[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vz_c[n]);

	fclose(otptfl);
}

void BreastCancer_3D::SaveInv_Growth_OD(string posname)
{
	FILE*	otptfl;
	int 	n;

	if ((otptfl = fopen(posname.c_str(), "a")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}

	fprintf(otptfl, "%.16e\n", Lx);
	fprintf(otptfl, "%.16e\n", Ly);
	fprintf(otptfl, "%.16e\n", Lz);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", xa[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", ya[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", za[n]);
	fprintf(otptfl, "%d\n", Nc);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Rc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", xc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", yc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", zc[n]);

	fclose(otptfl);
}

void BreastCancer_3D::SaveInv_Growth_LipidLoss_OD(string posname,
													vector<double> &rsc_all)
{
	FILE*	otptfl;
	int 	n;

	if ((otptfl = fopen(posname.c_str(), "a")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}

	fprintf(otptfl, "%.16e\n", Lx);
	fprintf(otptfl, "%.16e\n", Ly);
	fprintf(otptfl, "%.16e\n", Lz);
	for (n = 0; n < Na; n++)
		fprintf(otptfl, "%.16e\n", rsc_all[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", xa[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", ya[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", za[n]);
	fprintf(otptfl, "%d\n", Nc);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Rc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", xc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", yc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", zc[n]);

	fclose(otptfl);
}

void BreastCancer_3D::LoadUnitParameter(double &D0_unit,
										double &Alf_min,
										vector<double> &A0_unit,
										double &A_sum_unit)
{
	std::ifstream pfile;
	pfile.open(dva_name.c_str());
    if (pfile.good()){
    	pfile >> D0_unit;
        for(int i = 0; i < Nf; i++)
            pfile >> A0_unit[i];
		pfile >> A_sum_unit;
		double V0_unit;
		pfile >> V0_unit;
		pfile >> Alf_min;
    }
	else{
		pfile.close();
		std::cout << "Couldn't find dva: " << dva_name << std::endl;
		exit(-1);
	}
	pfile.close();

	pfile.open(f_unit_name.c_str());
    if (pfile.good()){
        for (int i = 0; i < Nf; i++){
        	for (int j = 0; j < 3; j++){
        		pfile >> f_unit[i][j];
        	}
        }
    }
	else{
		pfile.close();
		std::cout << "Couldn't find f_unit: " << f_unit_name << std::endl;
		exit(-1);
	}
	pfile.close();

	pfile.open(edgelist_name.c_str());
    if (pfile.good()){
        for (int i = 0; i < Ne; i++){
        	for (int j = 0; j < 6; j++){
        		pfile >> edgelist[i][j];
        	}
        }
    }
	else{
		pfile.close();
		std::cout << "Couldn't find edgelist: " << edgelist_name << std::endl;
		exit(-1);
	}
	pfile.close();

	pfile.open(nnn_list_name.c_str());
    if (pfile.good()){
        for (int i = 0; i < Nnn; i++){
        	for (int j = 0; j < 2; j++){
        		pfile >> nnlist[i][j];
        	}
        }
    }
	else{
		pfile.close();
		std::cout << "Couldn't find nn_list: " << nnn_list_name << std::endl;
		exit(-1);
	}
	pfile.close();

	pfile.open(theta_name.c_str());
    if (pfile.good()){
        for (int i = 0; i < Ne; i++){
			pfile >> theta0[i];
        }
    }
	else{
		pfile.close();
		std::cout << "Couldn't find theta: " << theta_name << std::endl;
		exit(-1);
	}
	pfile.close();
}

void BreastCancer_3D::LoadUnitPosition(vector<double> &xyz_unit)
{
	std::ifstream pfile(unit_pos_name.c_str());
    if (pfile.good()){
        for(int i = 0; i < 3 * Ns; i++){
            pfile >> xyz_unit[i];
        }
    }
	else{
		pfile.close();
		std::cout << "Couldn't find " << unit_pos_name << std::endl;
		exit(-1);
	}
	pfile.close();
}

void BreastCancer_3D::Load_BumpyToDPM()
{
	char 			Unit_file[200], rigiddir[200], rigidfile[200];
	int 			i, ns, na;
	double 			D_bump_unit, D_bump, D_max_unit, D_max;
	vector<double>	X_body(Ns, 0.0), Y_body(Ns, 0.0), Z_body(Ns, 0.0);

	int 			N_2 = 2 * Na;
	int 			N_3 = 3 * Na;
	int 			N_4 = 4 * Na;

	std::ifstream file1(rigid_unit_name.c_str());
    if (file1.good()){
		file1 >> D_bump_unit;
		file1 >> D_max_unit;
		for (ns = 0; ns < Ns; ns++)
			file1 >> X_body[ns];
		for (ns = 0; ns < Ns; ns++)
			file1 >> Y_body[ns];
		for (ns = 0; ns < Ns; ns++)
			file1 >> Z_body[ns];
		file1.close();
    }
	else {
		file1.close();
        exit(-1);
	}
	
	double 			Rc, L;
	vector<double> 	Ccen(N_3, 0.0), Q(N_4, 0.0);
	std::ifstream bumpyfile(rigid_name.c_str());
	if (bumpyfile.good()){
		bumpyfile >> Rc;
		bumpyfile >> Lx;
		Ly = Lx;
		Lz = Lx;
		for (i = 0; i < N_3; i++)
			bumpyfile >> Ccen[i];
		for (i = 0; i < N_4; i++)
			bumpyfile >> Q[i];
		bumpyfile.close();
	}
	else {
		bumpyfile.close();
		exit(-1);
	}
	for (ns = 0; ns < Ns; ns++){
		X_body[ns] *= Rc;
		Y_body[ns] *= Rc;
		Z_body[ns] *= Rc;
	}

	// Bumpy center and orientation to bump center
	int 		idx1, idx;
	double 		xcn, ycn, zcn;
	double 		q_norm;
	double 		q0, q1, q2, q3;
	double 		q00, q11, q22, q33;
	double		q01, q02, q03;
	double 		q12, q13, q23;
	double 		R_mat_11, R_mat_12, R_mat_13;
	double 		R_mat_21, R_mat_22, R_mat_23;
	double 		R_mat_31, R_mat_32, R_mat_33;
	for (na = 0; na < Na; na++){
		xcn = Ccen[na];
		ycn = Ccen[na + Na];
		zcn = Ccen[na + N_2];
		q0 = Q[na];
		q1 = Q[na + Na];
		q2 = Q[na + N_2];
		q3 = Q[na + N_3];

		q00 = q0 * q0;
		q11 = q1 * q1;
		q22 = q2 * q2;
		q33 = q3 * q3;
		q01 = q0 * q1;
		q02 = q0 * q2;
		q03 = q0 * q3;
		q12 = q1 * q2;
		q13 = q1 * q3;
		q23 = q2 * q3;

		R_mat_11 = q00 + q11 - q22 - q33;
		R_mat_12 = 2.0 * (q12 - q03);
		R_mat_13 = 2.0 * (q13 + q02);
		R_mat_21 = 2.0 * (q12 + q03);
		R_mat_22 = q00 - q11 + q22 - q33;
		R_mat_23 = 2.0 * (q23 - q01);
		R_mat_31 = 2.0 * (q13 - q02);
		R_mat_32 = 2.0 * (q23 + q01);
		R_mat_33 = q00 - q11 - q22 + q33;

		for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
			idx = ns - idx_list[na];
			xa[ns] = xcn + R_mat_11 * X_body[idx] + R_mat_12 * Y_body[idx] + R_mat_13 * Z_body[idx];
			ya[ns] = ycn + R_mat_21 * X_body[idx] + R_mat_22 * Y_body[idx] + R_mat_23 * Z_body[idx];
			za[ns] = zcn + R_mat_31 * X_body[idx] + R_mat_32 * Y_body[idx] + R_mat_33 * Z_body[idx];
		}
	}
}

void BreastCancer_3D::LoadConfig_Cancer(string filename)
{
    std::ifstream cpos(filename.c_str());
	if (cpos.good()){
		cpos >> Ly;
        for (int i = 0; i < Nc; i++){
            cpos >> Rc[i];
			Dc[i] = 2.0 * Rc[i];
			Rc_sq[i] = Rc[i] * Rc[i];
			Dc_cube[i] = 8.0 * Rc_sq[i] * Rc[i];
		}
        for (int i = 0; i < Nc; i++)
            cpos >> xc[i];
        for (int i = 0; i < Nc; i++)
            cpos >> yc[i];
		for (int i = 0; i < Nc; i++)
            cpos >> zc[i];
		cpos.close();
	}
	else {
		cpos.close();
		printf("Cancer File Not Found!\n");
		exit(-1);
	}
}

void BreastCancer_3D::LoadConfig_Adipocyte(string filename, bool full_load)
{
    std::ifstream apos(filename.c_str());
	if (apos.good()){
		apos >> Lx;
		Ly = Lx;
		Lz = Lx;
        if (full_load){
            for (int i = 0; i < Na; i++)
                apos >> Kb_all[i];
            for (int i = 0; i < Ns_tot; i++)
                apos >> xa[i];
            for (int i = 0; i < Ns_tot; i++)
                apos >> ya[i];
			for (int i = 0; i < Ns_tot; i++)
                apos >> za[i];
        }
		apos.close();
	}
	else {
		apos.close();
		printf("Adipocyte File Not Found!\n");
		exit(-1);
	}
}

void BreastCancer_3D::LoadConfig_All(string filename)
{
    std::ifstream Allpos(filename.c_str());
    if (Allpos.good()){
		Allpos >> Lx;
		Allpos >> Ly;
		Allpos >> Lz;

		for (int i = 0; i < Nc; i++){
			Allpos >> Rc[i];
			Dc[i] = 2.0 * Rc[i];
			Rc_sq[i] = Rc[i] * Rc[i];
			Dc_cube[i] = 8.0 * Rc_sq[i] * Rc[i];
		}
		for (int i = 0; i < Ns_tot; i++)
			Allpos >> xa[i];
		for (int i = 0; i < Ns_tot; i++)
			Allpos >> ya[i];
		for (int i = 0; i < Ns_tot; i++)
			Allpos >> za[i];
		for (int i = 0; i < Nc; i++)
			Allpos >> xc[i];
		for (int i = 0; i < Nc; i++)
			Allpos >> yc[i];
		for (int i = 0; i < Nc; i++)
			Allpos >> zc[i];

		Allpos.close();

		double  xa_cen_temp, ya_cen_temp, za_cen_temp;
		for (int na = 0; na < Na; na++){
			xa_cen_temp = 0.0;
			ya_cen_temp = 0.0;
			za_cen_temp = 0.0;
			for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				xa_cen_temp += xa[ns];
				ya_cen_temp += ya[ns];
				za_cen_temp += za[ns];
				xa_pin[ns] = xa[ns];
				ya_pin[ns] = ya[ns];
				za_pin[ns] = za[ns];
			}
			xa_cen_pin[na] = xa_cen_temp / Ns;
			ya_cen_pin[na] = ya_cen_temp / Ns;
			za_cen_pin[na] = za_cen_temp / Ns;
		}
	}
	else{
		printf("Mix File Not Found: %s!\n", filename.c_str());
		Allpos.close();
		exit(-1);
	}
}