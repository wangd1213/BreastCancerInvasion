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

#include "BreastCancer_2D.h"

using std::vector;
using std::string;

void BreastCancer_2D::SaveConfig_Adipocyte(string filename)
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

	fclose(otptfl);
}

void BreastCancer_2D::SaveConfig_Adipocyte_Comp(string filename)
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

	fclose(otptfl);
}

void BreastCancer_2D::SaveConfig_Cancer(string filename)
{
    FILE*   otptfl;
    if ((otptfl = fopen(filename.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s.\n", filename.c_str());
		exit(-1);
	}

	fprintf(otptfl, "%.32e\n", Lx);

	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", Rc[n]);
	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", xc[n]);
    for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", yc[n]);

	fclose(otptfl);
}

void BreastCancer_2D::SaveConfig_All(string filename)
{
    FILE*   otptfl;
    if ((otptfl = fopen(filename.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s.\n", filename.c_str());
		exit(-1);
	}

	fprintf(otptfl, "%.32e\n", Lx);
	fprintf(otptfl, "%.32e\n", Ly);

	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", Rc[n]);

	for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", xa[n]);
    for (int n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.32e\n", ya[n]);
	for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", xc[n]);
    for (int n = 0; n < Nc; n++)
		fprintf(otptfl, "%.32e\n", yc[n]);

	fclose(otptfl);
}

void BreastCancer_2D::SaveInv_Growth(string posname,
										string velname,
										vector<double> &Vx_a,
										vector<double> &Vy_a,
										vector<double> &Vx_c,
										vector<double> &Vy_c)
{
	FILE*	otptfl;
	int 	n;

	if ((otptfl = fopen(posname.c_str(), "a")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}

	fprintf(otptfl, "%.16e\n", Lx);
	fprintf(otptfl, "%.16e\n", Ly);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", xa[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", ya[n]);
	fprintf(otptfl, "%d\n", Nc);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Rc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", xc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", yc[n]);

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
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vx_c[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vy_c[n]);

	fclose(otptfl);
}

void BreastCancer_2D::SaveInv_Growth_LipidLoss(string posname,
												string velname,
												vector<double> &rsc_all,
												vector<double> &Vx_a,
												vector<double> &Vy_a,
												vector<double> &Vx_c,
												vector<double> &Vy_c)
{
	FILE*	otptfl;
	int 	n;

	if ((otptfl = fopen(posname.c_str(), "a")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}

	fprintf(otptfl, "%.16e\n", Lx);
	fprintf(otptfl, "%.16e\n", Ly);
	for (n = 0; n < Na; n++)
		fprintf(otptfl, "%.16e\n", rsc_all[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", xa[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", ya[n]);
	fprintf(otptfl, "%d\n", Nc);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Rc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", xc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", yc[n]);

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
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vx_c[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Vy_c[n]);

	fclose(otptfl);
}

void BreastCancer_2D::SaveInv_Growth_OD(string posname)
{
	FILE*	otptfl;
	int 	n;

	if ((otptfl = fopen(posname.c_str(), "a")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}

	fprintf(otptfl, "%.16e\n", Lx);
	fprintf(otptfl, "%.16e\n", Ly);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", xa[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", ya[n]);
	fprintf(otptfl, "%d\n", Nc);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Rc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", xc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", yc[n]);

	fclose(otptfl);
}

void BreastCancer_2D::SaveInv_Growth_LipidLoss_OD(string posname,
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
	for (n = 0; n < Na; n++)
		fprintf(otptfl, "%.16e\n", rsc_all[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", xa[n]);
	for (n = 0; n < Ns_tot; n++)
		fprintf(otptfl, "%.16e\n", ya[n]);
	fprintf(otptfl, "%d\n", Nc);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", Rc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", xc[n]);
	for (n = 0; n < Nc; n++)
		fprintf(otptfl, "%.16e\n", yc[n]);

	fclose(otptfl);
}

void BreastCancer_2D::LoadConfig_Cancer(string filename)
{
    std::ifstream cpos(filename.c_str());
	if (cpos.good()){
		cpos >> Lx;
        for (int i = 0; i < Nc; i++){
            cpos >> Rc[i];
			Dc[i] = 2.0 * Rc[i];
			Rc_sq[i] = Rc[i] * Rc[i];
		}
        for (int i = 0; i < Nc; i++)
            cpos >> xc[i];
        for (int i = 0; i < Nc; i++)
            cpos >> yc[i];
		cpos.close();
	}
	else {
		cpos.close();
		printf("Cancer File Not Found!\n");
		exit(-1);
	}
}

void BreastCancer_2D::LoadConfig_Adipocyte(string filename, bool full_load)
{
    std::ifstream apos(filename.c_str());
	if (apos.good()){
		apos >> Lx;
		Ly = Lx;
        if (full_load){
            for (int i = 0; i < Na; i++)
                apos >> Kb_all[i];
            for (int i = 0; i < Ns_tot; i++)
                apos >> xa[i];
            for (int i = 0; i < Ns_tot; i++)
                apos >> ya[i];
        }
		apos.close();
	}
	else {
		apos.close();
		printf("Adipocyte File Not Found!\n");
		exit(-1);
	}
}

void BreastCancer_2D::LoadConfig_Adipocyte_Comp(string filename)
{
    std::ifstream apos(filename.c_str());
	if (apos.good()){
		apos >> Lx;
		Ly = Lx;
		for (int i = 0; i < Ns_tot; i++)
			apos >> xa[i];
		for (int i = 0; i < Ns_tot; i++)
			apos >> ya[i];
		apos.close();
	}
	else {
		apos.close();
		printf("Adipocyte File Not Found!\n");
		exit(-1);
	}
}

void BreastCancer_2D::LoadConfig_All(string filename)
{
    std::ifstream Allpos(filename.c_str());
    if (Allpos.good()){
		Allpos >> Lx;
		Allpos >> Ly;

		for (int i = 0; i < Nc; i++){
			Allpos >> Rc[i];
			Dc[i] = 2.0 * Rc[i];
			Rc_sq[i] = Rc[i] * Rc[i];
		}
		for (int i = 0; i < Ns_tot; i++)
			Allpos >> xa[i];
		for (int i = 0; i < Ns_tot; i++)
			Allpos >> ya[i];
		for (int i = 0; i < Nc; i++)
			Allpos >> xc[i];
		for (int i = 0; i < Nc; i++)
			Allpos >> yc[i];

		Allpos.close();

		double  xa_cen_temp, ya_cen_temp;
		for (int na = 0; na < Na; na++){
			xa_cen_temp = 0.0;
			ya_cen_temp = 0.0;
			for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				xa_cen_temp += xa[ns];
				ya_cen_temp += ya[ns];
				xa_pin[ns] = xa[ns];
				ya_pin[ns] = ya[ns];
			}
			xa_cen_pin[na] = xa_cen_temp / (double)Ns[na];
			ya_cen_pin[na] = ya_cen_temp / (double)Ns[na];
		}
	}
	else{
		printf("Mix File Not Found: %s!\n", filename.c_str());
		Allpos.close();
		exit(-1);
	}
}

int BreastCancer_2D::LoadInvFile(string posname,
                                    string velname,
                                    vector<double> &Vx_a,
                                    vector<double> &Vy_a,
                                    vector<double> &Vx_c,
                                    vector<double> &Vy_c,
                                    vector<double> &theta_c)
{
	return 1;
	/*
    std::ifstream pfile(posname.c_str());
    std::ifstream vfile(velname.c_str());

    if (pfile.good() && vfile.good()){
        FILE*   otptfl;
        char cmd[1000], posfile_copy[300], velfile_copy[300];
		sprintf(posfile_copy, "%sPos_Inv_All_copy.txt", dirpre);
		sprintf(velfile_copy, "%sVel_Inv_All_copy.txt", dirpre);
        sprintf(cmd, "cp %s %s", posfile, posfile_copy);
        std::system(cmd);
        sprintf(cmd, "cp %s %s", velfile, velfile_copy);
        std::system(cmd);
        
        int n;
		int Npos = 2 * (1 + Ns_tot) + 3 * Nc;
		int Nvel = 2 * (Ns_tot + Nc) + 1;
		double temp;
		vector<double> pos_temp;
		vector<double> vel_temp;
		std::ifstream pfile(posfile_copy);
		std::ifstream vfile(velfile_copy);
		while (pfile >> temp)
			pos_temp.push_back(temp);
		while (vfile >> temp)
			vel_temp.push_back(temp);
		int nt_pos = pos_temp.size() / Npos;
		int nt_vel = vel_temp.size() / Nvel;
		int nt_min;
		if ((pos_temp.size() % Npos != 0) || (vel_temp.size() % Nvel != 0) || (nt_pos != nt_vel)){
			nt_min = (nt_pos > nt_vel) ? nt_vel : nt_pos;
			sprintf(cmd, "rm %s", posfile);
			std::system(cmd);
			sprintf(cmd, "rm %s", velfile);
			std::system(cmd);

			otptfl = fopen(posfile, "a");
			for (n = 0; n < nt_min * Npos; n++)
				fprintf(otptfl, "%.16e\n", pos_temp[n]);
			fclose(otptfl);
			
			otptfl = fopen(velfile, "a");
			for (n = 0; n < nt_min * Nvel; n++)
				fprintf(otptfl, "%.16e\n", vel_temp[n]);
			fclose(otptfl);
		}
		else
			nt_min = nt_pos;

		int n_start = (nt_min - 1) * Npos;
		Lx = pos_temp[n_start];
		n_start += 2;
		for (n = 0; n < Ns_tot; n++){
			xa[n] = pos_temp[n_start + 2 * n];
			ya[n] = pos_temp[n_start + 2 * n + 1];
		}
		n_start += 2 * Ns_tot;
		for (n = 0; n < Nc; n++){
			xc[n] = pos_temp[n_start + 3 * n];
			yc[n] = pos_temp[n_start + 3 * n + 1];
			theta_c[n] = pos_temp[n_start + 3 * n + 2];
		}

		n_start = (nt_min - 1) * Nvel;
		Vw = vel_temp[n_start];
		n_start++;
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] = vel_temp[n_start + 2 * n];
			Vy_a[n] = vel_temp[n_start + 2 * n + 1];
		}
		n_start += 2 * Ns_tot;
		for (n = 0; n < Nc; n++){
			Vx_c[n] = vel_temp[n_start + 2 * n];
			Vy_c[n] = vel_temp[n_start + 2 * n + 1];
		}

		sprintf(cmd, "rm %s", posfile_copy);
        std::system(cmd);
        sprintf(cmd, "rm %s", velfile_copy);
        std::system(cmd);

        return nt_min * save_every_step + 1;
    }
    else
        return 1;
	*/
}

void BreastCancer_2D::SaveInvFile(string posname,
                                    string velname,
                                    vector<double> &Vx_a,
                                    vector<double> &Vy_a,
                                    vector<double> &Vx_c,
                                    vector<double> &Vy_c,
                                    vector<double> &theta_c)
{
	/*
    FILE*   otptfl;
    if ((otptfl = fopen(posname.c_str(), "a")) == NULL) {
        fprintf(stderr, "Unable to open %s\n", posname.c_str());
        exit(1);
    }

    int n;

    fprintf(otptfl, "%.16e\n", Lx);
    fprintf(otptfl, "%.16e\n", Ly);
    for (n = 0; n < Ns_tot; n++){
        fprintf(otptfl, "%.16e\n", xa[n]);
        fprintf(otptfl, "%.16e\n", ya[n]);
    }
    for (n = 0; n < Nc; n++){
        fprintf(otptfl, "%.16e\n", xc[n]);
        fprintf(otptfl, "%.16e\n", yc[n]);
        fprintf(otptfl, "%.16e\n", theta_c[n]);
    }

    fclose(otptfl);

    if ((otptfl = fopen(velname.c_str(), "a")) == NULL) {
        fprintf(stderr, "Unable to open %s\n", velname.c_str());
        exit(1);
    }

    fprintf(otptfl, "%.16e\n", Vw);
    for (n = 0; n < Ns_tot; n++){
        fprintf(otptfl, "%.16e\n", Vx_a[n]);
        fprintf(otptfl, "%.16e\n", Vy_a[n]);
    }
    for (n = 0; n < Nc; n++){
        fprintf(otptfl, "%.16e\n", Vx_c[n]);
        fprintf(otptfl, "%.16e\n", Vy_c[n]);
    }

    fclose(otptfl);
	*/
}
