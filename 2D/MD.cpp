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

void BreastCancer_2D::SetDamping(double gamma_a_ratio,
                                    double gamma_c_ratio,
                                    double gamma_ac_ratio,
                                    double &gamma_a,
                                    double &gamma_c,
                                    double &gamma_ac)
{
    // damping coefficients
    double                  gamma_aa = Kaa / D0[Na - 1] / D0[Na - 1];
    double                  gamma_cc = Kcc / Rc[0] / Rc[0] / 4.0;
    gamma_ac = 4.0 * Kac / (D0[Na - 1] + 2.0 * Rc[0]) / (D0[Na - 1] + 2.0 * Rc[0]);
    gamma_a = gamma_a_ratio * std::sqrt(gamma_aa);
    gamma_c = gamma_c_ratio * std::sqrt(gamma_cc);
    gamma_ac = gamma_ac_ratio * std::sqrt(gamma_ac);
}

void BreastCancer_2D::SetOverDamping(double &gamma_a,
                                        double &gamma_c,
                                        double &gamma_w)
{
    // damping coefficients
    double      gamma_ratio = 2.0;
    double      gamma_aa = Kaa / D0[Na - 1] / D0[Na - 1];
    double      gamma_cc = Kcc / Rc[0] / Rc[0] / 4.0;
    gamma_a = gamma_ratio * std::sqrt(gamma_aa);
    gamma_c = gamma_ratio * std::sqrt(gamma_cc);
    gamma_w = gamma_a;
}

void BreastCancer_2D::Invasion_Growth_Smooth(int R_divide_c,
                                                int R_divide_e,
                                                double gamma_a_ratio,
                                                double gamma_c_ratio,
                                                double gamma_ac_ratio,
                                                double dt,
                                                int save_every_step,
                                                int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_" << Kr_str_c << "_" << Kr_str_ac
        << "_Rd_" << R_divide_c << "E" << R_divide_e 
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio
        << "_Smooth";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      R_divide = (double)R_divide_c * std::pow(10.0, R_divide_e);
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::exponential_distribution<> exp_dist(R_divide); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0);
    for (n = 0; n < Nc; n++){
        t_division[n] = exp_dist(gen_division);
        t_life[n] = uniform_dist(gen_uniform) * t_division[n];
    }

    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15; // increased shell for disk-line interaction

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            t_life[n] += dt;
            if (t_life[n] > t_division[n]){
                // time to divide
                // reset mother cell time
                t_life[n] = 0.0;
                t_division[n] = exp_dist(gen_division);
                // generate daughter cell
                t_life.push_back(0.0);
                t_division.push_back(exp_dist(gen_division));
                while (true){
                    Rc_temp = 1 + normal_dist(gen_r);
                    if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                        Rc.push_back(Rc_temp * Rs); // Rs set in the header
                        break;
                    }
                }
                divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                Vx_c.push_back(0.0);
                Vy_c.push_back(0.0);
                Fx_c.push_back(0.0);
                Fy_c.push_back(0.0);
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth(posname, velname, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_Growth_Smooth(int periphery_size,
                                                        int T_divide_c,
                                                        int T_divide_e,
                                                        double gamma_a_ratio,
                                                        double gamma_c_ratio,
                                                        double gamma_ac_ratio,
                                                        double dt,
                                                        int save_every_step,
                                                        int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_" << periphery_size
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio
        << "_Smooth";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, -1.0), t_life(Nc, 0.0);

    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<int>             is_in_periphery(Nc, 0);
    double                  r_cut_pl = 2.0 * periphery_size * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15; // increased shell for disk-line interaction

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

    Periphery_List_Initial(r_cut_pl, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell, check if close to periphery and reset time
                    t_life[n] = 0.0;
                    Periphery_List_Check(r_cut_pl, is_in_periphery, n);
                    if (is_in_periphery[n] == 1)
                        t_division[n] = T_divide * gamma_dist(gen_division);

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);
            is_in_periphery.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);
            for (n = Nc_old; n < Nc; n++){
                Periphery_List_Check(r_cut_pl, is_in_periphery, n);
                if (is_in_periphery[n] == 1)
                    t_division[n] = T_divide * gamma_dist(gen_division);
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth(posname, velname, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth_Smooth(double phi_c,
                                                            int T_divide_c,
                                                            int T_divide_e,
                                                            double gamma_a_ratio,
                                                            double gamma_c_ratio,
                                                            double gamma_ac_ratio,
                                                            double dt,
                                                            int save_every_step,
                                                            int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio
        << "_Smooth";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, -1.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(100 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(200 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15; // increased shell for disk-line interaction

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell, check if close to periphery and reset time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth(posname, velname, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth_OD_Smooth(double phi_c,
                                                                int T_divide_c,
                                                                int T_divide_e,
                                                                double dt,
                                                                int save_every_step,
                                                                int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_OD_Smooth";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, -1.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    double      gamma_a, gamma_c, gamma_w;
    SetOverDamping(gamma_a, gamma_c, gamma_w);

    double                  Ft = Pt * Ly;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(100 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(200 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15; // increased shell for disk-line interaction

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			xa[n] += dt * Fx_a[n] / gamma_a;
			ya[n] += dt * Fy_a[n] / gamma_a;
		}
		for (n = 0; n < Nc; n++){
            xc[n] += dt * Fx_c[n] / gamma_c;
            yc[n] += dt * Fy_c[n] / gamma_c; 
		}
        Lx += dt * (Fw - Ft) / gamma_w;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell, check if close to periphery and reset time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		if (nt % save_every_step == 0)
			SaveInv_Growth_OD(posname);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Growth(int R_divide_c,
                                        int R_divide_e,
                                        double gamma_a_ratio,
                                        double gamma_c_ratio,
                                        double gamma_ac_ratio,
                                        double dt,
                                        int save_every_step,
                                        int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_" << Kr_str_c << "_" << Kr_str_ac
        << "_Rd_" << R_divide_c << "E" << R_divide_e 
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio;
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      R_divide = (double)R_divide_c * std::pow(10.0, R_divide_e);
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::exponential_distribution<> exp_dist(R_divide); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0);
    for (n = 0; n < Nc; n++){
        t_division[n] = exp_dist(gen_division);
        t_life[n] = uniform_dist(gen_uniform) * t_division[n];
    }

    // damping coefficients
    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0;

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            t_life[n] += dt;
            if (t_life[n] > t_division[n]){
                // time to divide
                // reset mother cell time
                t_life[n] = 0.0;
                t_division[n] = exp_dist(gen_division);
                // generate daughter cell
                t_life.push_back(0.0);
                t_division.push_back(exp_dist(gen_division));
                while (true){
                    Rc_temp = 1 + normal_dist(gen_r);
                    if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                        Rc.push_back(Rc_temp * Rs); // Rs set in the header
                        break;
                    }
                }
                divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                Vx_c.push_back(0.0);
                Vy_c.push_back(0.0);
                Fx_c.push_back(0.0);
                Fy_c.push_back(0.0);
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth(posname, velname, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_Growth(int periphery_size,
                                                int T_divide_c,
                                                int T_divide_e,
                                                double gamma_a_ratio,
                                                double gamma_c_ratio,
                                                double gamma_ac_ratio,
                                                double dt,
                                                int save_every_step,
                                                int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_" << periphery_size
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio;
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0);

    // damping coefficients
    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    vector<int>             is_in_periphery(Nc, 0);
    double                  r_cut_pl = 2.0 * periphery_size * Rs;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0;

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

    Periphery_List_Initial(r_cut_pl, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell time
                    t_life[n] = 0.0;
                    Periphery_List_Check(r_cut_pl, is_in_periphery, n);
                    if (is_in_periphery[n] == 1)
                        t_division[n] = T_divide * gamma_dist(gen_division);
                    
                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);
            is_in_periphery.resize(Nc, 0);
            for (n = Nc_old; n < Nc; n++){
                Periphery_List_Check(r_cut_pl, is_in_periphery, n);
                if (is_in_periphery[n])
                    t_division[n] = T_divide * gamma_dist(gen_division);
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth(posname, velname, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth(double phi_c,
                                                    int T_divide_c,
                                                    int T_divide_e,
                                                    double gamma_a_ratio,
                                                    double gamma_c_ratio,
                                                    double gamma_ac_ratio,
                                                    double dt,
                                                    int save_every_step,
                                                    int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio;
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, -1.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(100 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0; // increased shell for disk-line interaction

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell, check if close to periphery and reset time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth(posname, velname, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth_OD(double phi_c,
                                                        int T_divide_c,
                                                        int T_divide_e,
                                                        double dt,
                                                        int save_every_step,
                                                        int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_OD";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, -1.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    double      gamma_a, gamma_c, gamma_w;
    SetOverDamping(gamma_a, gamma_c, gamma_w);

    double                  Ft = Pt * Ly;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(100 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0; // increased shell for disk-line interaction

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			xa[n] += dt * Fx_a[n] / gamma_a;
			ya[n] += dt * Fy_a[n] / gamma_a;
		}
		for (n = 0; n < Nc; n++){
            xc[n] += dt * Fx_c[n] / gamma_c;
            yc[n] += dt * Fy_c[n] / gamma_c; 
		}
        Lx += dt * (Fw - Ft) / gamma_w;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell, check if close to periphery and reset time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		if (nt % save_every_step == 0)
			SaveInv_Growth_OD(posname);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Growth_LipidLoss(int R_divide_c,
                                                int R_divide_e,
                                                int R_lip_c,
                                                int R_lip_e,
                                                double gamma_a_ratio,
                                                double gamma_c_ratio,
                                                double gamma_ac_ratio,
                                                double dt,
                                                int save_every_step,
                                                int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_" << Kr_str_c << "_" << Kr_str_ac
        << "_Rd_" << R_divide_c << "E" << R_divide_e 
        << "_Rl_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio;
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      R_divide = (double)R_divide_c * std::pow(10.0, R_divide_e);
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::exponential_distribution<> exp_dist(R_divide); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0);
    for (n = 0; n < Nc; n++){
        t_division[n] = exp_dist(gen_division);
        t_life[n] = uniform_dist(gen_uniform) * t_division[n];
    }

    // damping coefficients
    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15;

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            t_life[n] += dt;
            if (t_life[n] > t_division[n]){
                // time to divide
                // reset mother cell time
                t_life[n] = 0.0;
                t_division[n] = exp_dist(gen_division);
                // generate daughter cell
                t_life.push_back(0.0);
                t_division.push_back(exp_dist(gen_division));
                while (true){
                    Rc_temp = 1 + normal_dist(gen_r);
                    if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                        Rc.push_back(Rc_temp * Rs); // Rs set in the header
                        break;
                    }
                }
                divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                Vx_c.push_back(0.0);
                Vy_c.push_back(0.0);
                Fx_c.push_back(0.0);
                Fy_c.push_back(0.0);
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                A0_temp = A0[n] - (double)CC[n] * d_lip;
                if (A0_temp > A0_min[n]){
                    A0[n] = A0_temp;
                    rsc = std::sqrt(A0_temp / A0_ori[n]);
                    L0[n] = rsc * L0_ori[n];
                    D0[n] = rsc * D0_ori[n];
                    R0[n] = 0.5 * D0[n];
                    Da[n] = rsc * Da_ori[n];
                    rsc_all[n] = rsc;
                }
                else
                    reach_min[n] = 1;
            }
		}

        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss(posname, velname, rsc_all, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Growth_LipidLoss_Smooth(int R_divide_c,
                                                        int R_divide_e,
                                                        int R_lip_c,
                                                        int R_lip_e,
                                                        double gamma_a_ratio,
                                                        double gamma_c_ratio,
                                                        double gamma_ac_ratio,
                                                        double dt,
                                                        int save_every_step,
                                                        int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_" << Kr_str_c << "_" << Kr_str_ac
        << "_Rd_" << R_divide_c << "E" << R_divide_e 
        << "_Rl_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio
        << "_Smooth";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      R_divide = (double)R_divide_c * std::pow(10.0, R_divide_e);
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::exponential_distribution<> exp_dist(R_divide); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0);
    for (n = 0; n < Nc; n++){
        t_division[n] = exp_dist(gen_division);
        t_life[n] = uniform_dist(gen_uniform) * t_division[n];
    }

    // damping coefficients
    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15;

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            t_life[n] += dt;
            if (t_life[n] > t_division[n]){
                // time to divide
                // reset mother cell time
                t_life[n] = 0.0;
                t_division[n] = exp_dist(gen_division);
                // generate daughter cell
                t_life.push_back(0.0);
                t_division.push_back(exp_dist(gen_division));
                while (true){
                    Rc_temp = 1 + normal_dist(gen_r);
                    if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                        Rc.push_back(Rc_temp * Rs); // Rs set in the header
                        break;
                    }
                }
                divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                Vx_c.push_back(0.0);
                Vy_c.push_back(0.0);
                Fx_c.push_back(0.0);
                Fy_c.push_back(0.0);
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                A0_temp = A0[n] - (double)CC[n] * d_lip;
                if (A0_temp > A0_min[n]){
                    A0[n] = A0_temp;
                    rsc = std::sqrt(A0_temp / A0_ori[n]);
                    L0[n] = rsc * L0_ori[n];
                    D0[n] = rsc * D0_ori[n];
                    R0[n] = 0.5 * D0[n];
                    Da[n] = rsc * Da_ori[n];
                    rsc_all[n] = rsc;
                }
                else
                    reach_min[n] = 1;
            }
		}

        Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
        // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss(posname, velname, rsc_all, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_Growth_LipidLoss(int periphery_size,
                                                            int T_divide_c,
                                                            int T_divide_e,
                                                            int R_lip_c,
                                                            int R_lip_e,
                                                            double gamma_a_ratio,
                                                            double gamma_c_ratio,
                                                            double gamma_ac_ratio,
                                                            double dt,
                                                            int save_every_step,
                                                            int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_" << periphery_size
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Rl_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio;
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1.0 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0);

    // damping coefficients
    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    vector<int>             is_in_periphery(Nc, 0);
    double                  r_cut_pl = 2.0 * periphery_size * Rs;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15;

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

    Periphery_List_Initial(r_cut_pl, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell time
                    t_life[n] = 0.0;
                    Periphery_List_Check(r_cut_pl, is_in_periphery, n);
                    if (is_in_periphery[n] == 1)
                        t_division[n] = T_divide * gamma_dist(gen_division);
                    
                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);
            is_in_periphery.resize(Nc, 0);
            for (n = Nc_old; n < Nc; n++){
                Periphery_List_Check(r_cut_pl, is_in_periphery, n);
                if (is_in_periphery[n])
                    t_division[n] = T_divide * gamma_dist(gen_division);
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                A0_temp = A0[n] - (double)CC[n] * d_lip;
                if (A0_temp > A0_min[n]){
                    A0[n] = A0_temp;
                    rsc = std::sqrt(A0_temp / A0_ori[n]);
                    L0[n] = rsc * L0_ori[n];
                    D0[n] = rsc * D0_ori[n];
                    R0[n] = 0.5 * D0[n];
                    Da[n] = rsc * Da_ori[n];
                    rsc_all[n] = rsc;
                }
                else
                    reach_min[n] = 1;
            }
		}

        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss(posname, velname, rsc_all, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_Growth_LipidLoss_Smooth(int periphery_size,
                                                                    int T_divide_c,
                                                                    int T_divide_e,
                                                                    int R_lip_c,
                                                                    int R_lip_e,
                                                                    double gamma_a_ratio,
                                                                    double gamma_c_ratio,
                                                                    double gamma_ac_ratio,
                                                                    double dt,
                                                                    int save_every_step,
                                                                    int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_" << periphery_size
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Rl_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio
        << "_Smooth";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1.0 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0);

    // damping coefficients
    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    vector<int>             is_in_periphery(Nc, 0);
    double                  r_cut_pl = 2.0 * periphery_size * Rs;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15;

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    Periphery_List_Initial(r_cut_pl, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell time
                    t_life[n] = 0.0;
                    Periphery_List_Check(r_cut_pl, is_in_periphery, n);
                    if (is_in_periphery[n] == 1)
                        t_division[n] = T_divide * gamma_dist(gen_division);

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);
            is_in_periphery.resize(Nc, 0);
            for (n = Nc_old; n < Nc; n++){
                Periphery_List_Check(r_cut_pl, is_in_periphery, n);
                if (is_in_periphery[n])
                    t_division[n] = T_divide * gamma_dist(gen_division);
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                A0_temp = A0[n] - (double)CC[n] * d_lip;
                if (A0_temp > A0_min[n]){
                    A0[n] = A0_temp;
                    rsc = std::sqrt(A0_temp / A0_ori[n]);
                    L0[n] = rsc * L0_ori[n];
                    D0[n] = rsc * D0_ori[n];
                    R0[n] = 0.5 * D0[n];
                    Da[n] = rsc * Da_ori[n];
                    rsc_all[n] = rsc;
                }
                else
                    reach_min[n] = 1;
            }
		}

        Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
        // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss(posname, velname, rsc_all, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth_LipidLoss(double phi_c,
                                                                int T_divide_c,
                                                                int T_divide_e,
                                                                int R_lip_c,
                                                                int R_lip_e,
                                                                double gamma_a_ratio,
                                                                double gamma_c_ratio,
                                                                double gamma_ac_ratio,
                                                                double dt,
                                                                int save_every_step,
                                                                int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Rl_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio;
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1.0 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    // damping coefficients
    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(100 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
    
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0;

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                A0_temp = A0[n] - (double)CC[n] * d_lip;
                if (A0_temp > A0_min[n]){
                    A0[n] = A0_temp;
                    rsc = std::sqrt(A0_temp / A0_ori[n]);
                    L0[n] = rsc * L0_ori[n];
                    D0[n] = rsc * D0_ori[n];
                    R0[n] = 0.5 * D0[n];
                    Da[n] = rsc * Da_ori[n];
                    rsc_all[n] = rsc;
                }
                else
                    reach_min[n] = 1;
            }
		}

        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
        // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss(posname, velname, rsc_all, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth_LipidLoss_Smooth(double phi_c,
                                                                    int T_divide_c,
                                                                    int T_divide_e,
                                                                    int R_lip_c,
                                                                    int R_lip_e,
                                                                    double gamma_a_ratio,
                                                                    double gamma_c_ratio,
                                                                    double gamma_ac_ratio,
                                                                    double dt,
                                                                    int save_every_step,
                                                                    int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Rl_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_ga_" << std::fixed << std::setprecision(2) << gamma_a_ratio
        << "_gc_" << std::fixed << std::setprecision(2) << gamma_c_ratio
        << "_gac_" << std::fixed << std::setprecision(2) << gamma_ac_ratio
        << "_Smooth";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_inv << "/Vel_Inv_All.txt";
    velname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);
    if ((otptfl = fopen(velname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", velname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1.0 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    // damping coefficients
    double      gamma_a, gamma_c, gamma_ac;
    SetDamping(gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio, gamma_a, gamma_c, gamma_ac);

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = std::sqrt(Na) * Ns[0] / 4.0;
    double                  Ft = Pt * Ly;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(100 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
    
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15; // increased shell for disk-line interaction

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
			Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Lx += dt * Vw;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                A0_temp = A0[n] - (double)CC[n] * d_lip;
                if (A0_temp > A0_min[n]){
                    A0[n] = A0_temp;
                    rsc = std::sqrt(A0_temp / A0_ori[n]);
                    L0[n] = rsc * L0_ori[n];
                    D0[n] = rsc * D0_ori[n];
                    R0[n] = 0.5 * D0[n];
                    Da[n] = rsc * Da_ori[n];
                    rsc_all[n] = rsc;
                }
                else
                    reach_min[n] = 1;
            }
		}

        Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
        // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss(posname, velname, rsc_all, Vx_a, Vy_a, Vx_c, Vy_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth_LipidLoss_OD(double phi_c,
                                                                int T_divide_c,
                                                                int T_divide_e,
                                                                int R_lip_c,
                                                                int R_lip_e,
                                                                double dt,
                                                                int save_every_step,
                                                                int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Rl_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_OD";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1.0 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    // damping coefficients
    double      gamma_a, gamma_c, gamma_w;
    SetOverDamping(gamma_a, gamma_c, gamma_w);

    double                  Ft = Pt * Ly;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(100 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
    
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0;

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			xa[n] += dt * Fx_a[n] / gamma_a;
			ya[n] += dt * Fy_a[n] / gamma_a;
		}
		for (n = 0; n < Nc; n++){
            xc[n] += dt * Fx_c[n] / gamma_c;
            yc[n] += dt * Fy_c[n] / gamma_c; 
		}
        Lx += dt * (Fw - Ft) / gamma_w;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                A0_temp = A0[n] - (double)CC[n] * d_lip;
                if (A0_temp > A0_min[n]){
                    A0[n] = A0_temp;
                    rsc = std::sqrt(A0_temp / A0_ori[n]);
                    L0[n] = rsc * L0_ori[n];
                    D0[n] = rsc * D0_ori[n];
                    R0[n] = 0.5 * D0[n];
                    Da[n] = rsc * Da_ori[n];
                    rsc_all[n] = rsc;
                }
                else
                    reach_min[n] = 1;
            }
		}

        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
        // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss_OD(posname, rsc_all);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth_LipidLoss_OD_Smooth(double phi_c,
                                                                        int T_divide_c,
                                                                        int T_divide_e,
                                                                        int R_lip_c,
                                                                        int R_lip_e,
                                                                        double dt,
                                                                        int save_every_step,
                                                                        int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Rl_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_OD_Smooth";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1.0 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    // damping coefficients
    double      gamma_a, gamma_c, gamma_w;
    SetOverDamping(gamma_a, gamma_c, gamma_w);

    double                  Ft = Pt * Ly;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(100 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
    
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15; // increased shell for disk-line interaction

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			xa[n] += dt * Fx_a[n] / gamma_a;
			ya[n] += dt * Fy_a[n] / gamma_a;
		}
		for (n = 0; n < Nc; n++){
            xc[n] += dt * Fx_c[n] / gamma_c;
            yc[n] += dt * Fy_c[n] / gamma_c; 
		}
        Lx += dt * (Fw - Ft) / gamma_w;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                A0_temp = A0[n] - (double)CC[n] * d_lip;
                if (A0_temp > A0_min[n]){
                    A0[n] = A0_temp;
                    rsc = std::sqrt(A0_temp / A0_ori[n]);
                    L0[n] = rsc * L0_ori[n];
                    D0[n] = rsc * D0_ori[n];
                    R0[n] = 0.5 * D0[n];
                    Da[n] = rsc * Da_ori[n];
                    rsc_all[n] = rsc;
                }
                else
                    reach_min[n] = 1;
            }
		}

        Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
        // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss_OD(posname, rsc_all);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth_LipidLossLinear_OD(double phi_c,
                                                                        int T_divide_c,
                                                                        int T_divide_e,
                                                                        int R_lip_c,
                                                                        int R_lip_e,
                                                                        double dt,
                                                                        int save_every_step,
                                                                        int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Rl_Linear_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_OD";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1.0 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    // damping coefficients
    double      gamma_a, gamma_c, gamma_w;
    SetOverDamping(gamma_a, gamma_c, gamma_w);

    double                  Ft = Pt * Ly;

    // lipid loss parameters
    int                     CC_thresh = 10;
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), is_ll(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(100 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
    
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0;

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			xa[n] += dt * Fx_a[n] / gamma_a;
			ya[n] += dt * Fy_a[n] / gamma_a;
		}
		for (n = 0; n < Nc; n++){
            xc[n] += dt * Fx_c[n] / gamma_c;
            yc[n] += dt * Fy_c[n] / gamma_c; 
		}
        Lx += dt * (Fw - Ft) / gamma_w;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min, is_ll);
		for (n = 0; n < Na; n++){
            if (reach_min[n] == 0){
                if ((is_ll[n] == 1) || (CC[n] > CC_thresh)){
                    if (is_ll[n] == 0)
                        is_ll[n] = 1;
                    A0_temp = A0[n] - d_lip;
                    if (A0_temp > A0_min[n]){
                        A0[n] = A0_temp;
                        rsc = std::sqrt(A0_temp / A0_ori[n]);
                        L0[n] = rsc * L0_ori[n];
                        D0[n] = rsc * D0_ori[n];
                        R0[n] = 0.5 * D0[n];
                        Da[n] = rsc * Da_ori[n];
                        rsc_all[n] = rsc;
                    }
                    else
                        reach_min[n] = 1;
                }
            }
		}

        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
        // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss_OD(posname, rsc_all);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_2D::Invasion_Periphery_PF_Growth_LipidLossLinear_OD_Smooth(double phi_c,
                                                                                int T_divide_c,
                                                                                int T_divide_e,
                                                                                int R_lip_c,
                                                                                int R_lip_e,
                                                                                double dt,
                                                                                int save_every_step,
                                                                                int Nt)
{
    string      dir_inv, posname, velname;
    string      cmd;

    std::ostringstream oss;
    oss << dir_run << "Growth_Periphery_PF_"
        << std::fixed << std::setprecision(2) << phi_c
        << "_" << Kr_str_c << "_" << Kr_str_ac
        << "_Td_" << T_divide_c << "E" << T_divide_e 
        << "_Rl_Linear_" << R_lip_c << "E" << R_lip_e
        << "_Kpin_" << std::fixed << std::setprecision(4) << Kpin
        << "_OD_Smooth";
    dir_inv = oss.str();
    oss.str("");
    oss.clear();

	oss << "mkdir -p " << dir_inv;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_inv << "/Pos_Inv_All.txt";
    posname = oss.str();
    oss.str("");
    oss.clear();

    FILE*    otptfl;
    if ((otptfl = fopen(posname.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s\n", posname.c_str());
		exit(1);
	}
    fclose(otptfl);

    LoadConfig_All(all_name);

    int 				    n, nt = 0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle;
    std::mt19937 gen_division(pid);
    std::mt19937 gen_uniform(pid + 1000);
    std::mt19937 gen_r(pid + 2000);
    std::gamma_distribution<> gamma_dist(1.0 / T_divide_std / T_divide_std, T_divide_std * T_divide_std); // division time
    std::uniform_real_distribution<> uniform_dist(0.0, 1.0);
    std::normal_distribution<> normal_dist(0.0, 0.15); // radius for new cancer cells
    vector<double>  t_division(Nc, 0.0), t_life(Nc, 0.0), phi(Nc, 0.0);
    vector<int>     is_in_periphery(Nc, 0), is_dividing(Nc, 0), is_not_dividing_ori(Nc, 1);

    // damping coefficients
    double      gamma_a, gamma_c, gamma_w;
    SetOverDamping(gamma_a, gamma_c, gamma_w);

    double                  Ft = Pt * Ly;

    // lipid loss parameters
    int                     CC_thresh = 10;
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  A0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), is_ll(Na, 0), reach_min(Na, 0);
    vector<double>          Da_ori(Na, 0.0), D0_ori(Na, 0.0), R0_ori(Na, 0.0), L0_ori(Na, 0.0), A0_ori(Na, 0.0), A0_min(Na, 0.0), rsc_all(Na, 1.0);
    for (n = 0; n < Na; n++){
        Da_ori[n] = Da[n];
        D0_ori[n] = D0[n];
        R0_ori[n] = R0[n];
        L0_ori[n] = L0[n];
        A0_ori[n] = A0[n];
        A0_min[n] = rsc_min * rsc_min * A0[n];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(40 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(200 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(40 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_phi_cc(100 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_phi_count;
    
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0);
    vector<double>          xc_phi_save(Nc, 0.0), yc_phi_save(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_phi_cut = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
    r_cut *= a_cut_ac;
	r_cut = (r_cut > a_cut_c * Rs_ratio_max * Rs) ? r_cut : a_cut_c * Rs_ratio_max * Rs;
    r_cut *= 2.0 * 1.15; // increased shell for disk-line interaction

	// Verlet list initialization
	VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
    LocalPhi(xc, yc, r_shell, VL_phi_cc, VL_phi_count, phi);
    
    Periphery_List_PF_Initial(phi_c, phi, is_in_periphery);
    for (n = 0; n < Nc; n++){
        if (is_in_periphery[n]){
            is_dividing[n] = 1;
            is_not_dividing_ori[n] = 0;
            t_division[n] = T_divide * gamma_dist(gen_division);
        }
    }

	first_call = 0;

	for (nt = 1; nt <= Nt; nt++){
		for (n = 0; n < Ns_tot; n++){
			xa[n] += dt * Fx_a[n] / gamma_a;
			ya[n] += dt * Fy_a[n] / gamma_a;
		}
		for (n = 0; n < Nc; n++){
            xc[n] += dt * Fx_c[n] / gamma_c;
            yc[n] += dt * Fy_c[n] / gamma_c; 
		}
        Lx += dt * (Fw - Ft) / gamma_w;

        // Check for division
        for (n = 0; n < Nc; n++){
            if (is_in_periphery[n]){
                t_life[n] += dt;
                if (t_life[n] > t_division[n]){
                    // time to divide
                    // reset mother cell time
                    t_life[n] = 0.0;
                    is_in_periphery[n] = 0;
                    is_dividing[n] = 0;

                    while (true){
                        Rc_temp = 1 + normal_dist(gen_r);
                        if ((Rc_temp >= Rs_ratio_min) && (Rc_temp <= Rs_ratio_max)){
                            Rc.push_back(Rc_temp * Rs); // Rs set in the header
                            Dc.push_back(2.0 * Rc_temp * Rs);
                            Rc_sq.push_back(Rc_temp * Rs * Rc_temp * Rs);
                            break;
                        }
                    }
                    divide_angle = 2.0 * PI * uniform_dist(gen_uniform);
                    xc.push_back(xc[n] + Rc[n] * std::cos(divide_angle));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle));
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save.resize(Nc, 0.0);
            yc_save.resize(Nc, 0.0);
	        VL_cc.resize(40 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut, xc_save, yc_save, VL_cc, VL_ac, VL_count_all);

            xc_phi_save.resize(Nc, 0.0);
            yc_phi_save.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_phi_cc.resize(100 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count);
            VerletList_Disk(r_phi_cut, xc_phi_save, yc_phi_save, VL_phi_cc, VL_phi_count, first_call);
            LocalPhi_Update(xc, yc, r_shell, VL_phi_cc, VL_phi_count, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM_Disk(r_cut, xa_save, ya_save, xc_save, yc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min, is_ll);
		for (n = 0; n < Na; n++){
            if (reach_min[n] == 0){
                if ((is_ll[n] == 1) || (CC[n] > CC_thresh)){
                    if (is_ll[n] == 0)
                        is_ll[n] = 1;
                    A0_temp = A0[n] - d_lip;
                    if (A0_temp > A0_min[n]){
                        A0[n] = A0_temp;
                        rsc = std::sqrt(A0_temp / A0_ori[n]);
                        L0[n] = rsc * L0_ori[n];
                        D0[n] = rsc * D0_ori[n];
                        R0[n] = 0.5 * D0[n];
                        Da[n] = rsc * Da_ori[n];
                        rsc_all[n] = rsc;
                    }
                    else
                        reach_min[n] = 1;
                }
            }
		}

        Force_DPM_Disk_Smooth_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
        // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss_OD(posname, rsc_all);

        if (CheckInvasionBoundary())
            break;
	}
}
