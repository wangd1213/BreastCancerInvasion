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

void BreastCancer_3D::SetDamping(double gamma_a_ratio,
                                    double gamma_c_ratio,
                                    double gamma_ac_ratio,
                                    double &gamma_a,
                                    double &gamma_c,
                                    double &gamma_ac)
{
    // damping coefficients
    double                  gamma_aa = Kaa / R0[Na - 1] / R0[Na - 1] / 4.0;
    double                  gamma_cc = Kcc / Rs / Rs / 4.0;
    gamma_ac = Kac / (R0[Na - 1] + Rs) / (R0[Na - 1] + Rs);
    gamma_a = gamma_a_ratio * std::sqrt(gamma_aa);
    gamma_c = gamma_c_ratio * std::sqrt(gamma_cc);
    gamma_ac = gamma_ac_ratio * std::sqrt(gamma_ac);
}

void BreastCancer_3D::Invasion_Periphery_PF_Growth(double phi_c,
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
    double      Rc_temp, divide_angle_theta, divide_angle_phi;
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

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0), Vz_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0), Vz_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = 300.0;
    double                  Ft = Pt * Lx * Lz;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(400 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(800 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(200 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_cc_phi(400 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_count_aa, VL_count_cc, VL_count_ac, VL_count_phi;
	vector<double>			xa_save_aa(Ns_tot, 0.0), ya_save_aa(Ns_tot, 0.0), za_save_aa(Ns_tot, 0.0);
    vector<double>          xc_save_cc(Nc, 0.0), yc_save_cc(Nc, 0.0), zc_save_cc(Nc, 0.0);
    vector<double>			xa_save_ac(Ns_tot, 0.0), ya_save_ac(Ns_tot, 0.0), za_save_ac(Ns_tot, 0.0);
    vector<double>          xc_save_ac(Nc, 0.0), yc_save_ac(Nc, 0.0), zc_save_ac(Nc, 0.0);
    vector<double>          xc_save_phi(Nc, 0.0), yc_save_phi(Nc, 0.0), zc_save_phi(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_cut_phi = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut_aa = R0[0], r_cut_cc = Rc[0];
	for (n = 1; n < Na; n++)
		r_cut_aa = (r_cut_aa > R0[n]) ? r_cut_aa : R0[n];
    for (n = 1; n < Nc; n++)
        r_cut_cc = (r_cut_cc > Rc[n]) ? r_cut_cc : Rc[n];
    r_cut_cc *= 2.0 * a_cut_c;
    r_cut_aa *= 2.0;
    double                  r_cut_ac = (r_cut_aa > r_cut_cc) ? r_cut_aa : r_cut_cc;

	// Verlet list initialization
    VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
    VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
	VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
    VL_count_all[0] = VL_count_aa;
    VL_count_all[1] = VL_count_cc;
    VL_count_all[2] = VL_count_ac;
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
    LocalPhi(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, phi);
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
            Vz_a[n] += dt_half * (Fz_a[n] - gamma_a * Vz_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
            za[n] += dt * Vz_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            Vz_c[n] += dt_half * (Fz_c[n] - gamma_c * Vz_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
            zc[n] += dt * Vz_c[n];
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Ly += dt * Vw;

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
                            Dc_cube.push_back(8.0 * Rc_temp * Rc_temp * Rc_temp * Rs * Rs * Rs);
                            break;
                        }
                    }
                    divide_angle_theta = 2.0 * PI * uniform_dist(gen_uniform);
                    divide_angle_phi = std::acos(2.0 * uniform_dist(gen_uniform) - 1.0);
                    xc.push_back(xc[n] + Rc[n] * std::sin(divide_angle_phi) * std::cos(divide_angle_theta));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle_phi) * std::sin(divide_angle_theta));
                    zc.push_back(zc[n] + Rc[n] * std::cos(divide_angle_phi));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Vz_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                    Fz_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save_cc.resize(Nc, 0.0);
            yc_save_cc.resize(Nc, 0.0);
            zc_save_cc.resize(Nc, 0.0);
            xc_save_ac.resize(Nc, 0.0);
            yc_save_ac.resize(Nc, 0.0);
            zc_save_ac.resize(Nc, 0.0);
	        VL_cc.resize(200 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac);
            VerletList_Disk_Append(Nc_old, r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc);
            VL_count_all[1] = VL_count_cc;
            VL_count_all[2] = VL_count_ac;

            xc_save_phi.resize(Nc, 0.0);
            yc_save_phi.resize(Nc, 0.0);
            zc_save_phi.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_cc_phi.resize(400 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi);
            VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
            LocalPhi_Update(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
        VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
        VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
        VL_count_all[0] = VL_count_aa;
        VL_count_all[1] = VL_count_cc;
        VL_count_all[2] = VL_count_ac;
        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * (Fx_a[n] - gamma_a * Vx_a[n]);
            Vy_a[n] += dt_half * (Fy_a[n] - gamma_a * Vy_a[n]);
            Vz_a[n] += dt_half * (Fz_a[n] - gamma_a * Vz_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            Vz_c[n] += dt_half * (Fz_c[n] - gamma_c * Vz_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth(posname, velname, Vx_a, Vy_a, Vz_a, Vx_c, Vy_c, Vz_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_3D::Invasion_Periphery_PF_Growth_LipidLoss(double phi_c,
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

    int 				    n, nf, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle_theta, divide_angle_phi;
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

	vector<double> 		    Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0), Vz_a(Ns_tot, 0.0);
	vector<double>	 	    Vx_c(Nc, 0.0), Vy_c(Nc, 0.0), Vz_c(Nc, 0.0);
    Vw = 0.0;
    double                  mw = 300.0;
    double                  Ft = Pt * Lx * Lz;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  V0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Ra_ori(Na, 0.0), R0_ori(Na, 0.0), D0_intra_ori(Na, 0.0), D0_intra_sq_ori(Na, 0.0), V0_ori(Na, 0.0), V0_min(Na, 0.0), rsc_all(Na, 1.0);
    vector<vector <double>> A0_ori(Na, vector<double> (Nf, 0.0));
    for (n = 0; n < Na; n++){
        V0_ori[n] = V0[n];
        Ra_ori[n] = Ra[n];
        R0_ori[n] = R0[n];
        D0_intra_ori[n] = D0_intra[n];
        D0_intra_sq_ori[n] = D0_intra_sq[n];
        V0_min[n] = rsc_min * rsc_min * rsc_min * V0[n];
        for (nf = 0; nf < Nf; nf++)
            A0_ori[n][nf] = A0[n][nf];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(400 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(800 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(200 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_cc_phi(400 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_count_aa, VL_count_cc, VL_count_ac, VL_count_phi;
	vector<double>			xa_save_aa(Ns_tot, 0.0), ya_save_aa(Ns_tot, 0.0), za_save_aa(Ns_tot, 0.0);
    vector<double>          xc_save_cc(Nc, 0.0), yc_save_cc(Nc, 0.0), zc_save_cc(Nc, 0.0);
    vector<double>			xa_save_ac(Ns_tot, 0.0), ya_save_ac(Ns_tot, 0.0), za_save_ac(Ns_tot, 0.0);
    vector<double>          xc_save_ac(Nc, 0.0), yc_save_ac(Nc, 0.0), zc_save_ac(Nc, 0.0);
    vector<double>          xc_save_phi(Nc, 0.0), yc_save_phi(Nc, 0.0), zc_save_phi(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_cut_phi = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut_aa = R0[0], r_cut_cc = Rc[0];
	for (n = 1; n < Na; n++)
		r_cut_aa = (r_cut_aa > R0[n]) ? r_cut_aa : R0[n];
    for (n = 1; n < Nc; n++)
        r_cut_cc = (r_cut_cc > Rc[n]) ? r_cut_cc : Rc[n];
    r_cut_cc *= 2.0 * a_cut_c;
    r_cut_aa *= 2.0;
    double                  r_cut_ac = (r_cut_aa > r_cut_cc) ? r_cut_aa : r_cut_cc;

	// Verlet list initialization
	VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
    VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
	VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
    VL_count_all[0] = VL_count_aa;
    VL_count_all[1] = VL_count_cc;
    VL_count_all[2] = VL_count_ac;
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
    LocalPhi(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, phi);
    
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
            Vz_a[n] += dt_half * (Fz_a[n] - gamma_a * Vz_a[n]);
			xa[n] += dt * Vx_a[n];
			ya[n] += dt * Vy_a[n];
            za[n] += dt * Vz_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            Vz_c[n] += dt_half * (Fz_c[n] - gamma_c * Vz_c[n]);
            xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n]; 
            zc[n] += dt * Vz_c[n];
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);
        Ly += dt * Vw;

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
                            Dc_cube.push_back(8.0 * Rc_temp * Rc_temp * Rc_temp * Rs * Rs * Rs);
                            break;
                        }
                    }
                    divide_angle_theta = 2.0 * PI * uniform_dist(gen_uniform);
                    divide_angle_phi = std::acos(2.0 * uniform_dist(gen_uniform) - 1.0);
                    xc.push_back(xc[n] + Rc[n] * std::sin(divide_angle_phi) * std::cos(divide_angle_theta));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle_phi) * std::sin(divide_angle_theta));
                    zc.push_back(zc[n] + Rc[n] * std::cos(divide_angle_phi));
                    Vx_c.push_back(0.0);
                    Vy_c.push_back(0.0);
                    Vz_c.push_back(0.0);
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                    Fz_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save_cc.resize(Nc, 0.0);
            yc_save_cc.resize(Nc, 0.0);
            zc_save_cc.resize(Nc, 0.0);
            xc_save_ac.resize(Nc, 0.0);
            yc_save_ac.resize(Nc, 0.0);
            zc_save_ac.resize(Nc, 0.0);
	        VL_cc.resize(200 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac);
            VerletList_Disk_Append(Nc_old, r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc);
            VL_count_all[1] = VL_count_cc;
            VL_count_all[2] = VL_count_ac;

            xc_save_phi.resize(Nc, 0.0);
            yc_save_phi.resize(Nc, 0.0);
            zc_save_phi.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_cc_phi.resize(400 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi);
            VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
            LocalPhi_Update(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
        VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
        VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
        VL_count_all[0] = VL_count_aa;
        VL_count_all[1] = VL_count_cc;
        VL_count_all[2] = VL_count_ac;

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                V0_temp = V0[n] - (double)CC[n] * d_lip;
                if (V0_temp > V0_min[n]){
                    V0[n] = V0_temp;
                    rsc = std::pow(V0_temp / V0_ori[n], 1.0 / 3.0);
                    R0[n] = rsc * R0_ori[n];
                    D0_intra[n] = rsc * D0_intra_ori[n];
                    D0_intra_sq[n] = rsc * rsc * D0_intra_sq_ori[n];
                    Ra[n] = rsc * Ra_ori[n];
                    rsc_all[n] = rsc;
                    for (nf = 0; nf < Nf; nf++)
                        A0[n][nf] = rsc * rsc * A0_ori[n][nf];
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
            Vz_a[n] += dt_half * (Fz_a[n] - gamma_a * Vz_a[n]);
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * (Fx_c[n] - gamma_c * Vx_c[n]);
			Vy_c[n] += dt_half * (Fy_c[n] - gamma_c * Vy_c[n]);
            Vz_c[n] += dt_half * (Fz_c[n] - gamma_c * Vz_c[n]);
		}
        Vw += dt_half * (Fw - Ft - gamma_a * Vw);

		if (nt % save_every_step == 0)
			SaveInv_Growth_LipidLoss(posname, velname, rsc_all, Vx_a, Vy_a, Vz_a, Vx_c, Vy_c, Vz_c);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_3D::SetOverDamping(double &gamma_a,
                                        double &gamma_c,
                                        double &gamma_w)
{
    // damping coefficients
    double      gamma_ratio = 2.0;
    double      gamma_aa = Kaa / R0[Na - 1] / R0[Na - 1] / 4.0;
    double      gamma_cc = Kcc / Rs / Rs / 4.0;
    gamma_a = gamma_ratio * std::sqrt(gamma_aa);
    gamma_c = gamma_ratio * std::sqrt(gamma_cc);
    gamma_w = gamma_a;
}

void BreastCancer_3D::Invasion_Periphery_PF_Growth_OD(double phi_c,
                                                        int T_divide_c,
                                                        int T_divide_e,
                                                        double dt,
                                                        int save_every_step,
                                                        int Nt)
{
    string      dir_inv, posname;
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
    double      Rc_temp, divide_angle_theta, divide_angle_phi;
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

    double                  Ft = Pt * Lx * Lz;

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(400 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(800 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(200 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_cc_phi(400 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_count_aa, VL_count_cc, VL_count_ac, VL_count_phi;
	vector<double>			xa_save_aa(Ns_tot, 0.0), ya_save_aa(Ns_tot, 0.0), za_save_aa(Ns_tot, 0.0);
    vector<double>          xc_save_cc(Nc, 0.0), yc_save_cc(Nc, 0.0), zc_save_cc(Nc, 0.0);
    vector<double>			xa_save_ac(Ns_tot, 0.0), ya_save_ac(Ns_tot, 0.0), za_save_ac(Ns_tot, 0.0);
    vector<double>          xc_save_ac(Nc, 0.0), yc_save_ac(Nc, 0.0), zc_save_ac(Nc, 0.0);
    vector<double>          xc_save_phi(Nc, 0.0), yc_save_phi(Nc, 0.0), zc_save_phi(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_cut_phi = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut_aa = R0[0], r_cut_cc = Rc[0];
	for (n = 1; n < Na; n++)
		r_cut_aa = (r_cut_aa > R0[n]) ? r_cut_aa : R0[n];
    for (n = 1; n < Nc; n++)
        r_cut_cc = (r_cut_cc > Rc[n]) ? r_cut_cc : Rc[n];
    r_cut_cc *= 2.0 * a_cut_c;
    r_cut_aa *= 2.0;
    double                  r_cut_ac = (r_cut_aa > r_cut_cc) ? r_cut_aa : r_cut_cc;

	// Verlet list initialization
    VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
    VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
	VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
    VL_count_all[0] = VL_count_aa;
    VL_count_all[1] = VL_count_cc;
    VL_count_all[2] = VL_count_ac;
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
    LocalPhi(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, phi);
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
            za[n] += dt * Fz_a[n] / gamma_a;
		}
		for (n = 0; n < Nc; n++){
            xc[n] += dt * Fx_c[n] / gamma_c;
            yc[n] += dt * Fy_c[n] / gamma_c; 
            zc[n] += dt * Fz_c[n] / gamma_c;
		}
        Ly += dt * (Fw - Ft) / gamma_w;

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
                            Dc_cube.push_back(8.0 * Rc_temp * Rc_temp * Rc_temp * Rs * Rs * Rs);
                            break;
                        }
                    }
                    divide_angle_theta = 2.0 * PI * uniform_dist(gen_uniform);
                    divide_angle_phi = std::acos(2.0 * uniform_dist(gen_uniform) - 1.0);
                    xc.push_back(xc[n] + Rc[n] * std::sin(divide_angle_phi) * std::cos(divide_angle_theta));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle_phi) * std::sin(divide_angle_theta));
                    zc.push_back(zc[n] + Rc[n] * std::cos(divide_angle_phi));
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                    Fz_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save_cc.resize(Nc, 0.0);
            yc_save_cc.resize(Nc, 0.0);
            zc_save_cc.resize(Nc, 0.0);
            xc_save_ac.resize(Nc, 0.0);
            yc_save_ac.resize(Nc, 0.0);
            zc_save_ac.resize(Nc, 0.0);
	        VL_cc.resize(200 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac);
            VerletList_Disk_Append(Nc_old, r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc);
            VL_count_all[1] = VL_count_cc;
            VL_count_all[2] = VL_count_ac;

            xc_save_phi.resize(Nc, 0.0);
            yc_save_phi.resize(Nc, 0.0);
            zc_save_phi.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_cc_phi.resize(400 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi);
            VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
            LocalPhi_Update(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
        VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
        VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
        VL_count_all[0] = VL_count_aa;
        VL_count_all[1] = VL_count_cc;
        VL_count_all[2] = VL_count_ac;
        Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		if (nt % save_every_step == 0)
			SaveInv_Growth_OD(posname);

        if (CheckInvasionBoundary())
            break;
	}
}

void BreastCancer_3D::Invasion_Periphery_PF_Growth_LipidLoss_OD(double phi_c,
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

    int 				    n, nf, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle_theta, divide_angle_phi;
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

    double                  Ft = Pt * Lx * Lz;

    // lipid loss parameters
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  V0_temp, rsc, rsc_min = 0.1;
    vector<int>             CC(Na, 0), reach_min(Na, 0);
    vector<double>          Ra_ori(Na, 0.0), R0_ori(Na, 0.0), D0_intra_ori(Na, 0.0), D0_intra_sq_ori(Na, 0.0), V0_ori(Na, 0.0), V0_min(Na, 0.0), rsc_all(Na, 1.0);
    vector<vector <double>> A0_ori(Na, vector<double> (Nf, 0.0));
    for (n = 0; n < Na; n++){
        V0_ori[n] = V0[n];
        Ra_ori[n] = Ra[n];
        R0_ori[n] = R0[n];
        D0_intra_ori[n] = D0_intra[n];
        D0_intra_sq_ori[n] = D0_intra_sq[n];
        V0_min[n] = rsc_min * rsc_min * rsc_min * V0[n];
        for (nf = 0; nf < Nf; nf++)
            A0_ori[n][nf] = A0[n][nf];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(400 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(800 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(200 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_cc_phi(400 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_count_aa, VL_count_cc, VL_count_ac, VL_count_phi;
	vector<double>			xa_save_aa(Ns_tot, 0.0), ya_save_aa(Ns_tot, 0.0), za_save_aa(Ns_tot, 0.0);
    vector<double>          xc_save_cc(Nc, 0.0), yc_save_cc(Nc, 0.0), zc_save_cc(Nc, 0.0);
    vector<double>			xa_save_ac(Ns_tot, 0.0), ya_save_ac(Ns_tot, 0.0), za_save_ac(Ns_tot, 0.0);
    vector<double>          xc_save_ac(Nc, 0.0), yc_save_ac(Nc, 0.0), zc_save_ac(Nc, 0.0);
    vector<double>          xc_save_phi(Nc, 0.0), yc_save_phi(Nc, 0.0), zc_save_phi(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_cut_phi = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut_aa = R0[0], r_cut_cc = Rc[0];
	for (n = 1; n < Na; n++)
		r_cut_aa = (r_cut_aa > R0[n]) ? r_cut_aa : R0[n];
    for (n = 1; n < Nc; n++)
        r_cut_cc = (r_cut_cc > Rc[n]) ? r_cut_cc : Rc[n];
    r_cut_cc *= 2.0 * a_cut_c;
    r_cut_aa *= 2.0;
    double                  r_cut_ac = (r_cut_aa > r_cut_cc) ? r_cut_aa : r_cut_cc;

	// Verlet list initialization
	VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
    VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
	VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
    VL_count_all[0] = VL_count_aa;
    VL_count_all[1] = VL_count_cc;
    VL_count_all[2] = VL_count_ac;
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
    LocalPhi(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, phi);
    
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
            za[n] += dt * Fz_a[n] / gamma_a;
		}
		for (n = 0; n < Nc; n++){
            xc[n] += dt * Fx_c[n] / gamma_c;
            yc[n] += dt * Fy_c[n] / gamma_c; 
            zc[n] += dt * Fz_c[n] / gamma_c;
		}
        Ly += dt * (Fw - Ft) / gamma_w;

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
                            Dc_cube.push_back(8.0 * Rc_temp * Rc_temp * Rc_temp * Rs * Rs * Rs);
                            break;
                        }
                    }
                    divide_angle_theta = 2.0 * PI * uniform_dist(gen_uniform);
                    divide_angle_phi = std::acos(2.0 * uniform_dist(gen_uniform) - 1.0);
                    xc.push_back(xc[n] + Rc[n] * std::sin(divide_angle_phi) * std::cos(divide_angle_theta));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle_phi) * std::sin(divide_angle_theta));
                    zc.push_back(zc[n] + Rc[n] * std::cos(divide_angle_phi));
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                    Fz_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save_cc.resize(Nc, 0.0);
            yc_save_cc.resize(Nc, 0.0);
            zc_save_cc.resize(Nc, 0.0);
            xc_save_ac.resize(Nc, 0.0);
            yc_save_ac.resize(Nc, 0.0);
            zc_save_ac.resize(Nc, 0.0);
	        VL_cc.resize(200 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac);
            VerletList_Disk_Append(Nc_old, r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc);
            VL_count_all[1] = VL_count_cc;
            VL_count_all[2] = VL_count_ac;

            xc_save_phi.resize(Nc, 0.0);
            yc_save_phi.resize(Nc, 0.0);
            zc_save_phi.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_cc_phi.resize(400 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi);
            VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
            LocalPhi_Update(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
        VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
        VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
        VL_count_all[0] = VL_count_aa;
        VL_count_all[1] = VL_count_cc;
        VL_count_all[2] = VL_count_ac;

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min);
		for (n = 0; n < Na; n++){
            if ((reach_min[n] == 0) && (CC[n] > 0)){
                V0_temp = V0[n] - (double)CC[n] * d_lip;
                if (V0_temp > V0_min[n]){
                    V0[n] = V0_temp;
                    rsc = std::pow(V0_temp / V0_ori[n], 1.0 / 3.0);
                    R0[n] = rsc * R0_ori[n];
                    D0_intra[n] = rsc * D0_intra_ori[n];
                    D0_intra_sq[n] = rsc * rsc * D0_intra_sq_ori[n];
                    Ra[n] = rsc * Ra_ori[n];
                    rsc_all[n] = rsc;
                    for (nf = 0; nf < Nf; nf++)
                        A0[n][nf] = rsc * rsc * A0_ori[n][nf];
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

void BreastCancer_3D::Invasion_Periphery_PF_Growth_LipidLossLinear_OD(double phi_c,
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

    int 				    n, nf, nt = 0;
	double 				    dt_half = dt / 2.0;

    // division setup
    int         Nc_old = Nc;
    double      T_divide = (double)T_divide_c * std::pow(10.0, T_divide_e);
    double      T_divide_std = 0.4;
    double      Rc_temp, divide_angle_theta, divide_angle_phi;
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

    double                  Ft = Pt * Lx * Lz;

    // lipid loss parameters
    int                     CC_thresh = 10;
    double                  d_lip = (double)R_lip_c * std::pow(10.0, (double)R_lip_e) * dt;
    double                  V0_temp, rsc, rsc_min = 0.1;
    vector<int>             is_ll(Na, 0), CC(Na, 0), reach_min(Na, 0);
    vector<double>          Ra_ori(Na, 0.0), R0_ori(Na, 0.0), D0_intra_ori(Na, 0.0), D0_intra_sq_ori(Na, 0.0), V0_ori(Na, 0.0), V0_min(Na, 0.0), rsc_all(Na, 1.0);
    vector<vector <double>> A0_ori(Na, vector<double> (Nf, 0.0));
    for (n = 0; n < Na; n++){
        V0_ori[n] = V0[n];
        Ra_ori[n] = Ra[n];
        R0_ori[n] = R0[n];
        D0_intra_ori[n] = D0_intra[n];
        D0_intra_sq_ori[n] = D0_intra_sq[n];
        V0_min[n] = rsc_min * rsc_min * rsc_min * V0[n];
        for (nf = 0; nf < Nf; nf++)
            A0_ori[n][nf] = A0[n][nf];
    }

	// Verlet list parameters
    int                     first_call = 1;
	vector<vector <int>>	VL_aa(400 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(800 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(200 * Nc, vector<int> (2, -1));
    vector<vector <int>>    VL_cc_phi(400 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
    int                     VL_count_aa, VL_count_cc, VL_count_ac, VL_count_phi;
	vector<double>			xa_save_aa(Ns_tot, 0.0), ya_save_aa(Ns_tot, 0.0), za_save_aa(Ns_tot, 0.0);
    vector<double>          xc_save_cc(Nc, 0.0), yc_save_cc(Nc, 0.0), zc_save_cc(Nc, 0.0);
    vector<double>			xa_save_ac(Ns_tot, 0.0), ya_save_ac(Ns_tot, 0.0), za_save_ac(Ns_tot, 0.0);
    vector<double>          xc_save_ac(Nc, 0.0), yc_save_ac(Nc, 0.0), zc_save_ac(Nc, 0.0);
    vector<double>          xc_save_phi(Nc, 0.0), yc_save_phi(Nc, 0.0), zc_save_phi(Nc, 0.0);
    double                  r_shell_ratio = 2.0;
    double                  r_shell = r_shell_ratio * 2.0 * Rs_ratio_max * Rs;
    double                  r_cut_phi = (r_shell_ratio + 1.0) * 2.0 * Rs_ratio_max * Rs;
	double 		            r_cut_aa = R0[0], r_cut_cc = Rc[0];
	for (n = 1; n < Na; n++)
		r_cut_aa = (r_cut_aa > R0[n]) ? r_cut_aa : R0[n];
    for (n = 1; n < Nc; n++)
        r_cut_cc = (r_cut_cc > Rc[n]) ? r_cut_cc : Rc[n];
    r_cut_cc *= 2.0 * a_cut_c;
    r_cut_aa *= 2.0;
    double                  r_cut_ac = (r_cut_aa > r_cut_cc) ? r_cut_aa : r_cut_cc;

	// Verlet list initialization
	VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
    VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
	VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
    VL_count_all[0] = VL_count_aa;
    VL_count_all[1] = VL_count_cc;
    VL_count_all[2] = VL_count_ac;
	Force_DPM_Disk_PinCen_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
    // Force_DPM_Disk_Smooth_Rigid_PinCen_LipidLoss_VL(VL_aa, VL_cc, VL_ac, VL_count_all, reach_min, rsc_all, rsc_min);
    double  Pt = Pressure_DPM_Disk();

    VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
    LocalPhi(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, phi);
    
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
            za[n] += dt * Fz_a[n] / gamma_a;
		}
		for (n = 0; n < Nc; n++){
            xc[n] += dt * Fx_c[n] / gamma_c;
            yc[n] += dt * Fy_c[n] / gamma_c; 
            zc[n] += dt * Fz_c[n] / gamma_c;
		}
        Ly += dt * (Fw - Ft) / gamma_w;

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
                            Dc_cube.push_back(8.0 * Rc_temp * Rc_temp * Rc_temp * Rs * Rs * Rs);
                            break;
                        }
                    }
                    divide_angle_theta = 2.0 * PI * uniform_dist(gen_uniform);
                    divide_angle_phi = std::acos(2.0 * uniform_dist(gen_uniform) - 1.0);
                    xc.push_back(xc[n] + Rc[n] * std::sin(divide_angle_phi) * std::cos(divide_angle_theta));
                    yc.push_back(yc[n] + Rc[n] * std::sin(divide_angle_phi) * std::sin(divide_angle_theta));
                    zc.push_back(zc[n] + Rc[n] * std::cos(divide_angle_phi));
                    Fx_c.push_back(0.0);
                    Fy_c.push_back(0.0);
                    Fz_c.push_back(0.0);
                }
            }
        }
        Nc = int(xc.size());
        if (Nc_old < Nc){
            xc_save_cc.resize(Nc, 0.0);
            yc_save_cc.resize(Nc, 0.0);
            zc_save_cc.resize(Nc, 0.0);
            xc_save_ac.resize(Nc, 0.0);
            yc_save_ac.resize(Nc, 0.0);
            zc_save_ac.resize(Nc, 0.0);
	        VL_cc.resize(200 * Nc, vector<int> (2, -1));
            VerletList_DPM_Disk_Append(Nc_old, r_cut_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac);
            VerletList_Disk_Append(Nc_old, r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc);
            VL_count_all[1] = VL_count_cc;
            VL_count_all[2] = VL_count_ac;

            xc_save_phi.resize(Nc, 0.0);
            yc_save_phi.resize(Nc, 0.0);
            zc_save_phi.resize(Nc, 0.0);
            phi.resize(Nc, 0.0);
	        VL_cc_phi.resize(400 * Nc, vector<int> (2, -1));
            is_in_periphery.resize(Nc, 0);
            is_dividing.resize(Nc, 0);
            is_not_dividing_ori.resize(Nc, 0);
            t_division.resize(Nc, -1.0);
            t_life.resize(Nc, 0.0);

            VerletList_Disk_Append(Nc_old, r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi);
            VerletList_Disk(r_cut_phi, xc_save_phi, yc_save_phi, zc_save_phi, VL_cc_phi, VL_count_phi, first_call);
            LocalPhi_Update(xc, yc, zc, r_shell, VL_cc_phi, VL_count_phi, is_dividing, is_not_dividing_ori, phi);
            Periphery_List_PF_Update(phi_c, phi, is_in_periphery, is_dividing, is_not_dividing_ori);
            for (n = 0; n < Nc; n++){
                if ((is_in_periphery[n] == 1) && (!is_dividing[n])){
                    t_division[n] = T_divide * gamma_dist(gen_division);
                    is_dividing[n] = 1;
                }
            }
            Nc_old = Nc;
        }

		VerletList_DPM(r_cut_aa, xa_save_aa, ya_save_aa, za_save_aa, VL_aa, VL_count_aa, first_call);
        VerletList_Disk(r_cut_cc, xc_save_cc, yc_save_cc, zc_save_cc, VL_cc, VL_count_cc, first_call);
        VerletList_DPM_Disk(r_cut_ac, xa_save_ac, ya_save_ac, za_save_ac, xc_save_ac, yc_save_ac, zc_save_ac, VL_ac, VL_count_ac, first_call);
        VL_count_all[0] = VL_count_aa;
        VL_count_all[1] = VL_count_cc;
        VL_count_all[2] = VL_count_ac;

        Contact_AC(VL_ac, VL_count_all[2], CC, reach_min, is_ll);
		for (n = 0; n < Na; n++){
            if (reach_min[n] == 0){
                if ((is_ll[n] == 1) || (CC[n] > CC_thresh)){
                    if (is_ll[n] == 0)
                        is_ll[n] = 1;
                    V0_temp = V0[n] - d_lip;
                    if (V0_temp > V0_min[n]){
                        V0[n] = V0_temp;
                        rsc = std::pow(V0_temp / V0_ori[n], 1.0 / 3.0);
                        R0[n] = rsc * R0_ori[n];
                        D0_intra[n] = rsc * D0_intra_ori[n];
                        D0_intra_sq[n] = rsc * rsc * D0_intra_sq_ori[n];
                        Ra[n] = rsc * Ra_ori[n];
                        rsc_all[n] = rsc;
                        for (nf = 0; nf < Nf; nf++)
                            A0[n][nf] = rsc * rsc * A0_ori[n][nf];
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
