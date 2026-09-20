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

// Constructor implementation
BreastCancer_2D::BreastCancer_2D(int Na_set,
									int Nss_set,
									int Nsb_set,
									int Nc_set,
									double delta_set,
									int dAlf_c,
									int dAlf_e,
									double Kp_set,
									double KA_set,
									double Kaa_set,
									double Kb_set,
									double Kb_std_set,
									double Kw_a_set,
									double Kpin_set,
									double Kcc_set,
									double a_cut_c_set,
									double Kcc_attr_ratio,
									double a_cut_ac_set,
									double Kac_attr_ratio,
									double Kac_set,
									double Kw_c_set,
									int Pt_c_set,
									int Pt_e_set,
									int is_smooth_set,
									int seed,
									string dir_pre)
    : Na(Na_set), Nss(Nss_set), Nsb(Nsb_set), Nc(Nc_set),
	  delta(delta_set), Kp(Kp_set), KA(KA_set), Kaa(Kaa_set),
	  Kb(Kb_set), Kb_std(Kb_std_set), Kw_a(Kw_a_set),
	  Kpin(Kpin_set), Kcc(Kcc_set), a_cut_c(a_cut_c_set),
	  a_cut_ac(a_cut_ac_set), Kac(Kac_set), Kw_c(Kw_c_set),
	  is_smooth(is_smooth_set), pid(seed),
	  Pt_c(Pt_c_set), Pt_e(Pt_e_set)
{
    // Resize vectors to accommodate n disks
    Ns_tot = (Nss + Nsb) * Na / 2;
    xa.assign(Ns_tot, 0.0);
    ya.assign(Ns_tot, 0.0);
	xa_pin.assign(Ns_tot, 0.0);
    ya_pin.assign(Ns_tot, 0.0);
    xa_cen_pin.assign(Na, 0.0);
    ya_cen_pin.assign(Na, 0.0);
    xc.assign(Nc, 0.0);
    yc.assign(Nc, 0.0);
	Da.assign(Na, 0.0);
	D0.assign(Na, 0.0);
	R0.assign(Na, 0.0);
	L0.assign(Na, 0.0);
	A0.assign(Na, 0.0);
	alpha_0.assign(Na, 0.0);
	idx_list.assign(Na + 1, 0);
	ift.assign(Ns_tot, 0);
	jft.assign(Ns_tot, 0);
	Ns.assign(Na, 0);
	Rc.assign(Nc, 0.0);
	Dc.assign(Nc, 0.0);
	Rc_sq.assign(Nc, 0.0);
	lx.assign(Ns_tot, 0.0);
	ly.assign(Ns_tot, 0.0);
	lk.assign(Ns_tot, 0.0);
	Fx_a.assign(Ns_tot, 0.0);
	Fy_a.assign(Ns_tot, 0.0);
	Fx_c.assign(Nc, 0.0);
	Fy_c.assign(Nc, 0.0);
	Kb_all.assign(Na, 0.0);

	Alf = 1.0 + (double)dAlf_c * std::pow(10.0, dAlf_e);

    Kcc_attr = Kcc * Kcc_attr_ratio;
    a_cut_h = 0.5 + 0.5 * a_cut_c;
	Kac_attr = Kac * Kac_attr_ratio;
	a_cut_ac_h = 0.5 + 0.5 * a_cut_ac;
	Pt = (double)Pt_c * std::pow(10.0, Pt_e) * KA;

	for (int na = 0; na < Na; na++){
		Kb_all[na] = Kb;
    	if (na < Na / 2) {
    		Ns[na] = Nsb;
    		Da[na] = D0_a / std::sin(PI / Nsb); // D0_a set in the header
    		D0[na] = delta * D0_a;
			R0[na] = 0.5 * D0[na];
    		L0[na] = D0_a;
    		A0[na] = Nsb * D0_a * D0_a / 4.0 / std::tan(PI / Nsb) / Alf;
			alpha_0[na] = 2.0 * PI / Nsb;
			idx_list[na + 1] = idx_list[na] + Nsb;
			
			for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				if (ns == idx_list[na]){
					ift[ns] = ns + 1;
					jft[ns] = idx_list[na + 1] - 1;
				}
				else if (ns == idx_list[na + 1] - 1){
					ift[ns] = idx_list[na];
					jft[ns] = ns - 1;
				}
				else{
					ift[ns] = ns + 1;
					jft[ns] = ns - 1;
				}
			}
    	}
    	else {
    		Ns[na] = Nss;
    		Da[na] = D0_a / std::sin(PI / Nss);
    		D0[na] = delta * D0_a;
			R0[na] = 0.5 * D0[na];
    		L0[na] = D0_a;
    		A0[na] = Nss * D0_a * D0_a / 4.0 / std::tan(PI / Nss) / Alf;
			alpha_0[na] = 2.0 * PI / Nss;
			idx_list[na + 1] = idx_list[na] + Nss;
			for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				if (ns == idx_list[na]){
					ift[ns] = ns + 1;
					jft[ns] = idx_list[na + 1] - 1;
				}
				else if (ns == idx_list[na + 1] - 1){
					ift[ns] = idx_list[na];
					jft[ns] = ns - 1;
				}
				else{
					ift[ns] = ns + 1;
					jft[ns] = ns - 1;
				}
			}
    	}
    }

	std::ostringstream oss;
	oss << "KA_" << std::fixed << std::setprecision(3) << KA
		<< "_Kp_" << std::fixed << std::setprecision(3) << Kp
		<< "_Kb_" << std::fixed << std::setprecision(5) << Kb
		<< "_" << std::fixed << std::setprecision(4) << Kb_std;
    K_adip_str = oss.str();
    oss.str("");
    oss.clear();

	oss << "P_" << Pt_c << "E" << Pt_e;
    P_str = oss.str();
    oss.str("");
    oss.clear();

	oss << "Kcc_attr_" << std::fixed << std::setprecision(2) << Kcc_attr_ratio
		<< "_acc_" << std::fixed << std::setprecision(2) << a_cut_c;
    Kr_str_c = oss.str();
    oss.str("");
    oss.clear();

	oss << "Kac_attr_" << std::fixed << std::setprecision(2) << Kac_attr_ratio
		<< "_aac_" << std::fixed << std::setprecision(2) << a_cut_ac;
    Kr_str_ac = oss.str();
    oss.str("");
    oss.clear();

	oss << "dAlf_" << dAlf_c << "E" << dAlf_e;
    string dAlf_str = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_pre << "Na_" << std::setfill('0') << std::setw(3) << Na
		<< "_Ns_" << std::setfill('0') << std::setw(2) << Nss
		<< "_Nc_" << std::setfill('0') << std::setw(4) << Nc
        << "/"  << K_adip_str << "/" << std::setfill('0') << std::setw(5) << pid
		<< "/" << dAlf_str << "/" << P_str << "/";
    dir_run = oss.str();
    oss.str("");
    oss.clear();

	string cmd;
	oss << "mkdir -p " << dir_run;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_run << "Pos_Adip.txt";
    adip_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_run << "Pos_Cancer_" << Kr_str_c << ".txt";
    cancer_name = oss.str();
    oss.str("");
    oss.clear();

	if (is_smooth == 0){
		oss << dir_run << "Pos_All_" << Kr_str_c << "_" << Kr_str_ac << ".txt";
		all_name = oss.str();
		oss.str("");
		oss.clear();
	}
	else {
		oss << dir_run << "Pos_All_Smooth_" << Kr_str_c << "_" << Kr_str_ac << ".txt";
		all_name = oss.str();
		oss.str("");
		oss.clear();
	}
}
