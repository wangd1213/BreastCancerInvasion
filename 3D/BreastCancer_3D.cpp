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

// Constructor implementation
BreastCancer_3D::BreastCancer_3D(int Na_set,
									int Ns_set,
									int Nc_set,
									double delta_set,
									int dAlf_c,
									int dAlf_e,
									double KV_set,
									double KA_set,
									double Kaa_set,
									int Kb_c_set,
									int Kb_e_set,
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
									int seed,
									string dir_pre)
    : Na(Na_set), Ns(Ns_set), Nc(Nc_set),
	  delta(delta_set), KV(KV_set), KA(KA_set), Kaa(Kaa_set),
	  Kb_std(Kb_std_set), Kw_a(Kw_a_set),
	  Kpin(Kpin_set), Kcc(Kcc_set), a_cut_c(a_cut_c_set),
	  a_cut_ac(a_cut_ac_set), Kac(Kac_set), Kw_c(Kw_c_set),
	  pid(seed), Pt_c(Pt_c_set), Pt_e(Pt_e_set)
{
    // Resize vectors to accommodate n disks
	Nf = 2 * Ns - 4;
	Ne = 3 * Ns - 6;
	Nnn = (Ns - 7) * Ns / 2 + 6;
    Ns_tot = Ns * Na;
    xa.assign(Ns_tot, 0.0);
    ya.assign(Ns_tot, 0.0);
	za.assign(Ns_tot, 0.0);
	xa_pin.assign(Ns_tot, 0.0);
    ya_pin.assign(Ns_tot, 0.0);
	za_pin.assign(Ns_tot, 0.0);
    xa_cen_pin.assign(Na, 0.0);
    ya_cen_pin.assign(Na, 0.0);
	za_cen_pin.assign(Na, 0.0);
    xc.assign(Nc, 0.0);
    yc.assign(Nc, 0.0);
	zc.assign(Nc, 0.0);
	Ra.assign(Na, 0.0);
	R0.assign(Na, 0.0);
	D0_intra.assign(Na, 0.0);
	D0_intra_sq.assign(Na, 0.0);
	V0.assign(Na, 0.0);
	A0.assign(Na, vector<double> (Nf, 0.0));
	theta0.assign(Ne, 0.0);
	idx_list.assign(Na + 1, 0);
	edgelist.assign(Ne, vector<int> (6, 0));
    f_unit.assign(Nf, vector<int> (3, 0));
    nnlist.assign(Nnn, vector<int> (2, 0));
	Rc.assign(Nc, 0.0);
	Dc.assign(Nc, 0.0);
	Rc_sq.assign(Nc, 0.0);
	Dc_cube.assign(Nc, 0.0);
	Fx_a.assign(Ns_tot, 0.0);
	Fy_a.assign(Ns_tot, 0.0);
	Fz_a.assign(Ns_tot, 0.0);
	Fx_c.assign(Nc, 0.0);
	Fy_c.assign(Nc, 0.0);
	Fz_c.assign(Nc, 0.0);
	Kb_all.assign(Ns_tot, 0.0);

	std::ostringstream oss;
	oss << dir_pre << "Polyhedron/Ns_" << Ns << "/DVA.txt";
    dva_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_pre << "Polyhedron/Ns_" << Ns << "/F_unit.txt";
    f_unit_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_pre << "Polyhedron/Ns_" << Ns << "/EdgeList.txt";
    edgelist_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_pre << "Polyhedron/Ns_" << Ns << "/NonNeighborList.txt";
    nnn_list_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_pre << "Polyhedron/Ns_" << Ns << "/Theta0.txt";
    theta_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_pre << "Polyhedron/Ns_" << Ns << "/XYZ_unit.txt";
    unit_pos_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_pre << "Polyhedron/Ns_" << Ns << "/RigidUnit.txt";
    rigid_unit_name = oss.str();
    oss.str("");
    oss.clear();
	
	oss << "KV_" << std::fixed << std::setprecision(3) << KV
		<< "_KA_" << std::fixed << std::setprecision(3) << KA
		<< "_Kb_" << Kb_c_set << "E" << Kb_e_set
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
		<< "_Ns_" << std::setfill('0') << std::setw(3) << Ns
		<< "_Nc_" << std::setfill('0') << std::setw(4) << Nc
		<< "/" << K_adip_str << "/Dr_" << std::fixed << std::setprecision(2) << delta
        << "/" << std::setfill('0') << std::setw(5) << pid
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

	oss << dir_pre << "Rigid/Na_" << std::setfill('0') << std::setw(3) << Na
		<< "_Ns_" << std::setfill('0') << std::setw(3) << Ns
        << "/Dr_" << std::fixed << std::setprecision(2) << delta
		<< "/Pos_Rigid_" << std::setfill('0') << std::setw(5) << pid << ".txt";
    rigid_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_run << "Pos_Adip.txt";
    adip_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_run << "Pos_Cancer_" << Kr_str_c << ".txt";
    cancer_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_run << "Pos_All_" << Kr_str_c << "_" << Kr_str_ac << ".txt";
	all_name = oss.str();
	oss.str("");
	oss.clear();

	Alf = 1.0 + (double)dAlf_c * std::pow(10.0, dAlf_e);

	for (int na = 1; na <= Na; na++)
		idx_list[na] = na * Ns;

	SetParticleParameters();

    Kaa *= 4.0 * R0[0] * R0[0];
    Kcc *= 4.0 * Rs * Rs;
    Kac *= (R0[0] + Rs) * (R0[0] + Rs);
    Kw_a = Kaa / 4.0;
    Kw_c = Kcc / 4.0;
    Kcc_attr = Kcc * Kcc_attr_ratio;
    a_cut_h = 0.5 + 0.5 * a_cut_c;
	Kac_attr = Kac * Kac_attr_ratio;
	a_cut_ac_h = 0.5 + 0.5 * a_cut_ac;
	Pt = (double)Pt_c * std::pow(10.0, Pt_e) * KV;
	Kb = (double)Kb_c_set * std::pow(10.0, Kb_e_set);
}
