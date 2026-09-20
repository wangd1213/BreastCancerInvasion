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
#include "Bumpy_3D.h"

using std::vector;

void BreastCancer_3D::Bumpy_Test(string dir_prefix)
{
    Bumpy_3D BumpyPacking(Na, Ns, delta, 1.0, 1.0, Pt_c, Pt_e, pid, dir_prefix);
	BumpyPacking.Packing_AboveJamming();
	BumpyPacking.ReturnBumpPos(xa, ya, za);
}

void BreastCancer_3D::Adipocyte_Check()
{
    LoadConfig_Adipocyte(adip_name);

	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL_aa(100 * Ns_tot, vector<int> (4, -1));
	int     		 		VL_counter = 0;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0), za_save(Ns_tot, 0.0);
	double 		r_cut = R0[0];
	for (int n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
	r_cut *= 2.0;
	
    VerletList_DPM(r_cut, xa_save, ya_save, za_save, VL_aa, VL_counter, first_call);
	Force_DPM_VL(VL_aa, VL_counter);
    // Force_DPM();

    printf("Max Force: %.5e\n", MaxForce_DPM());
}

void BreastCancer_3D::AdipocytePacking_FromRigid(string dir_pre)
{
    std::ifstream posfile(adip_name);
    if (posfile.good()){
        posfile.close();
        return;
    }
    else
        posfile.close();

    AdipocyteInitialization_FromRigid(dir_pre);

    FIRE_DPM_Enthalpy_VL(Pt);
	// double P = Pressure_DPM();
	// printf("P: %.5e   Lx: %.5e  Ly: %.5e\n", P, Lx, Ly);

    if (Kb_std > 0.000001){
        /*
        // Assign Kb now based on the x location
        vector<double> 	xa_cen(Na, 0.0);
        for (na = 0; na < Na; na++){
            for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
                xa_cen[na] += pos_new[ns];
            xa_cen[na] /= Ns[na];
        }
        vector<size_t> 	idx_sort = sort_indices(xa_cen);
        */

        std::mt19937 gen(pid); //Standard mersenne_twister_engine seeded with rd()
        std::normal_distribution<> dis_kb(0.0, 1.0);

        double      Kb_min = 0.25 * Kb;
        double      Kb_max = 1.75 * Kb;
        double      Kb_na;
        // vector<double> 	Kb_list(Na, 0.0);
        for (int na = 0; na < Na; na++){
            while (true){
                Kb_na = (1.0 + Kb_std * dis_kb(gen)) * Kb;
                if ((Kb_na >= Kb_min) && (Kb_na <= Kb_max))
                    break;
            }
            // Kb_list[na] = Kb_na;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
                Kb_all[ns] = Kb_na;
        }
        
        /*
        std::sort(Kb_list.begin(), Kb_list.end());
        for (int n = 0; n < Na; n++){
            na = idx_sort[n];
            for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
                Kb_all[ns] = Kb_list[n];
        }
        */

        // dt_fire = 0.001;
        FIRE_DPM_Enthalpy_VL(Pt);
	    // P = Pressure_DPM_RigidSquare();
    }

    double xa_cen, za_cen;
    for (int na = 0; na < Na; na++){
        xa_cen = 0.0;
        za_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa_cen += xa[ns];
            za_cen += za[ns];
        }
        xa_cen /= Ns;
        za_cen /= Ns;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa[ns] -= std::floor(xa_cen / Lx) * Lx;
            za[ns] -= std::floor(za_cen / Lz) * Lz;
        }
    }

	SaveConfig_Adipocyte(adip_name);
}

void BreastCancer_3D::AdipocytePacking_Compress(string dir_pre,
                                                int P_start_c,
                                                int P_start_e,
                                                int P_end_c,
                                                int P_end_e,
                                                int P_num)
{
    std::ifstream posfile(adip_name);
    if (!posfile.good())
        AdipocytePacking_FromRigid(dir_pre);
    posfile.close();

    LoadConfig_Adipocyte(adip_name); // load Lx and Ly, Lx will be changed in CancerInitialization()

    double P, P_start = (double)P_start_c * std::pow(10.0, P_start_e), P_end = (double)P_end_c * std::pow(10.0, P_end_e);
    double dP = (P_end - P_start) / (P_num - 1);

    std::ostringstream oss;
    string  adip_comp_name;

    for (int pn = 0; pn < P_num; pn++){
        P = P_start + dP * pn;
        FIRE_DPM_Enthalpy_VL(P);

        double xa_cen, za_cen;
        for (int na = 0; na < Na; na++){
            xa_cen = 0.0;
            za_cen = 0.0;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                xa_cen += xa[ns];
                za_cen += za[ns];
            }
            xa_cen /= Ns;
            za_cen /= Ns;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                xa[ns] -= std::floor(xa_cen / Lx) * Lx;
                za[ns] -= std::floor(za_cen / Lz) * Lz;
            }
        }

        oss << dir_run << "Pos_Adip_Comp_"
            << std::setfill('0') << std::setw(3) << pn << ".txt";
        adip_comp_name = oss.str();
        oss.str("");
        oss.clear();
	    SaveConfig_Adipocyte_Comp(adip_comp_name);
    }
}

void BreastCancer_3D::AdipocytePacking()
{
    std::ifstream posfile(adip_name);
    if (posfile.good()){
        posfile.close();
        return;
    }
    else
        posfile.close();

    AdipocyteInitialization();

    double P;
    double rsc = 1.01;

    while (true){
        FIRE_DPM_VL();
        P = Pressure_DPM();
		// printf("P: %.5e   L: %.5e\n", P, Lx);

		if (P < Pt){
	        for (int ns = 0; ns < Ns_tot; ns++){
	        	xa[ns] /= rsc;
                ya[ns] /= rsc;
                za[ns] /= rsc;
            }
			Lx /= rsc;
            Ly /= rsc;
            Lz /= rsc;
	    }
	    else
	        break;
	}

    FIRE_DPM_Enthalpy_VL(Pt);
	P = Pressure_DPM();
	// printf("P: %.5e   Lx: %.5e  Ly: %.5e\n", P, Lx, Ly);

    if (Kb_std > 0.000001){
        /*
        // Assign Kb now based on the x location
        vector<double> 	xa_cen(Na, 0.0);
        for (na = 0; na < Na; na++){
            for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
                xa_cen[na] += pos_new[ns];
            xa_cen[na] /= Ns[na];
        }
        vector<size_t> 	idx_sort = sort_indices(xa_cen);
        */

        std::mt19937 gen(pid); //Standard mersenne_twister_engine seeded with rd()
        std::normal_distribution<> dis_kb(0.0, 1.0);

        double      Kb_min = 0.25 * Kb;
        double      Kb_max = 1.75 * Kb;
        double      Kb_na;
        // vector<double> 	Kb_list(Na, 0.0);
        for (int na = 0; na < Na; na++){
            while (true){
                Kb_na = (1.0 + Kb_std * dis_kb(gen)) * Kb;
                if ((Kb_na >= Kb_min) && (Kb_na <= Kb_max))
                    break;
            }
            // Kb_list[na] = Kb_na;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
                Kb_all[ns] = Kb_na;
        }
        
        /*
        std::sort(Kb_list.begin(), Kb_list.end());
        for (int n = 0; n < Na; n++){
            na = idx_sort[n];
            for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
                Kb_all[ns] = Kb_list[n];
        }
        */

        // dt_fire = 0.001;
        FIRE_DPM_Enthalpy_VL(Pt);
	    // P = Pressure_DPM_RigidSquare();
    }

    double xa_cen, za_cen;
    for (int na = 0; na < Na; na++){
        xa_cen = 0.0;
        za_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa_cen += xa[ns];
            za_cen += za[ns];
        }
        xa_cen /= Ns;
        za_cen /= Ns;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa[ns] -= std::floor(xa_cen / Lx) * Lx;
            za[ns] -= std::floor(za_cen / Lz) * Lz;
        }
    }

	SaveConfig_Adipocyte(adip_name);
}

void BreastCancer_3D::CancerPacking()
{
    std::ifstream posfile(cancer_name.c_str());
    if (posfile.good()){
		posfile.close();
	    return; // configuration already generated
    }
    else
        posfile.close();

    LoadConfig_Adipocyte(adip_name); // load Lx and Ly, Lx will be changed in CancerInitialization()
	CancerInitialization();
    
    double P;
    double rsc = 1.01;

    double a_cut_c_set = 1.0;

    if (std::abs(a_cut_c - 1.0) > std::pow(10.0, -7)){
        a_cut_c_set = a_cut_c;
        a_cut_c = 1.0;
        a_cut_h = 1.0;
    }

    while (true){
        FIRE_Disk_VL();
        P = Pressure_Disk();
        // printf("P: %.5e  Ly: %.5e\n", P, Ly);

        if (P > std::pow(10.0, -5))
            break;
        
        for (int nc = 0; nc < Nc; nc++)
            yc[nc] /= rsc;
        Ly /= rsc;
    }

    FIRE_Disk_Enthalpy_VL(Pt);

    if (std::abs(a_cut_c_set - 1.0) > std::pow(10.0, -7)){
        a_cut_c = a_cut_c_set;
        a_cut_h = 0.5 + 0.5 * a_cut_c;
        FIRE_Disk_Enthalpy_VL(Pt);
        // P = Pressure_Disk();
        // printf("P: %.5e  Rc_circle: %.5e\n", P, Rc_circle);
    }

    for (int nc = 0; nc < Nc; nc++){
        xc[nc] -= std::floor(xc[nc] / Lx) * Lx;
        zc[nc] -= std::floor(zc[nc] / Lz) * Lz;
    }

	SaveConfig_Cancer(cancer_name);
}

void BreastCancer_3D::AllPacking()
{
    std::ifstream posfile(all_name.c_str());
    if (posfile.good()){
		posfile.close();
	    return; // configuration already generated
    }
    else
        posfile.close();
    
    LoadConfig_Adipocyte(adip_name); // Lx, Ly, and Lz
    double  Ly_a = Ly;
	LoadConfig_Cancer(cancer_name); // Ly
    double  Ly_c = Ly;
    for (int na = 0; na < Na; na++)
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya[ns] += Ly_c;
    Ly = Ly_a + Ly_c;

    FIRE_DPM_Disk_Enthalpy_VL(Pt);
    // double  P = Pressure_DPM_Disk();
    // printf("P: %.5e\n", P);

    double  xa_cen, za_cen;
    for (int na = 0; na < Na; na++){
        xa_cen = 0.0;
        za_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa_cen += xa[ns];
            za_cen += za[ns];
        }
        xa_cen /= Ns;
        za_cen /= Ns;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa[ns] -= std::floor(xa_cen / Lx) * Lx;
            za[ns] -= std::floor(za_cen / Lz) * Lz;
        }
    }
    for (int nc = 0; nc < Nc; nc++){
        xc[nc] -= std::floor(xc[nc] / Lx) * Lx;
        zc[nc] -= std::floor(zc[nc] / Lz) * Lz;
    }

	SaveConfig_All(all_name);
}

void BreastCancer_3D::AllPacking_Immediate()
{
    std::ostringstream oss;
    string  all_immediate_name;

    oss << dir_run << "Pos_All_Immediate_" << Kr_str_c << "_" << Kr_str_ac << ".txt";
    all_immediate_name = oss.str();
    oss.str("");
    oss.clear();

    std::ifstream posfile(all_immediate_name.c_str());
    if (posfile.good()){
		posfile.close();
	    return; // configuration already generated
    }
    else
        posfile.close();

    LoadConfig_Adipocyte(adip_name); // Lx, Ly, and Lz
    double  Ly_a = Ly;
	LoadConfig_Cancer(cancer_name); // Ly_c
    double  Ly_c = Ly;
    for (int na = 0; na < Na; na++)
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya[ns] += Ly_c;
    Ly = Ly_a + Ly_c;

    // FIRE_DPM_Disk_VL();
    // double  P = Pressure_DPM_Disk();
    // printf("P: %.5e\n", P);

    for (int nc = 0; nc < Nc; nc++)
        yc[nc] = Ly - yc[nc];
    for (int ns = 0; ns < Ns_tot; ns++)
        ya[ns] = Ly - ya[ns];
    FIRE_DPM_Disk_Enthalpy_VL(Pt);
    for (int nc = 0; nc < Nc; nc++)
        yc[nc] = Ly - yc[nc];
    for (int ns = 0; ns < Ns_tot; ns++)
        ya[ns] = Ly - ya[ns];
    // P = Pressure_DPM_Disk();
    // printf("P: %.5e\n", P);

    double  xa_cen, za_cen;
    for (int na = 0; na < Na; na++){
        xa_cen = 0.0;
        za_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa_cen += xa[ns];
            za_cen += za[ns];
        }
        xa_cen /= Ns;
        za_cen /= Ns;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa[ns] -= std::floor(xa_cen / Lx) * Lx;
            za[ns] -= std::floor(za_cen / Lz) * Lz;
        }
    }
    for (int nc = 0; nc < Nc; nc++){
        xc[nc] -= std::floor(xc[nc] / Lx) * Lx;
        zc[nc] -= std::floor(zc[nc] / Lz) * Lz;
    }

	SaveConfig_All(all_immediate_name);
}

void BreastCancer_3D::AllPacking_ThreeSteps()
{
    std::ostringstream oss;
    string  all_immediate_name;

    oss << dir_run << "Pos_All_AfterFix_" << Kr_str_c << "_" << Kr_str_ac << ".txt";
    all_immediate_name = oss.str();
    oss.str("");
    oss.clear();

    std::ifstream posfile(all_immediate_name.c_str());
    if (posfile.good()){
		posfile.close();
	    return; // configuration already generated
    }
    else
        posfile.close();

    LoadConfig_Adipocyte(adip_name); // Lx and Ly
    double  Ly_a = Ly;
    LoadConfig_Cancer(cancer_name); // Lx_c
    double  Ly_c = Ly;
    for (int na = 0; na < Na; na++)
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya[ns] += Ly_c;
    Ly = Ly_a + Ly_c;

    FIRE_DPM_FixDisk_Enthalpy_VL(Pt);
    // double  P = Pressure_DPM_Disk();
    // printf("P: %.5e\n", P);

    double  xa_cen, za_cen;
    for (int na = 0; na < Na; na++){
        xa_cen = 0.0;
        za_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa_cen += xa[ns];
            za_cen += za[ns];
        }
        xa_cen /= Ns;
        za_cen /= Ns;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa[ns] -= std::floor(xa_cen / Lx) * Lx;
            za[ns] -= std::floor(za_cen / Lz) * Lz;
        }
    }

    oss << dir_run << "Pos_All_Fix_Adip_" << Kr_str_c << "_" << Kr_str_ac << ".txt";
    all_immediate_name = oss.str();
    oss.str("");
    oss.clear();

    SaveConfig_All(all_immediate_name);

    for (int nc = 0; nc < Nc; nc++)
        yc[nc] = Ly - yc[nc];
    for (int ns = 0; ns < Ns_tot; ns++)
        ya[ns] = Ly - ya[ns];
    FIRE_Disk_FixDPM_Enthalpy_VL(Pt);
    // double  P = Pressure_DPM_Disk();
    // printf("P: %.5e\n", P);

    for (int nc = 0; nc < Nc; nc++){
        xc[nc] -= std::floor(xc[nc] / Lx) * Lx;
        zc[nc] -= std::floor(zc[nc] / Lz) * Lz;
    }

    oss << dir_run << "Pos_All_Fix_Cancer_" << Kr_str_c << "_" << Kr_str_ac << ".txt";
    all_immediate_name = oss.str();
    oss.str("");
    oss.clear();

    SaveConfig_All(all_immediate_name);

    FIRE_DPM_Disk_VL();
    FIRE_DPM_Disk_Enthalpy_VL(Pt);
    for (int nc = 0; nc < Nc; nc++)
        yc[nc] = Ly - yc[nc];
    for (int ns = 0; ns < Ns_tot; ns++)
        ya[ns] = Ly - ya[ns];
    // P = Pressure_DPM_Disk();
    // printf("P: %.5e\n", P);

    for (int na = 0; na < Na; na++){
        xa_cen = 0.0;
        za_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa_cen += xa[ns];
            za_cen += za[ns];
        }
        xa_cen /= Ns;
        za_cen /= Ns;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa[ns] -= std::floor(xa_cen / Lx) * Lx;
            za[ns] -= std::floor(za_cen / Lz) * Lz;
        }
    }
    for (int nc = 0; nc < Nc; nc++){
        xc[nc] -= std::floor(xc[nc] / Lx) * Lx;
        zc[nc] -= std::floor(zc[nc] / Lz) * Lz;
    }

    oss << dir_run << "Pos_All_AfterFix_" << Kr_str_c << "_" << Kr_str_ac << ".txt";
    all_immediate_name = oss.str();
    oss.str("");
    oss.clear();

	SaveConfig_All(all_immediate_name);
}

void BreastCancer_3D::AllPacking_Compress()
{
    std::ostringstream oss;
    string  all_immediate_name;

    LoadConfig_Adipocyte(adip_name); // Lx and Ly
    double  Ly_a = Ly;
    LoadConfig_Cancer(cancer_name); // Lx_c
    double  Ly_c = Ly;
    for (int na = 0; na < Na; na++)
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya[ns] += Ly_c;
    Ly = Ly_a + Ly_c;

    for (int nc = 0; nc < Nc; nc++)
        yc[nc] = Ly - yc[nc];
    for (int ns = 0; ns < Ns_tot; ns++)
        ya[ns] = Ly - ya[ns];

    FIRE_DPM_Disk_VL();
    FIRE_DPM_Disk_Enthalpy_VL(Pt);
    for (int nc = 0; nc < Nc; nc++)
        yc[nc] = Ly - yc[nc];
    for (int ns = 0; ns < Ns_tot; ns++)
        ya[ns] = Ly - ya[ns];
    // P = Pressure_DPM_Disk();
    // printf("P: %.5e\n", P);

    double xa_cen, za_cen;
    for (int na = 0; na < Na; na++){
        xa_cen = 0.0;
        za_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa_cen += xa[ns];
            za_cen += za[ns];
        }
        xa_cen /= Ns;
        za_cen /= Ns;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            xa[ns] -= std::floor(xa_cen / Lx) * Lx;
            za[ns] -= std::floor(za_cen / Lz) * Lz;
        }
    }
    for (int nc = 0; nc < Nc; nc++){
        xc[nc] -= std::floor(xc[nc] / Lx) * Lx;
        zc[nc] -= std::floor(zc[nc] / Lz) * Lz;
    }

    SaveConfig_All(all_name);

    double  Pt_ori = Pt;
    int     P_num = (int)std::ceil(0.1 / Pt_ori);
    for (int pn = 1; pn < P_num; pn++){
        Pt += Pt_ori;
        for (int nc = 0; nc < Nc; nc++)
            yc[nc] = Ly - yc[nc];
        for (int ns = 0; ns < Ns_tot; ns++)
            ya[ns] = Ly - ya[ns];
        FIRE_DPM_Disk_Enthalpy_VL(Pt);
        for (int nc = 0; nc < Nc; nc++)
            yc[nc] = Ly - yc[nc];
        for (int ns = 0; ns < Ns_tot; ns++)
            ya[ns] = Ly - ya[ns];

        for (int na = 0; na < Na; na++){
            xa_cen = 0.0;
            za_cen = 0.0;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                xa_cen += xa[ns];
                za_cen += za[ns];
            }
            xa_cen /= Ns;
            za_cen /= Ns;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                xa[ns] -= std::floor(xa_cen / Lx) * Lx;
                za[ns] -= std::floor(za_cen / Lz) * Lz;
            }
        }
        for (int nc = 0; nc < Nc; nc++){
            xc[nc] -= std::floor(xc[nc] / Lx) * Lx;
            zc[nc] -= std::floor(zc[nc] / Lz) * Lz;
        }

        oss << dir_run << "Pos_All_" << Kr_str_c << "_" << Kr_str_ac
            << "_" << std::setfill('0') << std::setw(3) << pn << ".txt";
        all_immediate_name = oss.str();
        oss.str("");
        oss.clear();

        SaveConfig_All(all_immediate_name);
    }
}

void BreastCancer_3D::AllPacking_DeCompress()
{
    std::ostringstream oss;
    string  all_immediate_name;
    double  xa_cen, za_cen;

    std::ifstream Allpos;
    Allpos.open(all_name.c_str());
    if (!Allpos.good()){
        LoadConfig_Adipocyte(adip_name); // Lx and Ly
        double  Ly_a = Ly;
        LoadConfig_Cancer(cancer_name); // Ly_c
        double  Ly_c = Ly;
        for (int na = 0; na < Na; na++)
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
                ya[ns] += Ly_c;
        Ly = Ly_a + Ly_c;
        
        /*
        for (int nc = 0; nc < Nc; nc++)
            yc[nc] = Ly - yc[nc];
        for (int ns = 0; ns < Ns_tot; ns++)
            ya[ns] = Ly - ya[ns];

        // FIRE_DPM_Disk_VL();
        FIRE_DPM_Disk_Enthalpy_VL(Pt);
        for (int nc = 0; nc < Nc; nc++)
            yc[nc] = Ly - yc[nc];
        for (int ns = 0; ns < Ns_tot; ns++)
            ya[ns] = Ly - ya[ns];
        */
        FIRE_DPM_Disk_Enthalpy_VL(Pt);
        // P = Pressure_DPM_Disk();
        // printf("P: %.5e\n", P);

        for (int na = 0; na < Na; na++){
            xa_cen = 0.0;
            za_cen = 0.0;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                xa_cen += xa[ns];
                za_cen += za[ns];
            }
            xa_cen /= Ns;
            za_cen /= Ns;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                xa[ns] -= std::floor(xa_cen / Lx) * Lx;
                za[ns] -= std::floor(za_cen / Lz) * Lz;
            }
        }
        for (int nc = 0; nc < Nc; nc++){
            xc[nc] -= std::floor(xc[nc] / Lx) * Lx;
            zc[nc] -= std::floor(zc[nc] / Lz) * Lz;
        }

        SaveConfig_All(all_name);
    }
    else {
        LoadConfig_All(all_name);
    }
    Allpos.close();

    double  dPt = std::pow(10.0, -3);
    int     P_num = (int)std::floor(Pt / dPt);

    for (int pn = 1; pn < P_num; pn++){
        Pt -= dPt;

        oss << dir_run << "Pos_All_DeComp_" << Kr_str_c << "_" << Kr_str_ac
            << "_" << std::setfill('0') << std::setw(3) << pn << ".txt";
        all_immediate_name = oss.str();
        oss.str("");
        oss.clear();

        Allpos.open(all_immediate_name);
        if (Allpos.good()){
            Allpos.close();
            LoadConfig_All(all_immediate_name);
            continue;
        }
        else {
            Allpos.close();

            /*
            for (int nc = 0; nc < Nc; nc++)
                yc[nc] = Ly - yc[nc];
            for (int ns = 0; ns < Ns_tot; ns++)
                ya[ns] = Ly - ya[ns];
            FIRE_DPM_Disk_Enthalpy_VL(Pt);
            for (int nc = 0; nc < Nc; nc++)
                yc[nc] = Ly - yc[nc];
            for (int ns = 0; ns < Ns_tot; ns++)
                ya[ns] = Ly - ya[ns];
            */
            FIRE_DPM_Disk_Enthalpy_VL(Pt);

            for (int na = 0; na < Na; na++){
                xa_cen = 0.0;
                za_cen = 0.0;
                for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                    xa_cen += xa[ns];
                    za_cen += za[ns];
                }
                xa_cen /= Ns;
                za_cen /= Ns;
                for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                    xa[ns] -= std::floor(xa_cen / Lx) * Lx;
                    za[ns] -= std::floor(za_cen / Lz) * Lz;
                }
            }
            for (int nc = 0; nc < Nc; nc++){
                xc[nc] -= std::floor(xc[nc] / Lx) * Lx;
                zc[nc] -= std::floor(zc[nc] / Lz) * Lz;
            }

            SaveConfig_All(all_immediate_name);
        }
    }
}