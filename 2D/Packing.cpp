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

void BreastCancer_2D::AdipocytePacking()
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
            }
			Lx /= rsc;
            Ly /= rsc;
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

    double ya_cen;
    for (int na = 0; na < Na; na++){
        ya_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya_cen += ya[ns];
        ya_cen /= Ns[na];
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya[ns] -= std::floor(ya_cen / Ly) * Ly;
    }

	SaveConfig_Adipocyte(adip_name);
}

void BreastCancer_2D::AdipocytePacking_Compress(int P_start_c,
                                                int P_start_e,
                                                int P_end_c,
                                                int P_end_e,
                                                int P_num)
{
    double P;
    double ya_cen;

    std::ifstream posfile;
    posfile.open(adip_name);
    if (!posfile.good()){
        AdipocyteInitialization();

        double rsc = 1.01;

        while (true){
            FIRE_DPM_VL();
            P = Pressure_DPM();
            // printf("P: %.5e   L: %.5e\n", P, Lx);

            if (P < Pt){
                for (int ns = 0; ns < Ns_tot; ns++){
                    xa[ns] /= rsc;
                    ya[ns] /= rsc;
                }
                Lx /= rsc;
                Ly /= rsc;
            }
            else
                break;
        }

        FIRE_DPM_Enthalpy_VL(Pt);
        // double P = Pressure_DPM();
        // printf("P: %.5e   Lx: %.5e  Ly: %.5e\n", P, Lx, Ly);
        
        for (int na = 0; na < Na; na++){
            ya_cen = 0.0;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
                ya_cen += ya[ns];
            ya_cen /= Ns[na];
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
                ya[ns] -= std::floor(ya_cen / Ly) * Ly;
        }

        SaveConfig_Adipocyte(adip_name);
    }
    else
        LoadConfig_Adipocyte(adip_name);

    posfile.close();

    double P_start = (double)P_start_c * std::pow(10.0, P_start_e), P_end = (double)P_end_c * std::pow(10.0, P_end_e);
    double dP = (P_end - P_start) / (P_num - 1);

    std::ostringstream oss;
    string  adip_comp_name;

    for (int pn = 0; pn < P_num; pn++){
        oss << dir_run << "Pos_Adip_Comp_"
            << std::setfill('0') << std::setw(3) << pn << ".txt";
        adip_comp_name = oss.str();
        oss.str("");
        oss.clear();

        posfile.clear();
        posfile.open(adip_comp_name);
        if (posfile.good()){
            posfile.close();
            LoadConfig_Adipocyte_Comp(adip_comp_name);
            continue;
        }

        P = P_start + dP * pn;
        FIRE_DPM_Enthalpy_VL(P);

        for (int na = 0; na < Na; na++){
            ya_cen = 0.0;
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
                ya_cen += ya[ns];
            ya_cen /= Ns[na];
            for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
                ya[ns] -= std::floor(ya_cen / Ly) * Ly;
        }

	    SaveConfig_Adipocyte_Comp(adip_comp_name);
    }
}

void BreastCancer_2D::CancerPacking()
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
    while (true){
        FIRE_Disk_VL();
        P = Pressure_Disk();
        // printf("P: %.5e  Rc_circle: %.5e\n", P, Rc_circle);

        if (P > std::pow(10.0, -5))
            break;
        
        for (int n = 0; n < Nc; n++){
            xc[n] /= rsc;
        }
        Lx /= rsc;
    }

    if (std::abs(a_cut_c - 1.0) > std::pow(10.0, -7)){
        double  a_cut_c_set = a_cut_c;
        a_cut_c = 1.0;
        a_cut_h = 1.0;
        FIRE_Disk_Enthalpy_VL(Pt);
        P = Pressure_Disk();
        // printf("P: %.5e  Rc_circle: %.5e\n", P, Rc_circle);
        a_cut_c = a_cut_c_set;
        a_cut_h = 0.5 + 0.5 * a_cut_c;
        FIRE_Disk_Enthalpy_VL(Pt);
        // P = Pressure_Disk();
        // printf("P: %.5e  Rc_circle: %.5e\n", P, Rc_circle);
    }
    else {
        FIRE_Disk_Enthalpy_VL(Pt);
        P = Pressure_Disk();
        // printf("P: %.5e  Rc_circle: %.5e\n", P, Rc_circle);
    }

    for (int nc = 0; nc < Nc; nc++)
        yc[nc] -= std::floor(yc[nc] / Ly) * Ly;

	SaveConfig_Cancer(cancer_name);
}

void BreastCancer_2D::AllPacking()
{
    std::ifstream posfile(all_name.c_str());
    if (posfile.good()){
		posfile.close();
	    return; // configuration already generated
    }
    else
        posfile.close();

    LoadConfig_Adipocyte(adip_name); // Lx and Ly
    double  Lx_a = Lx;
	LoadConfig_Cancer(cancer_name); // Lx_c
    double  Lx_c = Lx;
    for (int na = 0; na < Na; na++)
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            xa[ns] += Lx_c;
    Lx = Lx_a + Lx_c;

    FIRE_DPM_Disk_Enthalpy_VL(Pt);
    double  P = Pressure_DPM_Disk();
    // printf("P: %.5e\n", P);

    double ya_cen;
    for (int na = 0; na < Na; na++){
        ya_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya_cen += ya[ns];
        ya_cen /= Ns[na];
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya[ns] -= std::floor(ya_cen / Ly) * Ly;
    }
    for (int nc = 0; nc < Nc; nc++)
        yc[nc] -= std::floor(yc[nc] / Ly) * Ly;

	SaveConfig_All(all_name);
}

void BreastCancer_2D::AllPacking_Smooth()
{
    std::ifstream posfile(all_name.c_str());
    if (posfile.good()){
		posfile.close();
	    return; // configuration already generated
    }
    else
        posfile.close();

    LoadConfig_Adipocyte(adip_name); // Lx and Ly
    double  Lx_a = Lx;
	LoadConfig_Cancer(cancer_name); // Lx_c
    double  Lx_c = Lx;
    for (int na = 0; na < Na; na++)
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            xa[ns] += Lx_c;
    Lx = Lx_a + Lx_c;

    FIRE_DPM_Disk_Smooth_Enthalpy_VL(Pt);
    double  P = Pressure_DPM_Disk();
    // printf("P: %.5e\n", P);

    double ya_cen;
    for (int na = 0; na < Na; na++){
        ya_cen = 0.0;
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya_cen += ya[ns];
        ya_cen /= Ns[na];
        for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
            ya[ns] -= std::floor(ya_cen / Ly) * Ly;
    }
    for (int nc = 0; nc < Nc; nc++)
        yc[nc] -= std::floor(yc[nc] / Ly) * Ly;

	SaveConfig_All(all_name);
}