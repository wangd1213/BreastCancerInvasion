#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string>
#include <time.h>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <vector>

#include "BreastCancer_2D.h"

using std::vector;

int main(int argc, char *argv[])
{
    /*
    int 	Na = 64;
    int 	Nss = 20;
    int 	Nsb = 28;
    int 	Nc = 400;
    double 	delta = 2.0;
    double 	Kb = 0.001;
    double 	Kb_std = 0.0;
    double 	a_cut_c = 1.0;
	double  Kcc_attr_ratio = 1.0;
    int    	Pt_c = 1;
    int    	Pt_e = -2;
    double 	Kpin = 0.05;
    int 	pid = 2;
    */
    // double 	Kp = 0.5;
    // double 	KA = 1.0;
    double 	Kaa = 0.1;
    double 	Kw_a = Kaa / 4.0;
    double 	Kcc = 0.05;
    double 	Kac = 0.1;
    double 	Kw_c = Kcc / 4.0;
    int     P_start_c;
    int     P_start_e;
    int     P_end_c;
    int     P_end_e;
    int     P_num;
    int     R_divide_c;
    int     R_divide_e;
    int     R_lip_c;
    int     R_lip_e;
    int     T_divide_c;
    int     T_divide_e;
    int     periphery_size;
    double  phi_c;
    double  gamma_a_ratio;
    double  gamma_c_ratio;
    int     save_every_step;
    int     Nt;

    int     Na = (int)atoi(argv[1]);
    int     Nss = (int)atoi(argv[2]);
    int     Nsb = (int)atoi(argv[3]);
    double  KA = (double)atof(argv[4]);
    double  Kp = (double)atof(argv[5]);
    int     Nc = (int)atoi(argv[6]);
    double  delta = (double)atof(argv[7]);
    int     dAlf_c = (int)atoi(argv[8]);
    int     dAlf_e = (int)atoi(argv[9]);
    double  Kb = (double)atof(argv[10]);
    double  Kb_std = (double)atof(argv[11]);
    double  a_cut_c = (double)atof(argv[12]);
    double  Kcc_attr_ratio = (double)atof(argv[13]);
    double  a_cut_ac = (double)atof(argv[14]);
    double  Kac_attr_ratio = (double)atof(argv[15]);
    int     Pt_c = (int)atoi(argv[16]);
    int     Pt_e = (int)atoi(argv[17]);
    double  Kpin = (double)atof(argv[18]);
    int     pid = (int)atoi(argv[19]);
    int     is_smooth = (int)atoi(argv[20]);
    int     sim_type = (int)atoi(argv[21]);
    if (sim_type == -1){
        P_start_c = (int)atoi(argv[22]);
        P_start_e = (int)atoi(argv[23]);
        P_end_c = (int)atoi(argv[24]);
        P_end_e = (int)atoi(argv[25]);
        P_num = (int)atoi(argv[26]);
    }
    else if (sim_type == 1){
        R_divide_c = (int)atoi(argv[22]);
        R_divide_e = (int)atoi(argv[23]);
        gamma_a_ratio = (double)atof(argv[24]);
        gamma_c_ratio = (double)atof(argv[25]);
        save_every_step = (int)atoi(argv[26]);
        Nt = (int)atoi(argv[27]);
    }
    else if (sim_type == 2){
        R_divide_c = (int)atoi(argv[22]);
        R_divide_e = (int)atoi(argv[23]);
        R_lip_c = (int)atoi(argv[24]);
        R_lip_e = (int)atoi(argv[25]);
        gamma_a_ratio = (double)atof(argv[26]);
        gamma_c_ratio = (double)atof(argv[27]);
        save_every_step = (int)atoi(argv[28]);
        Nt = (int)atoi(argv[29]);
    }
    else if (sim_type == 3){
        periphery_size = (int)atoi(argv[22]);
        T_divide_c = (int)atoi(argv[23]);
        T_divide_e = (int)atoi(argv[24]);
        gamma_a_ratio = (double)atof(argv[25]);
        gamma_c_ratio = (double)atof(argv[26]);
        save_every_step = (int)atoi(argv[27]);
        Nt = (int)atoi(argv[28]);
    }
    else if (sim_type == 4){
        periphery_size = (int)atoi(argv[22]);
        T_divide_c = (int)atoi(argv[23]);
        T_divide_e = (int)atoi(argv[24]);
        R_lip_c = (int)atoi(argv[25]);
        R_lip_e = (int)atoi(argv[26]);
        gamma_a_ratio = (double)atof(argv[27]);
        gamma_c_ratio = (double)atof(argv[28]);
        save_every_step = (int)atoi(argv[29]);
        Nt = (int)atoi(argv[30]);
    }
    else if (sim_type == 5){
        phi_c = (double)atof(argv[22]);
        T_divide_c = (int)atoi(argv[23]);
        T_divide_e = (int)atoi(argv[24]);
        gamma_a_ratio = (double)atof(argv[25]);
        gamma_c_ratio = (double)atof(argv[26]);
        save_every_step = (int)atoi(argv[27]);
        Nt = (int)atoi(argv[28]);
    }
    else if (sim_type == 6){
        phi_c = (double)atof(argv[22]);
        T_divide_c = (int)atoi(argv[23]);
        T_divide_e = (int)atoi(argv[24]);
        R_lip_c = (int)atoi(argv[25]);
        R_lip_e = (int)atoi(argv[26]);
        gamma_a_ratio = (double)atof(argv[27]);
        gamma_c_ratio = (double)atof(argv[28]);
        save_every_step = (int)atoi(argv[29]);
        Nt = (int)atoi(argv[30]);
    }
    else if (sim_type == 7){
        phi_c = (double)atof(argv[22]);
        T_divide_c = (int)atoi(argv[23]);
        T_divide_e = (int)atoi(argv[24]);
        save_every_step = (int)atoi(argv[25]);
        Nt = (int)atoi(argv[26]);
    }
    else if (sim_type == 8){
        phi_c = (double)atof(argv[22]);
        T_divide_c = (int)atoi(argv[23]);
        T_divide_e = (int)atoi(argv[24]);
        R_lip_c = (int)atoi(argv[25]);
        R_lip_e = (int)atoi(argv[26]);
        save_every_step = (int)atoi(argv[27]);
        Nt = (int)atoi(argv[28]);
    }

	string			dir_prefix = "/gpfs/gibbs/pi/ohern/dw672/Adipocyte/RoughParticle/FlatInterface/Files/";

	BreastCancer_2D Model(Na, Nss, Nsb, Nc, delta, dAlf_c, dAlf_e, Kp, KA, Kaa, Kb, Kb_std, Kw_a,
                            Kpin, Kcc, a_cut_c, Kcc_attr_ratio, a_cut_ac,
                            Kac_attr_ratio, Kac, Kw_c, Pt_c, Pt_e,
                            is_smooth, pid, dir_prefix);

    if (sim_type == -1){
        Model.AdipocytePacking_Compress(P_start_c, P_start_e, P_end_c, P_end_e, P_num);
    }
    else if (sim_type == 0){
        Model.AdipocytePacking();
        // printf("Adipocyte Rigid done\n");
        Model.CancerPacking();
        // printf("Cancer rigid done\n");
        if (is_smooth == 0)
            Model.AllPacking();
        else
            Model.AllPacking_Smooth();
        // printf("All packing done\n");
    }
    else if(sim_type == 1){
        double  dt = 0.01, gamma_ac_ratio = 2.0;
        if (is_smooth == 0){
            Model.Invasion_Growth(R_divide_c, R_divide_e, gamma_a_ratio, gamma_c_ratio,
                                    gamma_ac_ratio, dt, save_every_step, Nt);
        }
        else {
            Model.Invasion_Growth_Smooth(R_divide_c, R_divide_e, gamma_a_ratio, gamma_c_ratio,
                                            gamma_ac_ratio, dt, save_every_step, Nt);
        }
    }
    else if (sim_type == 2){
        double  dt = 0.01, gamma_ac_ratio = 2.0;
        if (is_smooth == 0){
            Model.Invasion_Growth_LipidLoss(R_divide_c, R_divide_e, R_lip_c, R_lip_e, gamma_a_ratio,
                                            gamma_c_ratio, gamma_ac_ratio, dt, save_every_step, Nt);
        }
        else {
            Model.Invasion_Growth_LipidLoss_Smooth(R_divide_c, R_divide_e, R_lip_c, R_lip_e, gamma_a_ratio,
                                                    gamma_c_ratio, gamma_ac_ratio, dt, save_every_step, Nt);
        }
    }
    else if (sim_type == 3){
        double  dt = 0.01, gamma_ac_ratio = 2.0;
        if (is_smooth == 0){
            Model.Invasion_Periphery_Growth(periphery_size, T_divide_c, T_divide_e, gamma_a_ratio, gamma_c_ratio,
                                            gamma_ac_ratio, dt, save_every_step, Nt);
        }
        else {
            Model.Invasion_Periphery_Growth_Smooth(periphery_size, T_divide_c, T_divide_e, gamma_a_ratio, gamma_c_ratio,
                                                    gamma_ac_ratio, dt, save_every_step, Nt);
        }
    }
    else if (sim_type == 4){
        double  dt = 0.01, gamma_ac_ratio = 2.0;
        if (is_smooth == 0){
            Model.Invasion_Periphery_Growth_LipidLoss(periphery_size, T_divide_c, T_divide_e, R_lip_c, R_lip_e,
                                                        gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio,
                                                        dt, save_every_step, Nt);
        }
        else {
            Model.Invasion_Periphery_Growth_LipidLoss_Smooth(periphery_size, T_divide_c, T_divide_e, R_lip_c, R_lip_e,
                                                                gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio,
                                                                dt, save_every_step, Nt);
        }
    }
    else if (sim_type == 5){
        double  dt = 0.01, gamma_ac_ratio = 2.0;
        if (is_smooth == 0){
            Model.Invasion_Periphery_PF_Growth(phi_c, T_divide_c, T_divide_e, gamma_a_ratio, gamma_c_ratio,
                                            gamma_ac_ratio, dt, save_every_step, Nt);
        }
        else {
            Model.Invasion_Periphery_PF_Growth_Smooth(phi_c, T_divide_c, T_divide_e, gamma_a_ratio, gamma_c_ratio,
                                                    gamma_ac_ratio, dt, save_every_step, Nt);
        }
    }
    else if (sim_type == 6){
        double  dt = 0.01, gamma_ac_ratio = 2.0;
        if (is_smooth == 0){
            Model.Invasion_Periphery_PF_Growth_LipidLoss(phi_c, T_divide_c, T_divide_e, R_lip_c, R_lip_e,
                                                            gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio,
                                                            dt, save_every_step, Nt);
        }
        else {
            Model.Invasion_Periphery_PF_Growth_LipidLoss_Smooth(phi_c, T_divide_c, T_divide_e, R_lip_c, R_lip_e,
                                                                gamma_a_ratio, gamma_c_ratio, gamma_ac_ratio,
                                                                dt, save_every_step, Nt);
        }
    }
    else if (sim_type == 7){
        double  dt = 0.01;
        if (is_smooth == 0){
            Model.Invasion_Periphery_PF_Growth_OD(phi_c, T_divide_c, T_divide_e, dt, save_every_step, Nt);
        }
        else {
            Model.Invasion_Periphery_PF_Growth_OD_Smooth(phi_c, T_divide_c, T_divide_e, dt, save_every_step, Nt);
        }
    }
    else if (sim_type == 8){
        double  dt = 0.01, gamma_ac_ratio = 2.0;
        if (is_smooth == 0){
            Model.Invasion_Periphery_PF_Growth_LipidLoss_OD(phi_c, T_divide_c, T_divide_e, R_lip_c, R_lip_e,
                                                            dt, save_every_step, Nt);
        }
        else {
            Model.Invasion_Periphery_PF_Growth_LipidLoss_OD_Smooth(phi_c, T_divide_c, T_divide_e, R_lip_c, R_lip_e,
                                                                    dt, save_every_step, Nt);
        }
    }
    
    /*
    int     f0_c = 2, f0_e = -3, Nt = 400000;
    double  dt = 0.01, gamma_a_ratio = 0.2, gamma_ac_ratio = 2.0;
    Model.Invasion_ABP_Direction(f0_c, f0_e, gamma_a_ratio, gamma_ac_ratio, dt, Nt);
    */
	
    return 0;
}
