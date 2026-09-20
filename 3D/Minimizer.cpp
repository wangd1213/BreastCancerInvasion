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

void BreastCancer_3D::FIRE_Disk_VL()
{
    double      R_min = Rc[0];
    for (int nc = 1; nc < Nc; nc++)
        if (R_min > Rc[nc])
            R_min = Rc[nc];
    double      dt_fire = 0.01 * R_min / std::sqrt(Kcc);
    double 	    Fthresh = std::pow(10.0, -11) * R_min * R_min / Kcc;
    int         Nt_fire = 100000000;

	int 		n, nt = 0;
	double 		P, Vel_norm, Acc_norm;
    vector<double>  Vx_c(Nc, 0.0), Vy_c(Nc, 0.0), Vz_c(Nc, 0.0);
    // double      mw = 40.0 * std::sqrt((double)N / 400.0);
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;
	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL_pp(100 * Nc, vector<int> (2, 0));
	int			 			VL_count = 0;
	vector<double>			xc_save(Nc, 0.0), yc_save(Nc, 0.0), zc_save(Nc, 0.0);
	double 		r_cut = Rc[0];
	for (n = 1; n < Nc; n++)
		if (r_cut < Rc[n])
			r_cut = Rc[n];
	r_cut *= 2.0 * a_cut_c;

	// Verlet list initialization
	VerletList_Disk(r_cut, xc_save, yc_save, zc_save, VL_pp, VL_count, first_call);

	Force_Disk_VL(VL_pp, VL_count);
	//Pwall = Energy_Disk(N, Nall, Lx, Ly, Lz, R, pos, Acc);
	if (MaxForce_Disk() < Fthresh)
		return;
	
	first_call = 0;

	for (nt = 1; nt < Nt_fire; nt++){
		P = 0.0;
		for (n = 0; n < Nc; n++)
			P += Fx_c[n] * Vx_c[n] + Fy_c[n] * Vy_c[n] + Fz_c[n] * Vz_c[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max){
				break;
			}
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
			}
            for (n = 0; n < Nc; n++){
                xc[n] -= dt_half * Vx_c[n];
                yc[n] -= dt_half * Vy_c[n];
				zc[n] -= dt_half * Vz_c[n];
                Vx_c[n] = 0.0;
                Vy_c[n] = 0.0;
				Vz_c[n] = 0.0;
            }
		}

		// MD using Verlet method
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
            Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
        }
		Vel_norm = 0.0;
		Acc_norm = 0.0;
		for (n = 0; n < Nc; n++){
			Vel_norm += Vx_c[n] * Vx_c[n] + Vy_c[n] * Vy_c[n] + Vz_c[n] * Vz_c[n];
			Acc_norm += Fx_c[n] * Fx_c[n] + Fy_c[n] * Fy_c[n] + Fz_c[n] * Fz_c[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Nc; n++){
			Vx_c[n] = delta_a * Vx_c[n] + v_rsc * Fx_c[n];
            Vy_c[n] = delta_a * Vy_c[n] + v_rsc * Fy_c[n];
			Vz_c[n] = delta_a * Vz_c[n] + v_rsc * Fz_c[n];
			xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n];
			zc[n] += dt * Vz_c[n];
		}
		VerletList_Disk(r_cut, xc_save, yc_save, zc_save, VL_pp, VL_count, first_call);
		Force_Disk_VL(VL_pp, VL_count);
		//Pwall = Energy_Disk(N, Nall, Lx, Ly, Lz, R, pos, Acc);
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
            Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
        }

		// if (nt % 10000 == 0)
            // printf("nt: %d  Max Force: %.5e\n", nt, MaxForce_Disk());

		if (MaxForce_Disk() < Fthresh){
			//printf("nt: %d  Max Force: %.5e  P: %.5e\n", nt, Acc_max, Acc[Nall]);
			return;
		}
	}
}

void BreastCancer_3D::FIRE_Disk_Enthalpy_VL(double Pt)
{
    double      R_min = Rc[0];
    for (int nc = 1; nc < Nc; nc++)
        if (R_min > Rc[nc])
            R_min = Rc[nc];
    double      dt_fire = 0.01 * R_min / std::sqrt(Kcc);
    double 	    Fthresh = std::pow(10.0, -11) * R_min * R_min / Kcc;
    int         Nt_fire = 100000000;

	int 		n, nt = 0;
	double 		P, Vel_norm, Acc_norm;
    double      mw = 300.0;
    Vw = 0.0;
	double 		Ft = Pt * Lx * Lz;
    vector<double>  Vx_c(Nc, 0.0), Vy_c(Nc, 0.0), Vz_c(Nc, 0.0);
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;
	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL_pp(100 * Nc, vector<int> (2, 0));
	int			 			VL_count = 0;
	vector<double>			xc_save(Nc, 0.0), yc_save(Nc, 0.0), zc_save(Nc, 0.0);
	double 		r_cut = Rc[0];
	for (n = 1; n < Nc; n++)
		if (r_cut < Rc[n])
			r_cut = Rc[n];
	r_cut *= 2.0 * a_cut_c;

	// Verlet list initialization
	VerletList_Disk(r_cut, xc_save, yc_save, zc_save, VL_pp, VL_count, first_call);

	Force_Disk_VL(VL_pp, VL_count);
    Fw -= Ft;
	//Pwall = Energy_Disk(N, Nall, Lx, Ly, Lz, R, pos, Acc);
	double Acc_max = MaxForce_Disk();
	if ((Acc_max < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Ft))
		return;
	
	first_call = 0;

	for (nt = 1; nt < Nt_fire; nt++){
		P = Fw * Vw;
		for (n = 0; n < Nc; n++)
			P += Fx_c[n] * Vx_c[n] + Fy_c[n] * Vy_c[n] + Fz_c[n] * Vz_c[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max){
				break;
			}
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
			}
            for (n = 0; n < Nc; n++){
                xc[n] -= dt_half * Vx_c[n];
                yc[n] -= dt_half * Vy_c[n];
				zc[n] -= dt_half * Vz_c[n];
                Vx_c[n] = 0.0;
                Vy_c[n] = 0.0;
				Vz_c[n] = 0.0;
            }
			Ly -= dt_half * Vw;
			Vw = 0.0;
		}

		// MD using Verlet method
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
            Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
        }
        Vw += dt_half * Fw / mw;
		Vel_norm = Vw * Vw;
		Acc_norm = Fw * Fw;
		for (n = 0; n < Nc; n++){
			Vel_norm += Vx_c[n] * Vx_c[n] + Vy_c[n] * Vy_c[n] + Vz_c[n] * Vz_c[n];
			Acc_norm += Fx_c[n] * Fx_c[n] + Fy_c[n] * Fy_c[n] + Fz_c[n] * Fz_c[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Nc; n++){
			Vx_c[n] = delta_a * Vx_c[n] + v_rsc * Fx_c[n];
            Vy_c[n] = delta_a * Vy_c[n] + v_rsc * Fy_c[n];
			Vz_c[n] = delta_a * Vz_c[n] + v_rsc * Fz_c[n];
			xc[n] += dt * Vx_c[n];
            yc[n] += dt * Vy_c[n];
			zc[n] += dt * Vz_c[n];
		}
        Vw = delta_a * Vw + v_rsc * Fw;
        Ly += dt * Vw;
		VerletList_Disk(r_cut, xc_save, yc_save, zc_save, VL_pp, VL_count, first_call);
		Force_Disk_VL(VL_pp, VL_count);
        Fw -= Ft;

		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
            Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
        }
        Vw += dt_half * Fw / mw;

		// if (nt % 10000 == 0)
            // printf("nt: %d  Max Force: %.5e  Fw: %.5e  Ly: %.5e\n", nt, MaxForce_Disk(), Fw, Ly);

		if ((MaxForce_Disk() < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Ft))
			//printf("nt: %d  Max Force: %.5e  P: %.5e\n", nt, Acc_max, Acc[Nall]);
			return;
	}
}

void BreastCancer_3D::FIRE_DPM_VL()
{
    double      Fthresh = std::pow(10.0, -11);
    double      dt_fire = 0.01;
    int         Nt_fire = 100000000;

	int 		n, nt = 0;
	double 		P, Vel_norm, Acc_norm;
    vector<double>  Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0), Vz_a(Ns_tot, 0.0);
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;

	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL_aa(100 * Ns_tot, vector<int> (4, -1));
	int     		 		VL_counter = 0;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0), za_save(Ns_tot, 0.0);
	double 		r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
	r_cut *= 2.0;
	
    VerletList_DPM(r_cut, xa_save, ya_save, za_save, VL_aa, VL_counter, first_call);
	Force_DPM_VL(VL_aa, VL_counter);

	if (MaxForce_DPM() < Fthresh){
		// printf("Acc_max: %.4e  Iteration Number: %d\n", MaxForce_DPM(), nt);
		return;
	}

	first_call = 0;
	for (nt = 1; nt < Nt_fire; nt++){
		P = 0.0;
		for (n = 0; n < Ns_tot; n++)
			P += Fx_a[n] * Vx_a[n] + Fy_a[n] * Vy_a[n] + Fz_a[n] * Vz_a[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max)
				break;
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
			}
            for (n = 0; n < Ns_tot; n++){
                xa[n] -= dt_half * Vx_a[n];
                ya[n] -= dt_half * Vy_a[n];
				za[n] =- dt_half * Vz_a[n];
                Vx_a[n] = 0.0;
                Vy_a[n] = 0.0;
				Vz_a[n] = 0.0;
            }
		}

		// MD using Verlet method
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }
		Vel_norm = 0.0;
		Acc_norm = 0.0;
		for (n = 0; n < Ns_tot; n++){
			Vel_norm += Vx_a[n] * Vx_a[n] + Vy_a[n] * Vy_a[n] + Vz_a[n] * Vz_a[n];
			Acc_norm += Fx_a[n] * Fx_a[n] + Fy_a[n] * Fy_a[n] + Fz_a[n] * Fz_a[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] = delta_a * Vx_a[n] + v_rsc * Fx_a[n];
            Vy_a[n] = delta_a * Vy_a[n] + v_rsc * Fy_a[n];
			Vz_a[n] = delta_a * Vz_a[n] + v_rsc * Fz_a[n];
			xa[n] += dt * Vx_a[n];
            ya[n] += dt * Vy_a[n];
			za[n] += dt * Vz_a[n];
		}
		
		VerletList_DPM(r_cut, xa_save, ya_save, za_save, VL_aa, VL_counter, first_call);
	    Force_DPM_VL(VL_aa, VL_counter);
		
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }

		if (MaxForce_DPM() < Fthresh)
			break;
		// if (nt % 10 == 0) printf("nt: %d  Max Force: %.5e\n", nt, MaxForce_DPM());
	}
	// printf("Acc_max: %.4e  Iteration Number: %d\n", MaxForce_DPM(), nt);
	// printf("VL_count: %d\n", VL_counter);
}

void BreastCancer_3D::FIRE_DPM_Enthalpy_VL(double Pt)
{
    double      Fthresh = std::pow(10.0, -11);
    double      dt_fire = 0.01;
    int         Nt_fire = 100000000;

	int 		n, nt = 0;
    double      mw = 300.0;
    Vw = 0.0;
    vector<double>  Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0), Vz_a(Ns_tot, 0.0);
	double 		P, Vel_norm, Acc_norm;
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;

	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL_aa(100 * Ns_tot, vector<int> (4, -1));
	int     		 		VL_counter = 0;
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0), za_save(Ns_tot, 0.0);
	double 		r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
	r_cut *= 2.0;
	
    VerletList_DPM(r_cut, xa_save, ya_save, za_save, VL_aa, VL_counter, first_call);
	Force_DPM_VL(VL_aa, VL_counter);
	Fw -= Pt * Ly * Ly;

	double Acc_max = MaxForce_DPM();
	if ((Acc_max < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Pt * Ly * Ly)){
		//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
		return;
	}

	first_call = 0;
	for (nt = 1; nt < Nt_fire; nt++){
		P = Fw * Vw;
		for (n = 0; n < Ns_tot; n++)
			P += Fx_a[n] * Vx_a[n] + Fy_a[n] * Vy_a[n] + Fz_a[n] * Vz_a[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max)
				break;
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
			}
            for (n = 0; n < Ns_tot; n++){
                xa[n] -= dt_half * Vx_a[n];
                ya[n] -= dt_half * Vy_a[n];
				za[n] -= dt_half * Vz_a[n];
                Vx_a[n] = 0.0;
                Vy_a[n] = 0.0;
				Vz_a[n] = 0.0;
                Ly -= dt_half * Vw;
				Lx = Ly;
				Lz = Ly;
                Vw = 0.0;
            }
		}

		// MD using Verlet method
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }
        Vw += dt_half * Fw / mw;
		Vel_norm = Vw * Vw;
		Acc_norm = Fw * Fw;
		for (n = 0; n < Ns_tot; n++){
			Vel_norm += Vx_a[n] * Vx_a[n] + Vy_a[n] * Vy_a[n] + Vz_a[n] * Vz_a[n];
			Acc_norm += Fx_a[n] * Fx_a[n] + Fy_a[n] * Fy_a[n] + Fz_a[n] * Fz_a[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] = delta_a * Vx_a[n] + v_rsc * Fx_a[n];
            Vy_a[n] = delta_a * Vy_a[n] + v_rsc * Fy_a[n];
			Vz_a[n] = delta_a * Vz_a[n] + v_rsc * Fz_a[n];
			xa[n] += dt * Vx_a[n];
            ya[n] += dt * Vy_a[n];
			za[n] += dt * Vz_a[n];
		}
        Vw = delta_a * Vw + v_rsc * Fw;
        Ly += dt * Vw;
		Lx = Ly;
		Lz = Ly;
		
		VerletList_DPM(r_cut, xa_save, ya_save, za_save, VL_aa, VL_counter, first_call);
        Force_DPM_VL(VL_aa, VL_counter);
        Fw -= Pt * Ly * Ly;

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }
        Vw += dt_half * Fw / mw;

		Acc_max = MaxForce_DPM();
		if ((Acc_max < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Pt * Ly * Ly))
			return;
		// if (nt % 1 == 0) printf("nt: %d  Max Force: %.5e\n", nt, Acc_max);
	}
	//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
}

void BreastCancer_3D::FIRE_DPM_Disk_VL()
{
    double      Fthresh = std::pow(10.0, -11);
    double      dt_fire = 0.01;
    int         Nt_fire = 100000000;

	int 		n, nt = 0;
    vector<double>  Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0), Vz_a(Ns_tot, 0.0);
	vector<double> 	Vx_c(Nc, 0.0), Vy_c(Nc, 0.0), Vz_c(Nc, 0.0);
	double 		P, Vel_norm, Acc_norm;
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;

	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL_aa(400 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(400 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(200 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0), za_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0), zc_save(Nc, 0.0);
	double 					r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
	r_cut *= a_cut_ac;
	for (n = 0; n < Nc; n++)
		r_cut = (r_cut > a_cut_c * Rc[n]) ? r_cut : a_cut_c * Rc[n];
    r_cut *= 2.0;
	
    VerletList_DPM_Disk(r_cut, xa_save, ya_save, za_save, xc_save, yc_save, zc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

	if (MaxForce_All() < Fthresh){
		//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
		return;
	}

	first_call = 0;
	for (nt = 1; nt < Nt_fire; nt++){
		P = 0.0;
		for (n = 0; n < Ns_tot; n++)
			P += Fx_a[n] * Vx_a[n] + Fy_a[n] * Vy_a[n] + Fz_a[n] * Vz_a[n];
		for (n = 0; n < Nc; n++)
			P += Fx_c[n] * Vx_c[n] + Fy_c[n] * Vy_c[n] + Fz_c[n] * Vz_c[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max)
				break;
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
			}
            for (n = 0; n < Ns_tot; n++){
                xa[n] -= dt_half * Vx_a[n];
                ya[n] -= dt_half * Vy_a[n];
				za[n] -= dt_half * Vz_a[n];
                Vx_a[n] = 0.0;
                Vy_a[n] = 0.0;
				Vz_a[n] = 0.0;
            }
			for (n = 0; n < Nc; n++){
				xc[n] -= dt_half * Vx_c[n];
				yc[n] -= dt_half * Vy_c[n];
				zc[n] -= dt_half * Vz_c[n];
				Vx_c[n] = 0.0;
				Vy_c[n] = 0.0;
				Vz_c[n] = 0.0;
			}
		}

		// MD using Verlet method
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
			Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
		}
		Vel_norm = 0.0;
		Acc_norm = 0.0;
		for (n = 0; n < Ns_tot; n++){
			Vel_norm += Vx_a[n] * Vx_a[n] + Vy_a[n] * Vy_a[n] + Vz_a[n] * Vz_a[n];
			Acc_norm += Fx_a[n] * Fx_a[n] + Fy_a[n] * Fy_a[n] + Fz_a[n] * Fz_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vel_norm += Vx_c[n] * Vx_c[n] + Vy_c[n] * Vy_c[n] + Vz_c[n] * Vz_c[n];
			Acc_norm += Fx_c[n] * Fx_c[n] + Fy_c[n] * Fy_c[n] + Fz_c[n] * Fz_c[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] = delta_a * Vx_a[n] + v_rsc * Fx_a[n];
            Vy_a[n] = delta_a * Vy_a[n] + v_rsc * Fy_a[n];
			Vz_a[n] = delta_a * Vz_a[n] + v_rsc * Fz_a[n];
			xa[n] += dt * Vx_a[n];
            ya[n] += dt * Vy_a[n];
			za[n] += dt * Vz_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] = delta_a * Vx_c[n] + v_rsc * Fx_c[n];
			Vy_c[n] = delta_a * Vy_c[n] + v_rsc * Fy_c[n];
			Vz_c[n] = delta_a * Vz_c[n] + v_rsc * Fz_c[n];
			xc[n] += dt * Vx_c[n];
			yc[n] += dt * Vy_c[n];
			zc[n] += dt * Vz_c[n];
		}
		
		VerletList_DPM_Disk(r_cut, xa_save, ya_save, za_save, xc_save, yc_save, zc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_VL(VL_aa, VL_cc, VL_ac, VL_count_all);

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
			Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
		}

		if (MaxForce_All() < Fthresh)
			return;
		// if (nt % 10000 == 0) printf("nt: %d  Max Force: %.5e\n", nt, MaxForce_All());
	}
	//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
}

void BreastCancer_3D::FIRE_DPM_FixDisk_Enthalpy_VL(double Pt)
{
    double      Fthresh = std::pow(10.0, -11);
    double      dt_fire = 0.01;
    int         Nt_fire = 100000000;

	int 		n, nt = 0;
    double      mw = 300.0;
    Vw = 0.0;
	double 		Ft = Pt * Lx * Lz;
    vector<double>  Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0), Vz_a(Ns_tot, 0.0);
	double 		P, Vel_norm, Acc_norm;
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;

	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL_aa(400 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(400 * Ns_tot, vector<int> (3, -1));
	vector<int>		 		VL_count_all(2, 0);
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0), za_save(Ns_tot, 0.0);
	double 					r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
	r_cut *= a_cut_ac;
	for (n = 0; n < Nc; n++)
		r_cut = (r_cut > Rc[n]) ? r_cut : Rc[n];
    r_cut *= 2.0;
	
    VerletList_DPM_FixDisk(r_cut, xa_save, ya_save, za_save, VL_aa, VL_ac, VL_count_all, first_call);
	Force_DPM_FixDisk_VL(VL_aa, VL_ac, VL_count_all);
	Fw -= Ft;

	if ((MaxForce_DPM() < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Ft)){
		//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
		return;
	}

	first_call = 0;
	for (nt = 1; nt < Nt_fire; nt++){
		P = Fw * Vw;
		for (n = 0; n < Ns_tot; n++)
			P += Fx_a[n] * Vx_a[n] + Fy_a[n] * Vy_a[n] + Fz_a[n] * Vz_a[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max)
				break;
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
			}
            for (n = 0; n < Ns_tot; n++){
                xa[n] -= dt_half * Vx_a[n];
                ya[n] -= dt_half * Vy_a[n];
				za[n] -= dt_half * Vz_a[n];
                Vx_a[n] = 0.0;
                Vy_a[n] = 0.0;
				Vz_a[n] = 0.0;
            }
			Ly -= dt_half * Vw;
			Vw = 0.0;
		}

		// MD using Verlet method
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }
        Vw += dt_half * Fw / mw;
		Vel_norm = Vw * Vw;
		Acc_norm = Fw * Fw;
		for (n = 0; n < Ns_tot; n++){
			Vel_norm += Vx_a[n] * Vx_a[n] + Vy_a[n] * Vy_a[n] + Vz_a[n] * Vz_a[n];
			Acc_norm += Fx_a[n] * Fx_a[n] + Fy_a[n] * Fy_a[n] + Fz_a[n] * Fz_a[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] = delta_a * Vx_a[n] + v_rsc * Fx_a[n];
            Vy_a[n] = delta_a * Vy_a[n] + v_rsc * Fy_a[n];
			Vz_a[n] = delta_a * Vz_a[n] + v_rsc * Fz_a[n];
			xa[n] += dt * Vx_a[n];
            ya[n] += dt * Vy_a[n];
			za[n] += dt * Vz_a[n];
		}
        Vw = delta_a * Vw + v_rsc * Fw;
        Ly += dt * Vw;
		
		VerletList_DPM_FixDisk(r_cut, xa_save, ya_save, za_save, VL_aa, VL_ac, VL_count_all, first_call);
        Force_DPM_FixDisk_VL(VL_aa, VL_ac, VL_count_all);
        Fw -= Ft;

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }
        Vw += dt_half * Fw / mw;

		if ((MaxForce_DPM() < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Ft))
			return;
		// if (nt % 10000 == 0) printf("nt: %d  Max Force: %.5e\n", nt, MaxForce_All());
	}
	//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
}

void BreastCancer_3D::FIRE_Disk_FixDPM_Enthalpy_VL(double Pt)
{
    double      Fthresh = std::pow(10.0, -11);
    double      dt_fire = 0.01;
    int         Nt_fire = 100000000;

	int 		n, nt = 0;
    double      mw = 300.0;
    Vw = 0.0;
	double 		Ft = Pt * Lx * Lz;
	vector<double> 	Vx_c(Nc, 0.0), Vy_c(Nc, 0.0), Vz_c(Nc, 0.0);
	double 		P, Vel_norm, Acc_norm;
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;

	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>> 	VL_ac(400 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(200 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(2, 0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0), zc_save(Nc, 0.0);
	double 					r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
	r_cut *= a_cut_ac;
	for (n = 0; n < Nc; n++)
		r_cut = (r_cut > a_cut_c * Rc[n]) ? r_cut : a_cut_c * Rc[n];
    r_cut *= 2.0;
	
    VerletList_Disk_FixDPM(r_cut, xc_save, yc_save, zc_save, VL_cc, VL_ac, VL_count_all, first_call);
	Force_Disk_FixDPM_VL(VL_cc, VL_ac, VL_count_all);
	Fw -= Ft;

	if ((MaxForce_Disk() < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Ft)){
		//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
		return;
	}

	first_call = 0;
	for (nt = 1; nt < Nt_fire; nt++){
		P = Fw * Vw;
		for (n = 0; n < Nc; n++)
			P += Fx_c[n] * Vx_c[n] + Fy_c[n] * Vy_c[n] + Fz_c[n] * Vz_c[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max)
				break;
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
			}
			for (n = 0; n < Nc; n++){
				xc[n] -= dt_half * Vx_c[n];
				yc[n] -= dt_half * Vy_c[n];
				zc[n] -= dt_half * Vz_c[n];
				Vx_c[n] = 0.0;
				Vy_c[n] = 0.0;
				Vz_c[n] = 0.0;
			}
			Ly -= dt_half * Vw;
			Vw = 0.0;
		}

		// MD using Verlet method
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
			Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
		}
        Vw += dt_half * Fw / mw;
		Vel_norm = Vw * Vw;
		Acc_norm = Fw * Fw;
		for (n = 0; n < Nc; n++){
			Vel_norm += Vx_c[n] * Vx_c[n] + Vy_c[n] * Vy_c[n] + Vz_c[n] * Vz_c[n];
			Acc_norm += Fx_c[n] * Fx_c[n] + Fy_c[n] * Fy_c[n] + Fz_c[n] * Fz_c[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Nc; n++){
			Vx_c[n] = delta_a * Vx_c[n] + v_rsc * Fx_c[n];
			Vy_c[n] = delta_a * Vy_c[n] + v_rsc * Fy_c[n];
			Vz_c[n] = delta_a * Vz_c[n] + v_rsc * Fz_c[n];
			xc[n] += dt * Vx_c[n];
			yc[n] += dt * Vy_c[n];
			zc[n] += dt * Vz_c[n];
		}
        Vw = delta_a * Vw + v_rsc * Fw;
        Ly += dt * Vw;
		
		VerletList_Disk_FixDPM(r_cut, xc_save, yc_save, zc_save, VL_cc, VL_ac, VL_count_all, first_call);
        Force_Disk_FixDPM_VL(VL_cc, VL_ac, VL_count_all);
        Fw -= Ft;

		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
			Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
		}
        Vw += dt_half * Fw / mw;

		if ((MaxForce_Disk() < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Ft))
			return;
		// if (nt % 10000 == 0) printf("nt: %d  Max Force: %.5e\n", nt, MaxForce_All());
	}
	//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
}

void BreastCancer_3D::FIRE_DPM_Disk_Enthalpy_VL(double Pt)
{
    double      Fthresh = std::pow(10.0, -11);
    double      dt_fire = 0.001;
    int         Nt_fire = 100000000;

	int 		n, nt = 0;
    double      mw = 300.0;
    Vw = 0.0;
	double 		Ft = Pt * Lx * Lz;
    vector<double>  Vx_a(Ns_tot, 0.0), Vy_a(Ns_tot, 0.0), Vz_a(Ns_tot, 0.0);
	vector<double> 	Vx_c(Nc, 0.0), Vy_c(Nc, 0.0), Vz_c(Nc, 0.0);
	double 		P, Vel_norm, Acc_norm;
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;

	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL_aa(400 * Ns_tot, vector<int> (4, -1));
	vector<vector <int>> 	VL_ac(400 * Ns_tot, vector<int> (3, -1));
	vector<vector <int>> 	VL_cc(200 * Nc, vector<int> (2, -1));
	vector<int>		 		VL_count_all(3, 0);
	vector<double>			xa_save(Ns_tot, 0.0), ya_save(Ns_tot, 0.0), za_save(Ns_tot, 0.0);
    vector<double>          xc_save(Nc, 0.0), yc_save(Nc, 0.0), zc_save(Nc, 0.0);
	double 					r_cut = R0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > R0[n]) ? r_cut : R0[n];
	r_cut *= a_cut_ac;
	for (n = 0; n < Nc; n++)
		r_cut = (r_cut > a_cut_c * Rc[n]) ? r_cut : a_cut_c * Rc[n];
    r_cut *= 2.0;
	
    VerletList_DPM_Disk(r_cut, xa_save, ya_save, za_save, xc_save, yc_save, zc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
	Force_DPM_Disk_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
	Fw -= Ft;

	if ((MaxForce_All() < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Ft)){
		//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
		return;
	}

	first_call = 0;
	for (nt = 1; nt < Nt_fire; nt++){
		P = Fw * Vw;
		for (n = 0; n < Ns_tot; n++)
			P += Fx_a[n] * Vx_a[n] + Fy_a[n] * Vy_a[n] + Fz_a[n] * Vz_a[n];
		for (n = 0; n < Nc; n++)
			P += Fx_c[n] * Vx_c[n] + Fy_c[n] * Vy_c[n] + Fz_c[n] * Vz_c[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max)
				break;
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
			}
            for (n = 0; n < Ns_tot; n++){
                xa[n] -= dt_half * Vx_a[n];
                ya[n] -= dt_half * Vy_a[n];
				za[n] -= dt_half * Vz_a[n];
                Vx_a[n] = 0.0;
                Vy_a[n] = 0.0;
				Vz_a[n] = 0.0;
            }
			for (n = 0; n < Nc; n++){
				xc[n] -= dt_half * Vx_c[n];
				yc[n] -= dt_half * Vy_c[n];
				zc[n] -= dt_half * Vz_c[n];
				Vx_c[n] = 0.0;
				Vy_c[n] = 0.0;
				Vz_c[n] = 0.0;
			}
			Ly -= dt_half * Vw;
			Vw = 0.0;
		}

		// MD using Verlet method
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
			Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
		}
        Vw += dt_half * Fw / mw;
		Vel_norm = Vw * Vw;
		Acc_norm = Fw * Fw;
		for (n = 0; n < Ns_tot; n++){
			Vel_norm += Vx_a[n] * Vx_a[n] + Vy_a[n] * Vy_a[n] + Vz_a[n] * Vz_a[n];
			Acc_norm += Fx_a[n] * Fx_a[n] + Fy_a[n] * Fy_a[n] + Fz_a[n] * Fz_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vel_norm += Vx_c[n] * Vx_c[n] + Vy_c[n] * Vy_c[n] + Vz_c[n] * Vz_c[n];
			Acc_norm += Fx_c[n] * Fx_c[n] + Fy_c[n] * Fy_c[n] + Fz_c[n] * Fz_c[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] = delta_a * Vx_a[n] + v_rsc * Fx_a[n];
            Vy_a[n] = delta_a * Vy_a[n] + v_rsc * Fy_a[n];
			Vz_a[n] = delta_a * Vz_a[n] + v_rsc * Fz_a[n];
			xa[n] += dt * Vx_a[n];
            ya[n] += dt * Vy_a[n];
			za[n] += dt * Vz_a[n];
		}
		for (n = 0; n < Nc; n++){
			Vx_c[n] = delta_a * Vx_c[n] + v_rsc * Fx_c[n];
			Vy_c[n] = delta_a * Vy_c[n] + v_rsc * Fy_c[n];
			Vz_c[n] = delta_a * Vz_c[n] + v_rsc * Fz_c[n];
			xc[n] += dt * Vx_c[n];
			yc[n] += dt * Vy_c[n];
			zc[n] += dt * Vz_c[n];
		}
        Vw = delta_a * Vw + v_rsc * Fw;
        Ly += dt * Vw;
		
		VerletList_DPM_Disk(r_cut, xa_save, ya_save, za_save, xc_save, yc_save, zc_save, VL_aa, VL_cc, VL_ac, VL_count_all, first_call);
        Force_DPM_Disk_VL(VL_aa, VL_cc, VL_ac, VL_count_all);
        Fw -= Ft;

		for (n = 0; n < Ns_tot; n++){
			Vx_a[n] += dt_half * Fx_a[n];
            Vy_a[n] += dt_half * Fy_a[n];
			Vz_a[n] += dt_half * Fz_a[n];
        }
		for (n = 0; n < Nc; n++){
			Vx_c[n] += dt_half * Fx_c[n];
			Vy_c[n] += dt_half * Fy_c[n];
			Vz_c[n] += dt_half * Fz_c[n];
		}
        Vw += dt_half * Fw / mw;

		if ((MaxForce_All() < Fthresh) && (std::abs(Fw) < std::pow(10.0, -5) * Ft)){
			// printf("nt: %d\n", nt);
			return;
		}
		// if (nt % 10000 == 0) printf("nt: %d  Max Force: %.5e  Fw: %.5e\n", nt, MaxForce_All(), Fw);
	}
	//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
}
