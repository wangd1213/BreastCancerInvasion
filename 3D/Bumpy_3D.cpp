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

#include "Bumpy_3D.h"

using std::vector;
using std::string;

// Constructor implementation
Bumpy_3D::Bumpy_3D(int N_set,
					int Ns_set,
					double delta_set,
					double K_set,
					double Kw_set,
					int Pt_c_set,
					int Pt_e_set,
					int seed,
					string dir_pre)
    : N(N_set), Ns(Ns_set),
	  delta(delta_set), K(K_set), Kw(Kw_set),
	  pid(seed), Pt_c(Pt_c_set), Pt_e(Pt_e_set)
{
    Ns_tot = Ns * N;
    x.assign(N, 0.0);
    y.assign(N, 0.0);
	z.assign(N, 0.0);
	Q0.assign(N, 0.0);
	Q1.assign(N, 0.0);
	Q2.assign(N, 0.0);
	Q3.assign(N, 0.0);
    x_bump.assign(Ns_tot, 0.0);
    y_bump.assign(Ns_tot, 0.0);
    z_bump.assign(Ns_tot, 0.0);
    X_unit.assign(Ns, 0.0);
    Y_unit.assign(Ns, 0.0);
    Z_unit.assign(Ns, 0.0);
	idx_list.assign(N + 1, 0);
	Fx.assign(N, 0.0);
	Fy.assign(N, 0.0);
	Fz.assign(N, 0.0);
	Tx.assign(N, 0.0);
	Ty.assign(N, 0.0);
	Tz.assign(N, 0.0);

	for (int n = 1; n <= N; n++)
		idx_list[n] = n * Ns;

	Pt = (double)Pt_c * std::pow(10.0, Pt_e);

	std::ostringstream oss;
	oss << dir_pre << "Polyhedron/Ns_" << Ns << "/DVA.txt";
    dva_name = oss.str();
    oss.str("");
    oss.clear();

	oss << dir_pre << "Polyhedron/Ns_" << Ns << "/RigidUnit.txt";
    rigid_unit_name = oss.str();
    oss.str("");
    oss.clear();

	oss << "P_" << Pt_c << "E" << Pt_e;
    P_str = oss.str();
    oss.str("");
    oss.clear();

    oss << dir_pre << "Rigid/Na_" << std::setfill('0') << std::setw(3) << N
		<< "_Ns_" << std::setfill('0') << std::setw(3) << Ns << "/"
        << "/Dr_" << std::fixed << std::setprecision(2) << delta << "/";
    dir_run = oss.str();
    oss.str("");
    oss.clear();

	string cmd;
	oss << "mkdir -p " << dir_run;
    cmd = oss.str();
    oss.str("");
    oss.clear();
	std::system(cmd.c_str());

	oss << dir_run << "Pos_Rigid_" << std::setfill('0') << std::setw(5) << pid << ".txt";
    rigid_name = oss.str();
    oss.str("");
    oss.clear();

	SetParticleParameters();
}

void Bumpy_3D::SetParticleParameters()
{
	double 	D0_unit, Alf_min, A_sum_unit;
    LoadUnitParameter(D0_unit, Alf_min, A_sum_unit);
	double 	V0 = 1.0;
    R = std::pow(V0 * Alf_min / (std::pow(A_sum_unit, 1.5) / std::sqrt(PI) / 6.0), 1.0 / 3.0);
	D_bump = delta * D0_unit * R;
    R_bump = delta * 0.5 * D_bump;
    LoadUnitBumpy();
}

void Bumpy_3D::LoadUnitParameter(double &D0_unit,
                                    double &Alf_min,
                                    double &A_sum_unit)
{
	std::ifstream pfile;
	pfile.open(dva_name.c_str());
    if (pfile.good()){
    	pfile >> D0_unit;
        double A0_temp;
        for(int i = 0; i < 2 * Ns - 4; i++)
            pfile >> A0_temp;
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
}

void Bumpy_3D::LoadUnitBumpy()
{
	int 			ns;
	double 			D_bump_unit;

	std::ifstream file1(rigid_unit_name.c_str());
    if (file1.good()){
		file1 >> D_bump_unit;
		file1 >> D_max;
		for (ns = 0; ns < Ns; ns++)
			file1 >> X_unit[ns];
		for (ns = 0; ns < Ns; ns++)
			file1 >> Y_unit[ns];
		for (ns = 0; ns < Ns; ns++)
			file1 >> Z_unit[ns];
		file1.close();

        D_max *= R;
        for (ns = 0; ns < Ns; ns++){
            X_unit[ns] *= R;
            Y_unit[ns] *= R;
            Z_unit[ns] *= R;
        }
    }
	else {
		file1.close();
        exit(-1);
	}
}

void Bumpy_3D::Initialization(double ipf)
{
    L = 2.0 * std::pow(N / ipf, 1.0 / 3.0);
	
    std::mt19937 gen(pid);
    std::uniform_real_distribution<> dis_uniform(0.0, 1.0);

    vector<double>  Ccen(3 * N, 0.0);
    while (true){
		for (int n = 0; n < N; n++){
			x[n] = L * dis_uniform(gen);
			y[n] = (L - 2.0 * R) * dis_uniform(gen) + R;
			z[n] = L * dis_uniform(gen);
		}
		if (!overlap(x, y, z))
			break;
	}

    double q_theta, q_phi, q_psi; // Euler angles
	double sin_q_theta, cos_q_theta;
	double sin_q_phi, cos_q_phi;
	double sin_q_psi, cos_q_psi;

    for (int n = 0; n < N; n++){
		q_theta = PI * dis_uniform(gen);
		q_phi = 2.0 * PI * dis_uniform(gen);
		q_psi = PI * dis_uniform(gen);
		sin_q_theta = std::sin(q_theta);
		cos_q_theta = std::cos(q_theta);
		sin_q_phi = std::sin(q_phi);
		cos_q_phi = std::cos(q_phi);
		sin_q_psi = std::sin(q_psi);
		cos_q_psi = std::cos(q_psi);
		Q0[n] = cos_q_theta;
		Q1[n] = sin_q_theta * sin_q_psi * sin_q_phi;
		Q2[n] = sin_q_theta * sin_q_psi * cos_q_phi;
		Q3[n] = sin_q_theta * cos_q_psi;
    }
}

bool Bumpy_3D::overlap(vector<double> &x,
						vector<double> &y,
						vector<double> &z)
{
	double dx, dy, dz, dr;
	
	for (int n = 0; n < N - 1; n++){
		for (int m = n + 1; m < N; m++){
			dx = x[m] - x[n];
			dx -= std::round(dx / L) * L;
			if (std::abs(dx) < 4.0 * R){
				dy = y[m] - y[n];
				if (std::abs(dy) < 4.0 * R){
					dz = z[m] - z[n];
					dz -= std::round(dz / L) * L;
					if ((std::abs(dx) < 4.0 * R) && (dx * dx + dy * dy + dz * dz < 16.0 * R * R))
                        return true;
				}
			}
		}
	}

	for (int n = 0; n < N; n++)
		if ((y[n] > L - R) || (y[n] < R))
			return true;

	return false;

}

void Bumpy_3D::Packing_AboveJamming(double ipf)
{
	std::ifstream posfile(rigid_name);
    if (posfile.good()){
        posfile.close();
        return;
    }
    else
        posfile.close();

	Initialization(ipf);

    int 		n;
	double 		rsc = 1.001;
    double 		Et_l = N * std::pow(10.0, -14);

	while (true){
		FIRE_VL();
        for (n = 0; n < N; n++){
			x[n] -= std::floor(x[n] / L) * L;
			z[n] -= std::floor(z[n] / L) * L;
        }

		// printf("L: %.5e E: %.5e\n", L, Energy());

		if (Energy() < Et_l)
            L /= rsc;
		else
			break;

		for (n = 0; n < N; n++){
			x[n] /= rsc;
            y[n] /= rsc;
            z[n] /= rsc;
        }
	}

    SaveConfig();
}

void Bumpy_3D::FIRE_VL()
{
    double 		dt_fire = 0.01;
	int 		Nt_fire = 10000000;
    double 		Fthresh = std::pow(10.0, -13);

    int 		n, nt = 0;
	vector<double> 	vx(N, 0.0), vy(N, 0.0), vz(N, 0.0);
	vector<double> 	wx(N, 0.0), wy(N, 0.0), wz(N, 0.0);
	double 		P, Vel_norm_t, Acc_norm_t, Vel_norm_r, Acc_norm_r;
	double 		dw, w_norm;
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
	double 		v_rsc_t, v_rsc_r;
	double 		q0, q1, q2, q3, q_norm;
	double 		dq0, dq1, dq2, dq3;
	double 		sin_dw;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;
	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>> 		VL(100 * Ns_tot, vector<int>(4, -1));
	int 		VL_counter;
	vector<double> 		x_bump_save(Ns_tot, 0.0), y_bump_save(Ns_tot, 0.0), z_bump_save(Ns_tot, 0.0);
	double 		r_cut = D_bump;

	// Verlet list initialization
    GetBumpPos();
	VerletList(r_cut, x_bump_save, y_bump_save, z_bump_save, VL, VL_counter, first_call);
    Force_VL(VL, VL_counter);

	if (MaxForce() < Fthresh){
		//printf("Eval: %.4e  Acc_max: %.4e\n", Eval, Acc_max);
		return;
	}

	first_call = 0;

	for (int nt = 1; nt < Nt_fire; nt++){
		P = 0.0;
		for (n = 0; n < N; n++)
			P += Fx[n] * vx[n] + Fy[n] * vy[n] + Fz[n] * vz[n] +
				 Tx[n] * wx[n] + Ty[n] * wy[n] + Tz[n] * wz[n];
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
            for (n = 0; n < N; n++){ // take systems half time step back and reset velocities
                x[n] -= dt_half * vx[n];
                y[n] -= dt_half * vy[n];
                z[n] -= dt_half * vz[n];

                w_norm = std::sqrt(wx[n] * wx[n] + wy[n] * wy[n] + wz[n] * wz[n]);
                if (w_norm < std::pow(10.0, -14)){
                    vx[n] = 0.0; vy[n] = 0.0; vz[n] = 0.0;
                    wx[n] = 0.0; wy[n] = 0.0; wz[n] = 0.0;
                    continue;
                }
                dw = -0.5 * dt_half * w_norm;
                sin_dw = std::sin(dw);
                dq0 = std::cos(dw);
                dq1 = sin_dw * wx[n] / w_norm;
                dq2 = sin_dw * wy[n] / w_norm;
                dq3 = sin_dw * wz[n] / w_norm;
                q0 = Q0[n]; q1 = Q1[n]; q2 = Q2[n]; q3 = Q3[n];
                Q0[n] = dq0 * q0 - dq1 * q1 - dq2 * q2 - dq3 * q3;
                Q1[n] = dq1 * q0 + dq0 * q1 - dq3 * q2 + dq2 * q3;
                Q2[n] = dq2 * q0 + dq3 * q1 + dq0 * q2 - dq1 * q3;
                Q3[n] = dq3 * q0 - dq2 * q1 + dq1 * q2 + dq0 * q3;
                q_norm = std::sqrt(Q0[n] * Q0[n] + Q1[n] * Q1[n] + Q2[n] * Q2[n] + Q3[n] * Q3[n]);
                Q0[n] /= q_norm; // make sure Q is always unit
                Q1[n] /= q_norm;
                Q2[n] /= q_norm;
                Q3[n] /= q_norm;

                vx[n] = 0.0; vy[n] = 0.0; vz[n] = 0.0;
                wx[n] = 0.0; wy[n] = 0.0; wz[n] = 0.0;
            }
		}

		// MD using Verlet method
		for (n = 0; n < N; n++){
			vx[n] += dt_half * Fx[n];
			vy[n] += dt_half * Fy[n];
			vz[n] += dt_half * Fz[n];
			wx[n] += dt_half * Tx[n];
			wy[n] += dt_half * Ty[n];
			wz[n] += dt_half * Tz[n];
		}
		Vel_norm_t = 0.0;
		Acc_norm_t = 0.0;
		Vel_norm_r = 0.0;
		Acc_norm_r = 0.0;
		for (n = 0; n < N; n++){
			Vel_norm_t += (vx[n] * vx[n] + vy[n] * vy[n] + vz[n] * vz[n]);
			Vel_norm_r += (wx[n] * wx[n] + wy[n] * wy[n] + wz[n] * wz[n]);
			Acc_norm_t += (Fx[n] * Fx[n] + Fy[n] * Fy[n] + Fz[n] * Fz[n]);
			Acc_norm_r += (Tx[n] * Tx[n] + Ty[n] * Ty[n] + Tz[n] * Tz[n]);
		}
		v_rsc_t = a_fire * std::sqrt(Vel_norm_t / Acc_norm_t);
		v_rsc_r = a_fire * std::sqrt(Vel_norm_r / Acc_norm_r);
		for (n = 0; n < N; n++){
			vx[n] = delta_a * vx[n] + v_rsc_t * Fx[n];
			vy[n] = delta_a * vy[n] + v_rsc_t * Fy[n];
			vz[n] = delta_a * vz[n] + v_rsc_t * Fz[n];
			wx[n] = delta_a * wx[n] + v_rsc_r * Tx[n];
			wy[n] = delta_a * wy[n] + v_rsc_r * Ty[n];
			wz[n] = delta_a * wz[n] + v_rsc_r * Tz[n];

			x[n] += dt * vx[n];
			y[n] += dt * vy[n];
			z[n] += dt * vz[n];

			w_norm = std::sqrt(wx[n] * wx[n] + wy[n] * wy[n] + wz[n] * wz[n]);
			if (w_norm < std::pow(10.0, -14))
				continue; // no need to rotate particles
			dw = 0.5 * dt * w_norm;
			sin_dw = std::sin(dw);
			dq0 = std::cos(dw);
			dq1 = sin_dw * wx[n] / w_norm;
			dq2 = sin_dw * wy[n] / w_norm;
			dq3 = sin_dw * wz[n] / w_norm;
			q0 = Q0[n]; q1 = Q1[n]; q2 = Q2[n]; q3 = Q3[n];
			Q0[n] = dq0 * q0 - dq1 * q1 - dq2 * q2 - dq3 * q3;
			Q1[n] = dq1 * q0 + dq0 * q1 - dq3 * q2 + dq2 * q3;
			Q2[n] = dq2 * q0 + dq3 * q1 + dq0 * q2 - dq1 * q3;
			Q3[n] = dq3 * q0 - dq2 * q1 + dq1 * q2 + dq0 * q3;
			q_norm = std::sqrt(Q0[n] * Q0[n] + Q1[n] * Q1[n] + Q2[n] * Q2[n] + Q3[n] * Q3[n]);
			Q0[n] /= q_norm; // make sure Q is always unit
			Q1[n] /= q_norm;
			Q2[n] /= q_norm;
			Q3[n] /= q_norm;
		}

		GetBumpPos();
		VerletList(r_cut, x_bump_save, y_bump_save, z_bump_save, VL, VL_counter, first_call);
		Force_VL(VL, VL_counter);

		for (n = 0; n < N; n++){
			vx[n] += dt_half * Fx[n];
			vy[n] += dt_half * Fy[n];
			vz[n] += dt_half * Fz[n];
			wx[n] += dt_half * Tx[n];
			wy[n] += dt_half * Ty[n];
			wz[n] += dt_half * Tz[n];
		}

		// if (nt % 10 == 0) printf("nt: %d VL_count: %d Max Force: %.5e\n", nt, VL_counter, MaxForce());

		if (MaxForce() < Fthresh){
			break;
		}
	}
}

void Bumpy_3D::GetBumpPos()
{
	// Bumpy center and orientation to bump center
	int 		n, ns, idx;
	double 		q0, q1, q2, q3;
	double 		q00, q11, q22, q33;
	double		q01, q02, q03;
	double 		q12, q13, q23;
	double 		R_mat_11, R_mat_12, R_mat_13;
	double 		R_mat_21, R_mat_22, R_mat_23;
	double 		R_mat_31, R_mat_32, R_mat_33;
	for (n = 0; n < N; n++){
		q0 = Q0[n];
		q1 = Q1[n];
		q2 = Q2[n];
		q3 = Q3[n];

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

		for (ns = idx_list[n]; ns < idx_list[n + 1]; ns++){
			idx = ns - idx_list[n];
			x_bump[ns] = x[n] + R_mat_11 * X_unit[idx] + R_mat_12 * Y_unit[idx] + R_mat_13 * Z_unit[idx];
			y_bump[ns] = y[n] + R_mat_21 * X_unit[idx] + R_mat_22 * Y_unit[idx] + R_mat_23 * Z_unit[idx];
			z_bump[ns] = z[n] + R_mat_31 * X_unit[idx] + R_mat_32 * Y_unit[idx] + R_mat_33 * Z_unit[idx];
		}
	}
}

void Bumpy_3D::VerletList(double r_cut,
                            vector<double> &x_bump_save,
                            vector<double> &y_bump_save,
                            vector<double> &z_bump_save,
                            vector<vector <int>> &VL,
                            int &VL_counter,
                            int first_call)
{
    double 		r_factor = 1.2;
	double 		Rc_cut = D_max + 1.01 * r_factor * D_bump;
	double 		Rc_cut_sq = Rc_cut * Rc_cut;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list_sq = r_factor * r_factor * r_cut_sq;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dz, dr;
	int 		n, m, ns, ms;

	if (first_call == 0){
		double dr_max = 0.0;
		for (ns = 0; ns < Ns_tot; ns++){
			dx = x_bump[ns] - x_bump_save[ns];
			dy = y_bump[ns] - y_bump_save[ns];
			dz = z_bump[ns] - z_bump_save[ns];
			dx = dx - std::round(dx / L) * L;
			// dy = dy - std::round(dy / L) * L;
			dz = dz - std::round(dz / L) * L;
			dr = dx * dx + dy * dy + dz * dz;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq)
	        return;
	}


	double 		x_cen_m, y_cen_m, z_cen_m;
	double 		x_cen_n, y_cen_n, z_cen_n;
	double 		dx_cen, dy_cen, dz_cen;
	double 		dx_cen_shift, dy_cen_shift, dz_cen_shift;
	double 		x1, y1, z1, x3, y3, z3;
	double 		pc1, pc2, dp_n, dp_m;

	VL_counter = 0;

	for (n = 0; n < N - 1; n++){
	    x_cen_n = x[n];
	    y_cen_n = y[n];
		z_cen_n = z[n];
	    for (m = n + 1; m < N; m++){
	        x_cen_m = x[m];
	        y_cen_m = y[m];
			z_cen_m = z[m];
	        dx_cen = x_cen_m - x_cen_n;
	        dy_cen = y_cen_m - y_cen_n;
			dz_cen = z_cen_m - z_cen_n;
	        dx_cen_shift = std::round(dx_cen / L) * L;
	        // dy_cen_shift = std::round(dy_cen / L) * L;
			dz_cen_shift = std::round(dz_cen / L) * L;
	        dx_cen = dx_cen - dx_cen_shift;
	        // dy_cen = dy_cen - dy_cen_shift;
			dz_cen = dz_cen - dz_cen_shift;
	        if (dx_cen * dx_cen + dy_cen * dy_cen + dz_cen * dz_cen > Rc_cut_sq){
	            continue;
	        }
	        
	        // eliminating vertices that cannot possibly be in contact
	        x_cen_m = x_cen_m - dx_cen_shift;
	        // y_cen_m = y_cen_m - dy_cen_shift;
			z_cen_m = z_cen_m - dz_cen_shift;
	        pc1 = -dx_cen * x_cen_n - dy_cen * y_cen_n - dz_cen * z_cen_n; // plane 1 equation constant C
	        pc2 = -dx_cen * x_cen_m - dy_cen * y_cen_m - dz_cen * z_cen_m; // plane 2 equation constant C
	        
	        for (ns = idx_list[n]; ns < idx_list[n + 1]; ns++){
	        	x1 = x_bump[ns]; y1 = y_bump[ns]; z1 = z_bump[ns];
	            dp_n = dx_cen * x1 + dy_cen * y1 + dz_cen * z1;
	            if ((dp_n + pc1) * (dp_n + pc2) > 0){
	                continue;
	            } // ns on the side that's away from possible cell-cell contact

	            for (ms = idx_list[m]; ms < idx_list[m + 1]; ms++){
	            	x3 = x_bump[ms] - dx_cen_shift;
	            	y3 = y_bump[ms] - dy_cen_shift;
					z3 = z_bump[ms] - dz_cen_shift;
	                dp_m = dx_cen * x3 + dy_cen * y3 + dz_cen * z3;
		            if ((dp_m + pc1) * (dp_m + pc2) > 0){
		                continue;
		            } // ms on the side that's away from possible cell-cell contact

	                dx = x3 - x1;
	                dy = y3 - y1;
					dz = z3 - z1;
	                dr = dx * dx + dy * dy + dz * dz;
	                if (dr < r_list_sq){
	                    VL[VL_counter][0] = ns;
						VL[VL_counter][1] = ms;
	                    VL[VL_counter][2] = n;
						VL[VL_counter][3] = m;
	                    VL_counter++;
						// printf("VL_count: %d n: %d m: %d ns: %d ms: %d dx_cen: %.5e dy_cen: %.5e dz_cen: %.5e dx: %.5e dy: %.5e dz: %.5e xn: %.5e xm: %.5e\n", VL_counter, n, m, ns, ms, dx_cen, dy_cen, dz_cen, dx, dy, dz, x[n], x[m]);
	                }
	            }
	        }
	    }
	}	

	for (ns = 0; ns < Ns_tot; ns++){
		x_bump_save[ns] = x_bump[ns];
		y_bump_save[ns] = y_bump[ns];
		z_bump_save[ns] = z_bump[ns];
	}
}

void Bumpy_3D::Force_VL(vector<vector <int>> &VL,
                        int VL_count)
{
	clearForces();
    InterForce_VL(VL, VL_count);
    WallForce();
}

void Bumpy_3D::InterForce_VL(vector<vector <int>> &VL,
                                int VL_count)
{
    int 		n, m, ns, ms;

    double 		dx, dy, dz, dnm;
	double 		F, dFx, dFy, dFz;
	double 		dxn, dyn, dzn;
	double 		dxm, dym, dzm;
    int 		vl_idx;

	for (vl_idx = 0; vl_idx < VL_count; vl_idx++){
		ns = VL[vl_idx][0];
		ms = VL[vl_idx][1];
        n = VL[vl_idx][2];
		m = VL[vl_idx][3];

		dx = x_bump[ms] - x_bump[ns];
		dx -= std::round(dx / L) * L;
		if (std::abs(dx) < D_bump){
			dy = y_bump[ms] - y_bump[ns];
			// dy -= std::round(dy / L) * L;
			if (std::abs(dy) < D_bump){
				dz = z_bump[ms] - z_bump[ns];
				dz -= std::round(dz / L) * L;
				if (std::abs(dz) < D_bump){
					dnm = std::sqrt(dx * dx + dy * dy + dz * dz);
					if (dnm < D_bump){
						F = K * (D_bump / dnm - 1.0);
						dFx = F * dx;
						dFy = F * dy;
						dFz = F * dz;

						Fx[n] -= dFx;
						Fy[n] -= dFy;
						Fz[n] -= dFz;
						Fx[m] += dFx;
						Fy[m] += dFy;
						Fz[m] += dFz;

                        dxn = x_bump[ns] - x[n];
                        dyn = y_bump[ns] - y[n];
                        dzn = z_bump[ns] - z[n];
                        dxm = x_bump[ms] - x[m];
                        dym = y_bump[ms] - y[m];
                        dzm = z_bump[ms] - z[m];

						Tx[n] -= (dyn * dFz - dzn * dFy);
						Ty[n] -= (dzn * dFx - dxn * dFz);
						Tz[n] -= (dxn * dFy - dyn * dFx);
						Tx[m] += (dym * dFz - dzm * dFy);
						Ty[m] += (dzm * dFx - dxm * dFz);
						Tz[m] += (dxm * dFy - dym * dFx);
					}
				}
			}
		}
	}
}

void Bumpy_3D::WallForce()
{
    int n, ns;
    double dFy;
    for (n = 0; n < N; n++){
        for (ns = idx_list[n]; ns < idx_list[n + 1]; ns++){
			if (y_bump[ns] < R_bump){
                dFy = Kw * (R_bump - y_bump[ns]);
                Fy[n] += dFy;
                Tx[n] -= (z_bump[ns] - z[n]) * dFy;
                Tz[n] += (x_bump[ns] - x[n]) * dFy;
			}
			else if (y_bump[ns] > L - R_bump){
                dFy = Kw * (L - R_bump - y_bump[ns]);
                Fy[n] += dFy;
                Tx[n] -= (z_bump[ns] - z[n]) * dFy;
                Tz[n] += (x_bump[ns] - x[n]) * dFy;
			}
		}
	}
}

double Bumpy_3D::Energy()
{
    int 		n, m, ns, ms;
    double 		dx, dy, dz, dnm;
	double 		E = 0.0;

	for (n = 0; n < N - 1; n++){
        for (m = n + 1; m < N; m++){
            for (ns = idx_list[n]; ns < idx_list[n + 1]; ns++){
                for (ms = idx_list[m]; ms < idx_list[m + 1]; ms++){
                    dx = x_bump[ms] - x_bump[ns];
                    dx -= std::round(dx / L) * L;
                    if (std::abs(dx) < D_bump){
                        dy = y_bump[ms] - y_bump[ns];
                        // dy -= std::round(dy / L) * L;
                        if (std::abs(dy) < D_bump){
                            dz = z_bump[ms] - z_bump[ns];
                            dz -= std::round(dz / L) * L;
                            if (std::abs(dz) < D_bump){
                                dnm = std::sqrt(dx * dx + dy * dy + dz * dz);
                                if (dnm < D_bump)
                                    E += 0.5 * K * (D_bump - dnm) * (D_bump - dnm);
                            }
                        }
                    }
                }
            }
        }
	}

    for (n = 0; n < N; n++){
        for (ns = idx_list[n]; ns < idx_list[n + 1]; ns++){
			if (y_bump[ns] < R_bump)
                E += 0.5 * Kw * (R_bump - y_bump[ns]) * (R_bump - y_bump[ns]);
			else if (y_bump[ns] > L - R_bump)
                E += 0.5 * Kw * (L - R_bump - y_bump[ns]) * (L - R_bump - y_bump[ns]);
		}
	}

    return E;
}

void Bumpy_3D::clearForces()
{
    for (int n = 0; n < N; n++){
        Fx[n] = 0.0;
        Fy[n] = 0.0;
        Fz[n] = 0.0;
        Tx[n] = 0.0;
        Ty[n] = 0.0;
        Tz[n] = 0.0;
    }
}

double Bumpy_3D::MaxForce()
{
	double Fmax = -1.0;
	double Fabs;

	for (int n = 0; n < N; n++){
		Fabs = std::abs(Fx[n]);
		Fmax = (Fmax > Fabs) ? Fmax : Fabs;
		Fabs = std::abs(Fy[n]);
		Fmax = (Fmax > Fabs) ? Fmax : Fabs;
		Fabs = std::abs(Fz[n]);
		Fmax = (Fmax > Fabs) ? Fmax : Fabs;
		Fabs = std::abs(Tx[n]);
		Fmax = (Fmax > Fabs) ? Fmax : Fabs;
		Fabs = std::abs(Ty[n]);
		Fmax = (Fmax > Fabs) ? Fmax : Fabs;
		Fabs = std::abs(Tz[n]);
		Fmax = (Fmax > Fabs) ? Fmax : Fabs;
	}

	return Fmax;
}

void Bumpy_3D::SaveConfig()
{
    FILE*   otptfl;
    if ((otptfl = fopen(rigid_name.c_str(), "w")) == NULL) {
		fprintf(stderr, "Unable to open %s.\n", rigid_name.c_str());
		exit(-1);
	}

    fprintf(otptfl, "%.32e\n", R);
	fprintf(otptfl, "%.32e\n", L);

	for (int n = 0; n < N; n++)
		fprintf(otptfl, "%.32e\n", x[n]);
    for (int n = 0; n < N; n++)
		fprintf(otptfl, "%.32e\n", y[n]);
	for (int n = 0; n < N; n++)
		fprintf(otptfl, "%.32e\n", z[n]);
    for (int n = 0; n < N; n++)
		fprintf(otptfl, "%.32e\n", Q0[n]);
    for (int n = 0; n < N; n++)
		fprintf(otptfl, "%.32e\n", Q1[n]);
	for (int n = 0; n < N; n++)
		fprintf(otptfl, "%.32e\n", Q2[n]);
    for (int n = 0; n < N; n++)
		fprintf(otptfl, "%.32e\n", Q3[n]);

	fclose(otptfl);
}

void Bumpy_3D::LoadConfig()
{
	std::ifstream bumpyfile(rigid_name.c_str());
	if (bumpyfile.good()){
		bumpyfile >> R;
		bumpyfile >> L;
		for (int i = 0; i < N; i++)
			bumpyfile >> x[i];
		for (int i = 0; i < N; i++)
			bumpyfile >> y[i];
		for (int i = 0; i < N; i++)
			bumpyfile >> z[i];
		for (int i = 0; i < N; i++)
			bumpyfile >> Q0[i];
		for (int i = 0; i < N; i++)
			bumpyfile >> Q1[i];
		for (int i = 0; i < N; i++)
			bumpyfile >> Q2[i];
		for (int i = 0; i < N; i++)
			bumpyfile >> Q3[i];
		bumpyfile.close();
	}
	else {
		bumpyfile.close();
		exit(-1);
	}
}

double Bumpy_3D::ReturnBumpPos(vector<double> &x_bump_return,
								vector<double> &y_bump_return,
								vector<double> &z_bump_return)
{
	LoadConfig();
	GetBumpPos();
	for (int ns = 0; ns < Ns_tot; ns++){
		x_bump_return[ns] = x_bump[ns];
		y_bump_return[ns] = y_bump[ns];
		z_bump_return[ns] = z_bump[ns];
	}

	return L;
}
