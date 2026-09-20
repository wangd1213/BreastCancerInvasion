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

void BreastCancer_3D::InterForce_Disk()
{
    int 		n, m;

	// Disk-disk Force
    double 		Dnm, Dnm_a, Dnm_h, dnm, Dnm_sq, dd;
    double 		dx, dy, dz, dFx, dFy, dFz, F;

	for (n = 0; n < Nc - 1; n++){
		for (m = n + 1; m < Nc; m++){
			Dnm = Rc[n] + Rc[m];
			Dnm_a = a_cut_c * Dnm;
			Dnm_h = a_cut_h * Dnm;
			Dnm_sq = Dnm * Dnm;
			
			dx = xc[m] - xc[n];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < Dnm_a){
				dy = yc[m] - yc[n];
				// dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < Dnm_a){
					dz = zc[m] - zc[n];
					dz -= std::round(dz / Lz) * Lz;
					if (std::abs(dz) < Dnm_a){
						dnm = std::sqrt(dx * dx + dy * dy + dz * dz);
						if (dnm < Dnm_a){
							// dd = 1.0 - dnm / Dnm;
							if (dnm <= Dnm)
								F = Kcc * (Dnm / dnm - 1.0) / Dnm_sq;
							else if (dnm < Dnm_h)
								F = Kcc_attr * (Dnm / dnm - 1.0) / Dnm_sq;
							else
								F = Kcc_attr * (1.0 - Dnm_a / dnm) / Dnm_sq;
							// Eval += 0.5 * dd * dd;

							dFx = F * dx;
							dFy = F * dy;
							dFz = F * dz;
							Fx_c[n] -= dFx;
							Fy_c[n] -= dFy;
							Fz_c[n] -= dFz;
							Fx_c[m] += dFx; // 3rd law
							Fy_c[m]	+= dFy;
							Fz_c[m] += dFz;
						}
					}
				}
			}
		}
	}
}

void BreastCancer_3D::InterForce_Disk_VL(vector<vector <int>> &VL_cc,
                                            int VL_count)
{
    int 		n, m;

	// Disk-disk Force
    double 		Dnm, Dnm_a, Dnm_h, dnm, Dnm_sq, dd;
    double 		dx, dy, dz, dFx, dFy, dFz, F;
    int 		vl_idx;

	for (vl_idx = 0; vl_idx < VL_count; vl_idx++){
		n = VL_cc[vl_idx][0];
		m = VL_cc[vl_idx][1];

		Dnm = Rc[n] + Rc[m];
        Dnm_a = a_cut_c * Dnm;
        Dnm_h = a_cut_h * Dnm;
		Dnm_sq = Dnm * Dnm;
		
		dx = xc[m] - xc[n];
		dx -= std::round(dx / Lx) * Lx;
		if (std::abs(dx) < Dnm_a){
			dy = yc[m] - yc[n];
			// dy -= std::round(dy / Ly) * Ly;
			if (std::abs(dy) < Dnm_a){
				dz = zc[m] - zc[n];
				dz -= std::round(dz / Lz) * Lz;
				if (std::abs(dz) < Dnm_a){
					dnm = std::sqrt(dx * dx + dy * dy + dz * dz);
					if (dnm < Dnm_a){
						// dd = 1.0 - dnm / Dnm;
						if (dnm <= Dnm)
							F = Kcc * (Dnm / dnm - 1.0) / Dnm_sq;
						else if (dnm < Dnm_h)
							F = Kcc_attr * (Dnm / dnm - 1.0) / Dnm_sq;
						else
							F = Kcc_attr * (1.0 - Dnm_a / dnm) / Dnm_sq;
						// Eval += 0.5 * dd * dd;

						dFx = F * dx;
						dFy = F * dy;
						dFz = F * dz;
						Fx_c[n] -= dFx;
						Fy_c[n] -= dFy;
						Fz_c[n] -= dFz;
						Fx_c[m] += dFx; // 3rd law
						Fy_c[m]	+= dFy;
						Fz_c[m] += dFz;
					}
				}
			}
		}
	}
}

double BreastCancer_3D::WallForce_Disk()
{
    // Disk-disk Force
    int 		n;
    double 		Dnm;
    double 		dF;
    double      Fw_c = 0.0;

	for (n = 0; n < Nc; n++){
		Dnm = Rc[n];
		if (yc[n] < Dnm)
			Fy_c[n] += Kw_c * (Dnm - yc[n]) / Dnm / Dnm;
		else if (yc[n] > Ly - Dnm){
			dF = Kw_c * (Dnm - Ly + yc[n]) / Dnm / Dnm;
			Fy_c[n] -= dF;
            Fw_c += dF;
		}
	}

    return Fw_c;
}

void BreastCancer_3D::Force_Disk_VL(vector<vector <int>> &VL_cc,
                                    int VL_count)
{
    clearForces_Cancer();

	InterForce_Disk_VL(VL_cc, VL_count);
    Fw = WallForce_Disk();
}

void BreastCancer_3D::ShapeForce_DPM()
{
    int 		i, nf, ne;
    int 		v_idx_start;
	int 		ns1, ns2, ns3, ns4;
	int 		nf1, nf2;

	double 		x1, y1, z1;
	double 		x2, y2, z2;
	double 		x3, y3, z3;
	double 		x4, y4, z4;
	double 		dx12, dy12, dz12;
	double 		dx13, dy13, dz13;
	double 		dx14, dy14, dz14;
	double 		dx32, dy32, dz32;
	double 		dx34, dy34, dz34;
	double 		r13, r13_inv;
	double 		r_dot_1, r_dot_2, r_dot_3, r_dot_4;
	double 		FA2_x, FA2_y, FA2_z;
	double 		FA3_x, FA3_y, FA3_z;
	double 		dA, dV, V;
	double 		A_dot, A_sq, Ax, Ay, Az;
	double 		theta, dtheta; // face-face angle
	
	vector<double>	R_cross_32_x(Nf, 0.0), R_cross_32_y(Nf, 0.0), R_cross_32_z(Nf, 0.0);
	vector<double>	R_cross_13_x(Nf, 0.0), R_cross_13_y(Nf, 0.0), R_cross_13_z(Nf, 0.0);
	vector<double>	R_cross_21_x(Nf, 0.0), R_cross_21_y(Nf, 0.0), R_cross_21_z(Nf, 0.0);
	vector<double>	A(Nf, 0.0);
	vector<double>	A_vec_x(Nf, 0.0), A_vec_y(Nf, 0.0), A_vec_z(Nf, 0.0);
	vector<double>	A_rescale_x(Nf, 0.0), A_rescale_y(Nf, 0.0), A_rescale_z(Nf, 0.0);
	int 			face1[3], face2[3];

	for (int na = 0; na < Na; na++){
		v_idx_start = idx_list[na];
    	V = 0.0;
		// printf("na: %d\n", na);

    	// Area force
    	for (nf = 0; nf < Nf; nf++){
    		ns1 = f_unit[nf][0] + v_idx_start;
    		ns2 = f_unit[nf][1] + v_idx_start;
    		ns3 = f_unit[nf][2] + v_idx_start;

    		x1 = xa[ns1]; y1 = ya[ns1]; z1 = za[ns1];
    		x2 = xa[ns2]; y2 = ya[ns2]; z2 = za[ns2];
    		x3 = xa[ns3]; y3 = ya[ns3]; z3 = za[ns3];

    		// record Ri cross Rj for later volume force calculation
    		R_cross_32_x[nf] = y3 * z2 - z3 * y2;
    		R_cross_32_y[nf] = z3 * x2 - x3 * z2;
    		R_cross_32_z[nf] = x3 * y2 - y3 * x2;
    		R_cross_13_x[nf] = y1 * z3 - z1 * y3;
    		R_cross_13_y[nf] = z1 * x3 - x1 * z3;
    		R_cross_13_z[nf] = x1 * y3 - y1 * x3;
    		R_cross_21_x[nf] = y2 * z1 - z2 * y1;
    		R_cross_21_y[nf] = z2 * x1 - x2 * z1;
    		R_cross_21_z[nf] = x2 * y1 - y2 * x1;

    		dx12 = x2 - x1;
    		dy12 = y2 - y1;
    		dz12 = z2 - z1;
    		dx13 = x3 - x1;
    		dy13 = y3 - y1;
    		dz13 = z3 - z1;

    		Ax = (dy12 * dz13 - dz12 * dy13) * 0.5;
    		Ay = (dz12 * dx13 - dx12 * dz13) * 0.5;
    		Az = (dx12 * dy13 - dy12 * dx13) * 0.5;
    		A[nf] = std::sqrt(Ax * Ax + Ay * Ay + Az * Az);
			A_vec_x[nf] = Ax;
			A_vec_y[nf] = Ay;
			A_vec_z[nf] = Az;
			// printf("nf: %d A: %.5e\n", nf, A[nf]);

    		dA = (A[nf] - A0[na][nf]) * KA / A[nf];
    		FA2_x = dA * (Ay * dz13 - Az * dy13);
    		FA2_y = dA * (Az * dx13 - Ax * dz13);
    		FA2_z = dA * (Ax * dy13 - Ay * dx13);
    		FA3_x = dA * (dy12 * Az - dz12 * Ay);
    		FA3_y = dA * (dz12 * Ax - dx12 * Az);
    		FA3_z = dA * (dx12 * Ay - dy12 * Ax);

    		Fx_a[ns1] -= (FA2_x + FA3_x);
    		Fy_a[ns1] -= (FA2_y + FA3_y);
    		Fz_a[ns1] -= (FA2_z + FA3_z);
    		Fx_a[ns2] += FA2_x;
    		Fy_a[ns2] += FA2_y;
    		Fz_a[ns2] += FA2_z;
    		Fx_a[ns3] += FA3_x;
    		Fy_a[ns3] += FA3_y;
    		Fz_a[ns3] += FA3_z;

    		V += (Ax * x1 + Ay * y1 + Az * z1);
    	}
		V /= 3.0;
		// printf("V: %.5e\n", V);

    	// Volume force
    	dV = KV * (V - V0[na]) / 3.0;
    	for (nf = 0; nf < Nf; nf++){
    		ns1 = f_unit[nf][0] + v_idx_start;
    		ns2 = f_unit[nf][1] + v_idx_start;
    		ns3 = f_unit[nf][2] + v_idx_start;

    		Fx_a[ns1] += dV * R_cross_32_x[nf];
    		Fy_a[ns1] += dV * R_cross_32_y[nf];
    		Fz_a[ns1] += dV * R_cross_32_z[nf];
    		Fx_a[ns2] += dV * R_cross_13_x[nf];
    		Fy_a[ns2] += dV * R_cross_13_y[nf];
    		Fz_a[ns2] += dV * R_cross_13_z[nf];
    		Fx_a[ns3] += dV * R_cross_21_x[nf];
    		Fy_a[ns3] += dV * R_cross_21_y[nf];
    		Fz_a[ns3] += dV * R_cross_21_z[nf];
    	}

    	// Bending force
		if (Kb > std::pow(10.0, -8)){
			for (nf = 0; nf < Nf; nf++){
				A_sq = A[nf] * A[nf];
				A_rescale_x[nf] = A_vec_x[nf] / A_sq;
				A_rescale_y[nf] = A_vec_y[nf] / A_sq;
				A_rescale_z[nf] = A_vec_z[nf] / A_sq;
			}
			for (ne = 0; ne < Ne; ne++){
				nf1 = edgelist[ne][2];
				nf2 = edgelist[ne][3];
				for (i = 0; i < 3; i++){
					face1[i] = f_unit[nf1][i];
					face2[i] = f_unit[nf2][i];
				}
				ns1 = edgelist[ne][0] + v_idx_start;
				ns4 = face2[3 - edgelist[ne][5]] + v_idx_start;
				x1 = xa[ns1]; y1 = ya[ns1]; z1 = za[ns1];
				x4 = xa[ns4]; y4 = ya[ns4]; z4 = za[ns4];
				dx14 = x4 - x1;
				dy14 = y4 - y1;
				dz14 = z4 - z1;

				A_dot = (A_vec_x[nf1] * A_vec_x[nf2] + A_vec_y[nf1] * A_vec_y[nf2] + A_vec_z[nf1] * A_vec_z[nf2]) / A[nf1] / A[nf2];
				if (A_dot >= 1.0)
					theta = 0.0;
				else
					theta = std::acos(A_dot);
				if ((A_vec_x[nf1] * dx14 + A_vec_y[nf1] * dy14 + A_vec_z[nf1] * dz14) < 0.0)
					theta = -theta;
				dtheta = theta - theta0[ne];
				// printf("ne: %d  theta: %.5e\n", theta);

				ns2 = face1[3 - edgelist[ne][4]] + v_idx_start;
				ns3 = edgelist[ne][1] + v_idx_start;

				x2 = xa[ns2]; y2 = ya[ns2]; z2 = za[ns2];
				x3 = xa[ns3]; y3 = ya[ns3]; z3 = za[ns3];

				dx12 = x2 - x1;
				dy12 = y2 - y1;
				dz12 = z2 - z1;
				dx13 = x3 - x1;
				dy13 = y3 - y1;
				dz13 = z3 - z1;
				dx32 = x2 - x3;
				dy32 = y2 - y3;
				dz32 = z2 - z3;
				dx34 = x4 - x3;
				dy34 = y4 - y3;
				dz34 = z4 - z3;

				r13 = std::sqrt(dx13 * dx13 + dy13 * dy13 + dz13 * dz13);
				r13_inv = Kb * dtheta / r13;
				r13 *= dtheta * Kb;
				r_dot_1 = dx13 * dx32 + dy13 * dy32 + dz13 * dz32; // R13 * R32
				r_dot_2 = dx13 * dx34 + dy13 * dy34 + dz13 * dz34; // R13 * R34
				r_dot_3 = dx13 * dx12 + dy13 * dy12 + dz13 * dz12; // R13 * R12
				r_dot_4 = dx13 * dx14 + dy13 * dy14 + dz13 * dz14; // R13 * R14

				Fx_a[ns1] -= (r_dot_1 * A_rescale_x[nf1] + r_dot_2 * A_rescale_x[nf2]) * r13_inv;
				Fy_a[ns1] -= (r_dot_1 * A_rescale_y[nf1] + r_dot_2 * A_rescale_y[nf2]) * r13_inv;
				Fz_a[ns1] -= (r_dot_1 * A_rescale_z[nf1] + r_dot_2 * A_rescale_z[nf2]) * r13_inv;
				Fx_a[ns2] -= A_rescale_x[nf1] * r13;
				Fy_a[ns2] -= A_rescale_y[nf1] * r13;
				Fz_a[ns2] -= A_rescale_z[nf1] * r13;
				Fx_a[ns3] += (r_dot_3 * A_rescale_x[nf1] + r_dot_4 * A_rescale_x[nf2]) * r13_inv;
				Fy_a[ns3] += (r_dot_3 * A_rescale_y[nf1] + r_dot_4 * A_rescale_y[nf2]) * r13_inv;
				Fz_a[ns3] += (r_dot_3 * A_rescale_z[nf1] + r_dot_4 * A_rescale_z[nf2]) * r13_inv;
				Fx_a[ns4] -= A_rescale_x[nf2] * r13;
				Fy_a[ns4] -= A_rescale_y[nf2] * r13;
				Fz_a[ns4] -= A_rescale_z[nf2] * r13;
			}
		}
    }
}

void BreastCancer_3D::IntraForce_DPM()
{
    int 		i, ns, ms, na;
    double 		dF, dFx, dFy, dFz;
    double 		dx, dy, dz;
    double 		dnm, dd;

	// intra-cell interaction
	for (na = 0; na < Na; na++){
		for (i = 0; i < Nnn; i++){
			ns = nnlist[i][0] + idx_list[na];
			ms = nnlist[i][1] + idx_list[na];
			dx = xa[ms] - xa[ns];
			if (std::abs(dx) < D0_intra[na]){
				dy = ya[ms] - ya[ns];
				if (std::abs(dy) < D0_intra[na]){
					dz = za[ms] - za[ns];
					if (std::abs(dz) < D0_intra[na]){
						dnm = std::sqrt(dx * dx + dy * dy + dz * dz);
						if (dnm < D0_intra[na]){
							dF = Kaa * (1.0 - D0_intra[na] / dnm) / D0_intra_sq[na];
							dFx = dF * dx;
							dFy = dF * dy;
							dFz = dF * dz;
							Fx_a[ns] += dFx;
							Fy_a[ns] += dFy;
							Fz_a[ns] += dFz;
							Fx_a[ms] -= dFx; // 3rd law
							Fy_a[ms] -= dFy;
							Fz_a[ms] -= dFz;
						}
					}
				}
			}
	    }
	}
}

void BreastCancer_3D::InterForce_DPM()
{
    // Cell-cell Force
    int 		ns, ms, na, ma;
    double 		Dnm, Dnm_sq;
    double 		F, dFx, dFy, dFz;
    double 		dx, dy, dz;
    double 		dnm, dd;

	for (na = 0; na < Na - 1; na++){
		for (ma = na + 1; ma < Na; ma++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				for (ms = idx_list[ma]; ms < idx_list[ma + 1]; ms++){
					Dnm = R0[na] + R0[ma];
					Dnm_sq = Dnm * Dnm;
					
					dx = xa[ms] - xa[ns];
					dx -= std::round(dx / Lx) * Lx;
					if (std::abs(dx) < Dnm){
						dy = ya[ms] - ya[ns];
						// dy -= std::round(dy / Ly) * Ly;
						if (std::abs(dy) < Dnm){
							dz = za[ms] - za[ns];
							dz -= std::round(dz / Lz) * Lz;
							if (std::abs(dz) < Dnm){
								dnm = std::sqrt(dx * dx + dy * dy + dz * dz);
								if (dnm < Dnm){
									F = Kaa * (1.0 - Dnm / dnm) / Dnm_sq;
									dFx = F * dx;
									dFy = F * dy;
									dFz = F * dz;
									Fx_a[ns] += dFx;
									Fy_a[ns] += dFy;
									Fz_a[ns] += dFz;
									Fx_a[ms] -= dFx;
									Fy_a[ms] -= dFy;
									Fz_a[ms] -= dFz;
								}
							}
						}
					}
				}
			}
		}
	}
}

void BreastCancer_3D::InterForce_DPM_VL(vector<vector <int>> &VL_aa,
				                        int VL_count)
{
    // Cell-cell Force
    int 		ns, ms, na, ma;
    double 		Dnm, Dnm_sq;
    double 		F, dFx, dFy, dFz;
    double 		dx, dy, dz;
    double 		dnm, dd;
	int 		vl_idx;

	for (vl_idx = 0; vl_idx < VL_count; vl_idx++){
		ns = VL_aa[vl_idx][0];
		ms = VL_aa[vl_idx][1];
		Dnm = R0[VL_aa[vl_idx][2]] + R0[VL_aa[vl_idx][3]];
		Dnm_sq = Dnm * Dnm;
		
		dx = xa[ms] - xa[ns];
		dx -= std::round(dx / Lx) * Lx;
		if (std::abs(dx) < Dnm){
			dy = ya[ms] - ya[ns];
			// dy -= std::round(dy / Ly) * Ly;
			if (std::abs(dy) < Dnm){
				dz = za[ms] - za[ns];
				dz -= std::round(dz / Lz) * Lz;
				if (std::abs(dz) < Dnm){
					dnm = std::sqrt(dx * dx + dy * dy + dz * dz);
					if (dnm < Dnm){
						F = Kaa * (1.0 - Dnm / dnm) / Dnm_sq;
						dFx = F * dx;
						dFy = F * dy;
						dFz = F * dz;
						Fx_a[ns] += dFx;
						Fy_a[ns] += dFy;
						Fz_a[ns] += dFz;
						Fx_a[ms] -= dFx;
						Fy_a[ms] -= dFy;
						Fz_a[ms] -= dFz;
					}
				}
			}
		}
	}
}

double BreastCancer_3D::WallForce_DPM()
{
    // adipocyte-wall force
    int     na, ns;
    double  Dnm, Dnm_sq;
    double  Fw_a = 0.0, dF;

	for (na = 0; na < Na; na++){
	    Dnm = R0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
			if (ya[ns] < Dnm)
				Fy_a[ns] += Kw_a * (Dnm - ya[ns]) / Dnm_sq;
			else if (ya[ns] > Ly - Dnm){
				dF = Kw_a * (Dnm - Ly + ya[ns]) / Dnm_sq;
				Fy_a[ns] -= dF;
				Fw_a += dF;
			}
	    }
	}

    return Fw_a;
}

void BreastCancer_3D::Force_DPM()
{
    clearForces_DPM();

	ShapeForce_DPM();
    IntraForce_DPM();
	InterForce_DPM();
    Fw = WallForce_DPM();
}

void BreastCancer_3D::Force_DPM_VL(vector<vector <int>> &VL_aa,
									int VL_count)
{
    clearForces_DPM();

	ShapeForce_DPM();
    IntraForce_DPM();
	InterForce_DPM_VL(VL_aa, VL_count);
    Fw = WallForce_DPM();
}

void BreastCancer_3D::InterForce_DPM_Disk()
{
    int     ns, na, nc;
    double  Dnm, Dnm_sq, Dnm_a, Dnm_h;
    double  dx, dy, dz, dnm;
    double  F, dFx, dFy, dFz;

	for (nc = 0; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			Dnm = R0[na] + Rc[nc];
			Dnm_a = a_cut_ac * Dnm;
			Dnm_h = a_cut_ac_h * Dnm;
			Dnm_sq = Dnm * Dnm;
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xc[nc] - xa[ns];
				dx -= std::round(dx / Lx) * Lx;
				if (std::abs(dx) < Dnm_a){
					dy = yc[nc] - ya[ns];
					// dy -= std::round(dy / Ly) * Ly;
					if (std::abs(dy) < Dnm_a){
						dz = zc[nc] - za[ns];
						dz -= std::round(dz / Lz) * Lz;
						if (std::abs(dz) < Dnm_a){
							dnm = std::sqrt(dx * dx + dy * dy + dz * dz);
							if (dnm < Dnm_a){
								if (dnm <= Dnm)
									F = Kac * (Dnm / dnm - 1.0) / Dnm_sq;
								else if (dnm < Dnm_h)
									F = Kac_attr * (Dnm / dnm - 1.0) / Dnm_sq;
								else
									F = Kac_attr * (1.0 - Dnm_a / dnm) / Dnm_sq;
								// Eval += 0.5 * dd * dd;

								dFx = F * dx;
								dFy = F * dy;
								dFz = F * dz;
								Fx_a[ns] -= dFx;
								Fy_a[ns] -= dFy;
								Fz_a[ns] -= dFz;
								Fx_c[nc] += dFx; // 3rd law
								Fy_c[nc] += dFy;
								Fz_c[nc] += dFz;
							}
						}
					}
				}
			}
		}
	}
}

void BreastCancer_3D::InterForce_DPM_Disk_VL(vector<vector <int>> &VL_ac,
						                        int VL_count)
{
    int     ns, na, nc;
    double  Dnm, Dnm_sq, Dnm_a, Dnm_h;
    double  dx, dy, dz, dnm;
    double  F, dFx, dFy, dFz;

	for (int vl_idx = 0; vl_idx < VL_count; vl_idx++){
		ns = VL_ac[vl_idx][0];
		nc = VL_ac[vl_idx][2];
		Dnm = R0[VL_ac[vl_idx][1]] + Rc[nc];
		Dnm_a = a_cut_ac * Dnm;
		Dnm_h = a_cut_ac_h * Dnm;
		Dnm_sq = Dnm * Dnm;
		
		dx = xc[nc] - xa[ns];
		dx -= std::round(dx / Lx) * Lx;
		if (std::abs(dx) < Dnm_a){
			dy = yc[nc] - ya[ns];
			// dy -= std::round(dy / Ly) * Ly;
			if (std::abs(dy) < Dnm_a){
				dz = zc[nc] - za[ns];
				dz -= std::round(dz / Lz) * Lz;
				if (std::abs(dz) < Dnm_a){
					dnm = std::sqrt(dx * dx + dy * dy + dz * dz);
					if (dnm < Dnm_a){
						if (dnm <= Dnm)
							F = Kac * (Dnm / dnm - 1.0) / Dnm_sq;
						else if (dnm < Dnm_h)
							F = Kac_attr * (Dnm / dnm - 1.0) / Dnm_sq;
						else
							F = Kac_attr * (1.0 - Dnm_a / dnm) / Dnm_sq;
						// Eval += 0.5 * dd * dd;

						dFx = F * dx;
						dFy = F * dy;
						dFz = F * dz;
						Fx_a[ns] -= dFx;
						Fy_a[ns] -= dFy;
						Fz_a[ns] -= dFz;
						Fx_c[nc] += dFx; // 3rd law
						Fy_c[nc] += dFy;
						Fz_c[nc] += dFz;
					}
				}
			}
		}
	}
}

void BreastCancer_3D::Force_DPM_FixDisk_VL(vector<vector <int>> &VL_aa,
											vector<vector <int>> &VL_ac,
											vector<int> &VL_count)
{
    clearForces_DPM();
    clearForces_Cancer();
    Fw = 0.0;

	ShapeForce_DPM();
    IntraForce_DPM();
    InterForce_DPM_VL(VL_aa, VL_count[0]);

	InterForce_DPM_Disk_VL(VL_ac, VL_count[1]);

	Fw += WallForce_DPM();
	Fw += WallForce_Disk();
}

void BreastCancer_3D::Force_Disk_FixDPM_VL(vector<vector <int>> &VL_cc,
											vector<vector <int>> &VL_ac,
											vector<int> &VL_count)
{
    clearForces_DPM();
    clearForces_Cancer();
    Fw = 0.0;

	InterForce_Disk_VL(VL_cc, VL_count[0]);

	InterForce_DPM_Disk_VL(VL_ac, VL_count[1]);

	Fw += WallForce_DPM();
	Fw += WallForce_Disk();
}

void BreastCancer_3D::Force_DPM_Disk_VL(vector<vector <int>> &VL_aa,
										vector<vector <int>> &VL_cc,
										vector<vector <int>> &VL_ac,
										vector<int> &VL_count)
{
    clearForces_DPM();
    clearForces_Cancer();
    Fw = 0.0;

	ShapeForce_DPM();
    IntraForce_DPM();
    InterForce_DPM_VL(VL_aa, VL_count[0]);

	InterForce_Disk_VL(VL_cc, VL_count[1]);

	InterForce_DPM_Disk_VL(VL_ac, VL_count[2]);

	Fw += WallForce_DPM();
	Fw += WallForce_Disk();
}

void BreastCancer_3D::PinningForce_DPM_Center()
{
    double  xa_cen, ya_cen, za_cen;
    double  dFx, dFy, dFz;
    int     ns;
    if (Kpin > 0.000001){
        for (int na = 0; na < Na; na++){
            xa_cen = 0.0;
            ya_cen = 0.0;
			za_cen = 0.0;
            for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                xa_cen += xa[ns];
                ya_cen += ya[ns];
				za_cen += za[ns];
            }
            xa_cen /= Ns;
            ya_cen /= Ns;
			za_cen /= Ns;
            dFx = Kpin * (xa_cen_pin[na] - xa_cen) / Ns;
            dFy = Kpin * (ya_cen_pin[na] - ya_cen) / Ns;
			dFz = Kpin * (za_cen_pin[na] - za_cen) / Ns;
            for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                Fx_a[ns] += dFx;
                Fy_a[ns] += dFy;
				Fz_a[ns] += dFz;
            }
        }
    }
}

void BreastCancer_3D::PinningForce_DPM_Center_LipidLoss(vector<int> &reach_min,
														vector<double> &rsc_all,
														double rsc_min)
{
    double  xa_cen, ya_cen, za_cen;
    double  dFx, dFy, dFz, Kpin_eff;
    int     ns;
    if (Kpin > 0.000001){
        for (int na = 0; na < Na; na++){
			if (reach_min[na] == 0){
				xa_cen = 0.0;
				ya_cen = 0.0;
				za_cen = 0.0;
				for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
					xa_cen += xa[ns];
					ya_cen += ya[ns];
					za_cen += za[ns];
				}
				xa_cen /= Ns;
				ya_cen /= Ns;
				za_cen /= Ns;
				Kpin_eff = Kpin * (rsc_all[na] - rsc_min) / (1.0 - rsc_min) / Ns;
				dFx = Kpin_eff * (xa_cen_pin[na] - xa_cen);
				dFy = Kpin_eff * (ya_cen_pin[na] - ya_cen);
				dFz = Kpin_eff * (za_cen_pin[na] - za_cen);
				for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
					Fx_a[ns] += dFx;
					Fy_a[ns] += dFy;
					Fz_a[ns] += dFz;
				}
			}
        }
    }
}

void BreastCancer_3D::PinningForce_DPM_Vertex()
{
    int     ns;
    if (Kpin > 0.000001){
        for (int na = 0; na < Na; na++){
            for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                Fx_a[ns] += Kpin * (xa_pin[ns] - xa[ns]);
                Fy_a[ns] += Kpin * (ya_pin[ns] - ya[ns]);
				Fz_a[ns] += Kpin * (za_pin[ns] - za[ns]);
            }
        }
    }
}

void BreastCancer_3D::Force_DPM_Disk_Pin_VL(vector<vector <int>> &VL_aa,
											vector<vector <int>> &VL_cc,
											vector<vector <int>> &VL_ac,
											vector<int> &VL_count)
{
    clearForces_DPM();
    clearForces_Cancer();
    Fw = 0.0;

	ShapeForce_DPM();
    IntraForce_DPM();
    PinningForce_DPM_Vertex();
    InterForce_DPM_VL(VL_aa, VL_count[0]);

	InterForce_Disk_VL(VL_cc, VL_count[1]);

	InterForce_DPM_Disk_VL(VL_ac, VL_count[2]);

	Fw += WallForce_DPM();
	Fw += WallForce_Disk();
}

void BreastCancer_3D::Force_DPM_Disk_PinCen_VL(vector<vector <int>> &VL_aa,
												vector<vector <int>> &VL_cc,
												vector<vector <int>> &VL_ac,
												vector<int> &VL_count)
{
    clearForces_DPM();
    clearForces_Cancer();
    Fw = 0.0;

	ShapeForce_DPM();
    IntraForce_DPM();
    PinningForce_DPM_Center();
    InterForce_DPM_VL(VL_aa, VL_count[0]);

	InterForce_Disk_VL(VL_cc, VL_count[1]);

	InterForce_DPM_Disk_VL(VL_ac, VL_count[2]);

	Fw += WallForce_DPM();
	Fw += WallForce_Disk();
}

void BreastCancer_3D::clearForces_DPM()
{
    for (int ns = 0; ns < Ns_tot; ns++){
        Fx_a[ns] = 0.0;
        Fy_a[ns] = 0.0;
		Fz_a[ns] = 0.0;
    }
}

void BreastCancer_3D::clearForces_Cancer()
{
    for (int nc = 0; nc < Nc; nc++){
        Fx_c[nc] = 0.0;
        Fy_c[nc] = 0.0;
		Fz_c[nc] = 0.0;
    }
}

double BreastCancer_3D::Pressure_DPM()
{
	double 		Pwall = 0.0;
	double 		Dnm, Dnm_sq, dFy;

	// adipocyte-wall force
	for (int na = 0; na < Na; na++){
	    Dnm = R0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
	    	if (ya[ns] > Ly - Dnm)
				Pwall += Kw_a * (Dnm - Ly + ya[ns]) / Dnm_sq;
	}

	return Pwall / Lx / Lz;
}

double BreastCancer_3D::Pressure_Disk()
{
	// disk-wall force
    double  Fw = 0.0;
	for (int n = 0; n < Nc; n++)
		if (yc[n] > Ly - Rc[n])
            Fw += Kw_c * (Rc[n] - Ly + yc[n]) / Rc[n] / Rc[n];

    return Fw / Lx / Lz;
}

double BreastCancer_3D::Pressure_DPM_Disk()
{
	return Pressure_Disk() + Pressure_DPM();
}