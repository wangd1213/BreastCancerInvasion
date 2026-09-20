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

const double phi_mean_tumor = 0.65;

double BreastCancer_3D::LocalPF(double R_cut,
                                double R_cut_sq,
                                double R,
                                double D,
                                double R_sq,
                                double D_cube,
                                double dnm,
                                double dnm_sq)
{
	if (dnm < R_cut - R) // particle m fully inside of the probe circle
        return PI * D_cube / 6.0;
    else if (dnm < std::sqrt(R_cut_sq - R_sq)){ // particle m partially inside of the probe circle
        double d2 = (R_cut_sq - R_sq - dnm_sq) / 2.0 / dnm;
        double drh = R - d2;
        double drh_c = R_cut - d2 - dnm;
        return PI / 6.0 * (D_cube - 2.0 * (drh * drh * (D + d2) - drh_c * drh_c * (2.0 * R_cut + d2 + dnm)));
    }
    else {
        double d2 = (R_sq + dnm_sq - R_cut_sq) / 2.0 / dnm;
        double drh = R - d2;
        double drh_c = R_cut + d2 - dnm;
        return PI / 3.0 * (drh * drh * (D + d2) + drh_c * drh_c * (2.0 * R_cut - d2 + dnm));
    }
}

void BreastCancer_3D::LocalPhi(vector<double> &xc,
                                vector<double> &yc,
                                vector<double> &zc,
                                double r_shell,
                                vector<vector <int>> &VL_phi_cc,
                                int VL_phi_count,
                                vector<double> &phi)
{
    int         nc, mc;
    double      D_phi, D_phi_n, D_phi_m, D_phi_n_sq, D_phi_m_sq;
    double      dx, dy, dz, Dnm, dnm_sq, dnm;

    for (nc = 0; nc < Nc; nc++)
        phi[nc] = 0.0;

    for (int vl_idx = 0; vl_idx < VL_phi_count; vl_idx++){
		nc = VL_phi_cc[vl_idx][0];
		mc = VL_phi_cc[vl_idx][1];
        Dnm = Rc[nc] + Rc[mc];
        D_phi = Dnm + r_shell;
        D_phi_n = Rc[nc] + r_shell;
        D_phi_m = Rc[mc] + r_shell;
        D_phi_n_sq = D_phi_n * D_phi_n;
        D_phi_m_sq = D_phi_m * D_phi_m;
		dx = xc[mc] - xc[nc];
        dx -= std::round(dx / Lx) * Lx;
        if (std::abs(dx) < D_phi){
            dy = yc[mc] - yc[nc];
            // dy -= std::round(dy / Ly) * Ly;
            if (std::abs(dy) < D_phi){
                dz = zc[mc] - zc[nc];
                dz -= std::round(dz / Lz) * Lz;
                if (std::abs(dz) < D_phi){
                    dnm_sq = dx * dx + dy * dy + dz * dz;
                    dnm = std::sqrt(dnm_sq);
                    if (dnm < D_phi){ // local packing fraction only
                        phi[nc] += LocalPF(D_phi_n, D_phi_n_sq, Rc[mc], Dc[mc], Rc_sq[mc], Dc_cube[mc], dnm, dnm_sq);
                        phi[mc] += LocalPF(D_phi_m, D_phi_m_sq, Rc[nc], Dc[nc], Rc_sq[nc], Dc_cube[nc], dnm, dnm_sq);
                    }
                }
            }
        }
	}

    for (nc = 0; nc < Nc; nc++){
        phi[nc] += PI * Dc_cube[nc] / 6.0;
        D_phi_n = Rc[nc] + r_shell;
        if (yc[nc] < D_phi_n)
            phi[nc] += phi_mean_tumor * PI / 3.0 * (D_phi_n - yc[nc]) * (D_phi_n - yc[nc]) * (2.0 * D_phi_n + yc[nc]);
        phi[nc] = 0.75 * phi[nc] / PI / D_phi_n / D_phi_n / D_phi_n;
    }
}

void BreastCancer_3D::LocalPhi_Update(vector<double> &xc,
                                        vector<double> &yc,
                                        vector<double> &zc,
                                        double r_shell,
                                        vector<vector <int>> &VL_phi_cc,
                                        int VL_phi_count,
                                        vector<int> &is_dividing,
                                        vector<int> &is_not_dividing_ori,
                                        vector<double> &phi)
{
    int         nc, mc;
    double      D_phi, D_phi_n, D_phi_m, D_phi_n_sq, D_phi_m_sq;
    double      dx, dy, dz, Dnm, dnm_sq, dnm;

    for (nc = 0; nc < Nc; nc++)
        phi[nc] = 0.0;

    for (int vl_idx = 0; vl_idx < VL_phi_count; vl_idx++){
		nc = VL_phi_cc[vl_idx][0];
		mc = VL_phi_cc[vl_idx][1];
        if (((is_not_dividing_ori[nc]) && (is_not_dividing_ori[mc])) || ((is_dividing[nc]) && (is_dividing[mc])))
            continue;
        Dnm = Rc[nc] + Rc[mc];
        D_phi = Dnm + r_shell;
        D_phi_n = Rc[nc] + r_shell;
        D_phi_m = Rc[mc] + r_shell;
        D_phi_n_sq = D_phi_n * D_phi_n;
        D_phi_m_sq = D_phi_m * D_phi_m;
		dx = xc[mc] - xc[nc];
        dx -= std::round(dx / Lx) * Lx;
        if (std::abs(dx) < D_phi){
            dy = yc[mc] - yc[nc];
            // dy -= std::round(dy / Ly) * Ly;
            if (std::abs(dy) < D_phi){
                dz = zc[mc] - zc[nc];
                dz -= std::round(dz / Lz) * Lz;
                if (std::abs(dz) < D_phi){
                    dnm_sq = dx * dx + dy * dy + dz * dz;
                    dnm = std::sqrt(dnm_sq);
                    if (dnm < D_phi){
                        if ((!is_not_dividing_ori[nc]) && (!is_dividing[nc]))
                            phi[nc] += LocalPF(D_phi_n, D_phi_n_sq, Rc[mc], Dc[mc], Rc_sq[mc], Dc_cube[mc], dnm, dnm_sq);
                        if ((!is_not_dividing_ori[mc]) && (!is_dividing[mc]))
                            phi[mc] += LocalPF(D_phi_m, D_phi_m_sq, Rc[nc], Dc[nc], Rc_sq[nc], Dc_cube[mc], dnm, dnm_sq);
                    }
                }
            }
        }
	}
    
    for (nc = 0; nc < Nc; nc++){
        if ((!is_not_dividing_ori[nc]) && (!is_dividing[nc])){
            phi[nc] += PI * Dc_cube[nc] / 6.0;
            D_phi_n = Rc[nc] + r_shell;
            if (yc[nc] < D_phi_n)
                phi[nc] += phi_mean_tumor * PI / 3.0 * (D_phi_n - yc[nc]) * (D_phi_n - yc[nc]) * (2.0 * D_phi_n + yc[nc]);
            phi[nc] = 0.75 * phi[nc] / PI / D_phi_n / D_phi_n / D_phi_n;
        }
    }
}