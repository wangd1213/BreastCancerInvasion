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

double BreastCancer_2D::LocalPF(double R_cut,
                                double R_cut_sq,
                                double R,
                                double D,
                                double R_sq,
                                double dnm,
                                double dnm_sq)
{
	if (dnm < R_cut - R) // particle m fully inside of the probe circle
		return PI * R_sq;
	else if (dnm < std::sqrt(R_cut_sq - R_sq)) { // particle m partially inside of the probe circle
		double d2 = (R_cut_sq - R_sq - dnm_sq) / 2.0 / dnm;
		double d1 = std::sqrt(R_sq - d2 * d2);
		double theta = std::acos(d2 / 0.5 / D);
		double theta_small = std::asin(d1 / R_cut);
		return theta_small * R_cut_sq + (d2 - dnm) * d1 + (PI - theta) * R_sq;
	}
	else {
		double d2 = (R_sq + dnm_sq - R_cut_sq) / 2.0 / dnm;
		double d1 = std::sqrt(R_sq - d2 * d2);
		double theta = std::acos(d2 / 0.5 / D);
		double theta_small = std::asin(d1 / R_cut);
		return theta_small * R_cut_sq - dnm * d1 + theta * R_sq;
	}	
}

void BreastCancer_2D::LocalPhi(vector<double> &xc,
                                vector<double> &yc,
                                double r_shell,
                                vector<vector <int>> &VL_phi_cc,
                                int VL_phi_count,
                                vector<double> &phi)
{
    int         nc, mc;
    double      D_phi, D_phi_n, D_phi_m, D_phi_n_sq, D_phi_m_sq;
    double      dx, dy, Dnm, dnm_sq, dnm;

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
        if (std::abs(dx) < D_phi){
            dy = yc[mc] - yc[nc];
            dy -= std::round(dy / Ly) * Ly;
            if (std::abs(dy) < D_phi){
                dnm_sq = dx * dx + dy * dy;
                dnm = std::sqrt(dnm_sq);
                if (dnm < D_phi){ // local packing fraction only
                    phi[nc] += LocalPF(D_phi_n, D_phi_n_sq, Rc[mc], Dc[mc], Rc_sq[mc], dnm, dnm_sq);
                    phi[mc] += LocalPF(D_phi_m, D_phi_m_sq, Rc[nc], Dc[nc], Rc_sq[nc], dnm, dnm_sq);
                }
            }
        }
	}

    double phi_mean_tumor = 0.85;
    double theta;
    for (nc = 0; nc < Nc; nc++){
        phi[nc] += PI * Rc_sq[nc];
        D_phi_n = Rc[nc] + r_shell;
        if (xc[nc] < D_phi_n){
            theta = std::acos(xc[nc] / D_phi_n);
            phi[nc] += phi_mean_tumor * (theta * D_phi_n - xc[nc] * std::sin(theta)) * D_phi_n;
        }
        phi[nc] = phi[nc] / PI / D_phi_n / D_phi_n;
    }
}

void BreastCancer_2D::LocalPhi_Update(vector<double> &xc,
                                        vector<double> &yc,
                                        double r_shell,
                                        vector<vector <int>> &VL_phi_cc,
                                        int VL_phi_count,
                                        vector<int> &is_dividing,
                                        vector<int> &is_not_dividing_ori,
                                        vector<double> &phi)
{
    int         nc, mc;
    double      D_phi, D_phi_n, D_phi_m, D_phi_n_sq, D_phi_m_sq;
    double      dx, dy, Dnm, dnm_sq, dnm;

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
        if (std::abs(dx) < D_phi){
            dy = yc[mc] - yc[nc];
            dy -= std::round(dy / Ly) * Ly;
            if (std::abs(dy) < D_phi){
                dnm_sq = dx * dx + dy * dy;
                dnm = std::sqrt(dnm_sq);
                if (dnm < D_phi){
                    if ((!is_not_dividing_ori[nc]) && (!is_dividing[nc]))
                        phi[nc] += LocalPF(D_phi_n, D_phi_n_sq, Rc[mc], Dc[mc], Rc_sq[mc], dnm, dnm_sq);
                    if ((!is_not_dividing_ori[mc]) && (!is_dividing[mc]))
                        phi[mc] += LocalPF(D_phi_m, D_phi_m_sq, Rc[nc], Dc[nc], Rc_sq[nc], dnm, dnm_sq);
                }
            }
        }
	}

    double phi_mean_tumor = 0.85;
    double theta;
    for (nc = 0; nc < Nc; nc++){
        if ((!is_not_dividing_ori[nc]) && (!is_dividing[nc])){
            phi[nc] += PI * Rc_sq[nc];
            D_phi_n = Rc[nc] + r_shell;
            if (xc[nc] < D_phi_n){
                theta = std::acos(xc[nc] / D_phi_n);
                phi[nc] += phi_mean_tumor * (theta * D_phi_n - xc[nc] * std::sin(theta)) * D_phi_n;
            }
            phi[nc] = phi[nc] / PI / D_phi_n / D_phi_n;
        }
    }
}
