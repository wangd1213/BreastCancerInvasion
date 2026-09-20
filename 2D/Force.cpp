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

void BreastCancer_2D::InterForce_Disk_VL(vector<vector <int>> &VL_cc,
                                            int VL_count)
{
    int 		n, m;

	// Disk-disk Force
    double 		Dnm, Dnm_a, Dnm_h, dnm, Dnm_sq, dd;
    double 		x1, dx, dy, dFx, dFy, F;
    int 		vl_idx;

	for (vl_idx = 0; vl_idx < VL_count; vl_idx++){
		n = VL_cc[vl_idx][0];
		m = VL_cc[vl_idx][1];

		Dnm = Rc[n] + Rc[m];
        Dnm_a = a_cut_c * Dnm;
        Dnm_h = a_cut_h * Dnm;
		Dnm_sq = Dnm * Dnm;
		
		dx = xc[m] - xc[n];
		if (std::abs(dx) < Dnm_a){
			dy = yc[m] - yc[n];
			dy -= std::round(dy / Ly) * Ly;
			if (std::abs(dy) < Dnm_a){
				dnm = std::sqrt(dx * dx + dy * dy);
				if (dnm < Dnm_a){
					// dd = 1.0 - dnm / Dnm;
					if (dnm <= Dnm)
						F = Kcc * (Dnm / dnm - 1.0) / Dnm_sq;
					else if ((dnm > Dnm) && (dnm < Dnm_h))
						F = Kcc_attr * (Dnm / dnm - 1.0) / Dnm_sq;
					else
						F = Kcc_attr * (1.0 - Dnm_a / dnm) / Dnm_sq;
					// Eval += 0.5 * dd * dd;

					dFx = F * dx;
					dFy = F * dy;
					Fx_c[n] -= dFx;
					Fy_c[n] -= dFy;
					Fx_c[m] += dFx; // 3rd law
					Fy_c[m]	+= dFy;
				}
			}
		}
	}
}

double BreastCancer_2D::WallForce_Disk()
{
    // Disk-disk Force
    int 		n;
    double 		Dnm;
    double 		dF;
    double      Fw_c = 0.0;

	for (n = 0; n < Nc; n++){
		Dnm = Rc[n];
		if (xc[n] < Dnm)
			Fx_c[n] += Kw_c * (Dnm - xc[n]) / Dnm / Dnm;
		else if (xc[n] > Lx - Dnm){
			dF = Kw_c * (Dnm - Lx + xc[n]) / Dnm / Dnm;
			Fx_c[n] -= dF;
            Fw_c += dF;
		}
	}

    return Fw_c;
}

void BreastCancer_2D::Force_Disk_VL(vector<vector <int>> &VL_cc,
                                    int VL_count)
{
    clearForces_Cancer();

	InterForce_Disk_VL(VL_cc, VL_count);
    Fw = WallForce_Disk();
}

void BreastCancer_2D::ShapeForce_DPM()
{
    int 		ns, ms, na, ma;
    int     	ifts, jfts;
    double  	lxi, lyi, lki, dli, dlj;
    double  	F;
	double 		L0_a, L0_sq, Kp_na;
	double 		Ns_a;
    double      A0_a, KA_a, Aa, dAa;
    double 		alphai, d_alpha_i, d_alpha_j, Kb_s;

	vector<double> 	lk_sq(Ns_tot, 0.0);
	vector<double> 	dxA(Ns_tot, 0.0);
	vector<double> 	dyA(Ns_tot, 0.0);
	vector<double> 	theta(Ns_tot, 0.0);
	vector<double>	d_alpha(Ns_tot, 0.0);
	vector<double>	cos_theta(Ns_tot, 0.0);
	vector<double> 	sin_theta(Ns_tot, 0.0);
	
	for (na = 0; na < Na; na++){
		L0_a = L0[na];
		L0_sq = L0_a * L0_a;
		Kp_na = Kp / double(Ns[na]) / L0_sq;

		for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
			ifts = ift[ns];
			lxi = xa[ifts] - xa[ns];
			lyi = ya[ifts] - ya[ns];
			lx[ns] = lxi;
			ly[ns] = lyi;
			lk_sq[ns] = lxi * lxi + lyi * lyi;
			lk[ns] = std::sqrt(lk_sq[ns]);
			// dli = lki / L0_c - 1.0;
			// Eval += Kp * dli * dli / double(Ns[nc]) / 2.0;
		}
		for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
			jfts = jft[ns];
	    	dlj = L0_a / lk[jfts] - 1.0;
	    	dli = 1.0 - L0_a / lk[ns];
	    	Fx_a[ns] += Kp_na * (dlj * lx[jfts] + dli * lx[ns]);
    		Fy_a[ns] += Kp_na * (dlj * ly[jfts] + dli * ly[ns]);
		}

        Aa = 0.0;
    	A0_a = A0[na];
        KA_a = ((double)Ns[na] / 20.0) * ((double)Ns[na] / 20.0) * KA;
        // printf("na: %d  Ka: %.5e\n", na, KA_a);
    	for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
    		ifts = ift[ns];
    		jfts = jft[ns];
    		Aa += (xa[ns] * ya[ifts] - ya[ns] * xa[ifts]);
    		dxA[ns] = (xa[ifts] - xa[jfts]) / A0_a;
    		dyA[ns] = (ya[jfts] - ya[ifts]) / A0_a;
    	}
    	Aa /= 2.0;
    	dAa = 0.5 * KA_a * (Aa / A0_a - 1.0);
    	// Eval += 0.5 * dAc * dAc / KA;
    	for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
    		Fx_a[ns] += dAa * dyA[ns];
    		Fy_a[ns] += dAa * dxA[ns];
    	}

        Kb_s = ((double)Ns[na] / 20.0) * ((double)Ns[na] / 20.0) * Kb_all[na];
        for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
    		theta[ns] = std::atan2(ly[ns], lx[ns]);
            cos_theta[ns] = lx[ns] / lk_sq[ns];
            sin_theta[ns] = ly[ns] / lk_sq[ns];
    	}
        for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
    		alphai = theta[ns] - theta[jft[ns]];
			d_alpha[ns] = alphai - std::round(alphai / 2.0 / PI) * 2.0 * PI - alpha_0[na];
    	}
        for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
            // printf("ns: %d  Kb: %.5e\n", ns, Kb_s);
            ifts = ift[ns];
            jfts = jft[ns];
            d_alpha_i = d_alpha[ns] - d_alpha[ifts];
            d_alpha_j = d_alpha[ns] - d_alpha[jfts];
            Fx_a[ns] -= Kb_s * (d_alpha_i * sin_theta[ns] + d_alpha_j * sin_theta[jfts]);
            Fy_a[ns] += Kb_s * (d_alpha_i * cos_theta[ns] + d_alpha_j * cos_theta[jfts]);
    	}
	}
}

void BreastCancer_2D::IntraForce_DPM()
{
    int 		ns, ms, na, ma;
    double  	F;

	// Cell-cell Force
    double 		Dnm, Dnm_sq;
    double 		dFx, dFy;
    double 		dx, dy;
    double 		dnm, dd;

	// intra-cell interaction
	for (na = 0; na < Na; na++){
	    Dnm = L0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (ns = idx_list[na]; ns < idx_list[na + 1] - 2; ns++){
	        for (ms = ns + 2; ms < idx_list[na + 1]; ms++){
	        	if ((ns == idx_list[na]) && (ms == idx_list[na + 1] - 1)){
		    		continue;
		    	}
	            dx = xa[ms] - xa[ns];
	            if (std::abs(dx) < Dnm){
	                dy = ya[ms] - ya[ns];
					if (std::abs(dy) < Dnm){
						dnm = std::sqrt(dx * dx + dy * dy);
						if(dnm < Dnm){
							F = Kaa * (1.0 - Dnm / dnm) / Dnm_sq;
							// dd = 1.0 - dnm / Dnm;
							// Eval += KC_half * dd * dd;  // cell-cell PE

							dFx = F * dx;
							dFy = F * dy;
							Fx_a[ns] += dFx;
							Fy_a[ns] += dFy;
							Fx_a[ms] -= dFx; // 3rd law
							Fy_a[ms] -= dFy;
						}
					}
	            }
	        }
	    }
	}
}

void BreastCancer_2D::InterForce_DPM_VL(vector<vector <int>> &VL_aa,
				                        int VL_count)
{
    // Cell-cell Force
    int 		ns, ms, na, ma;
    double 		Dnm, Dnm_sq;
    double 		F, dFx, dFy;
    double 		dx, dy;
    double 		dnm, dd;
	int 		vl_idx;

	for (vl_idx = 0; vl_idx < VL_count; vl_idx++){
		ns = VL_aa[vl_idx][0];
		ms = VL_aa[vl_idx][1];
		Dnm = R0[VL_aa[vl_idx][2]] + R0[VL_aa[vl_idx][3]];
		Dnm_sq = Dnm * Dnm;
		
		dx = xa[ms] - xa[ns];
		if (std::abs(dx) < Dnm){
			dy = ya[ms] - ya[ns];
			dy -= std::round(dy / Ly) * Ly;
			if (std::abs(dy) < Dnm){
				dnm = std::sqrt(dx * dx + dy * dy);
				if (dnm < Dnm){
					F = Kaa * (1.0 - Dnm / dnm) / Dnm_sq;
					dFx = F * dx;
					dFy = F * dy;
					Fx_a[ns] += dFx;
					Fy_a[ns] += dFy;
					Fx_a[ms] -= dFx;
					Fy_a[ms] -= dFy;
				}
			}
		}
	}
}

double BreastCancer_2D::WallForce_DPM()
{
    // adipocyte-wall force
    int     na, ns;
    double  Dnm, Dnm_sq;
    double  Fw_a = 0.0, dF;

	for (na = 0; na < Na; na++){
	    Dnm = R0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
			if (xa[ns] < Dnm)
				Fx_a[ns] += Kw_a * (Dnm - xa[ns]) / Dnm_sq;
			else if (xa[ns] > Lx - Dnm){
				dF = Kw_a * (Dnm - Lx + xa[ns]) / Dnm_sq;
				Fx_a[ns] -= dF;
				Fw_a += dF;
			}
	    }
	}

    return Fw_a;
}

void BreastCancer_2D::Force_DPM_VL(vector<vector <int>> &VL_aa,
				                    			int VL_count)
{
    clearForces_DPM();

	ShapeForce_DPM();
    IntraForce_DPM();
	InterForce_DPM_VL(VL_aa, VL_count);
    Fw = WallForce_DPM();
}

void BreastCancer_2D::InterForce_DPM_Disk_VL(vector<vector <int>> &VL_ac,
						                        int VL_count)
{
    int     ns, na, nc;
    double  Dnm, Dnm_sq, Dnm_a, Dnm_h;
    double  dx, dy, dnm;
    double  F, dFx, dFy;

	for (int vl_idx = 0; vl_idx < VL_count; vl_idx++){
		ns = VL_ac[vl_idx][0];
		nc = VL_ac[vl_idx][2];
		Dnm = R0[VL_ac[vl_idx][1]] + Rc[nc];
		Dnm_a = a_cut_ac * Dnm;
		Dnm_h = a_cut_ac_h * Dnm;
		Dnm_sq = Dnm * Dnm;
		
		dx = xc[nc] - xa[ns];
		if (std::abs(dx) < Dnm_a){
			dy = yc[nc] - ya[ns];
			dy -= std::round(dy / Ly) * Ly;
			if (std::abs(dy) < Dnm_a){
				dnm = std::sqrt(dx * dx + dy * dy);
				if (dnm < Dnm_a){
					if (dnm <= Dnm)
						F = Kac * (Dnm / dnm - 1.0) / Dnm_sq;
					else if ((dnm > Dnm) && (dnm < Dnm_h))
						F = Kac_attr * (Dnm / dnm - 1.0) / Dnm_sq;
					else
						F = Kac_attr * (1.0 - Dnm_a / dnm) / Dnm_sq;
					// Eval += 0.5 * dd * dd;

					dFx = F * dx;
					dFy = F * dy;
					Fx_a[ns] -= dFx;
					Fy_a[ns] -= dFy;
					Fx_c[nc] += dFx; // 3rd law
					Fy_c[nc] += dFy;
				}
			}
		}
	}
}

void BreastCancer_2D::InterForce_DPM_Disk_Smooth()
{
    int     ns, ms, na, nc;
	double 	C1, C2;
    double  Dnm, Dnm_sq, Dnm_a, Dnm_h;
    double  dx13, dy13, dnm, dx_line, dy_line, dnm_line;
	double 	d_alpha;
    double  F, dFx, dFy;
	double 	yc_box;
	vector<double> ya_cen(Na, 0.0);

	for (na = 0; na < Na; na++){
		for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
			ya_cen[na] += ya[ns];
		ya_cen[na] /= Ns[na];
	}

	for (nc = 0; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			Dnm = R0[na] + Rc[nc];
			Dnm_a = a_cut_ac * Dnm;
			Dnm_h = a_cut_ac_h * Dnm;
			Dnm_sq = Dnm * Dnm;
			yc_box = yc[nc] - std::round((yc[nc] - ya_cen[na]) / Ly) * Ly;
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				ms = ift[ns];
				dx13 = xc[nc] - xa[ns];
				dy13 = yc_box - ya[ns];

				C1 = lx[ns] * xa[ns] + ly[ns] * ya[ns];
				C2 = lx[ns] * xa[ms] + ly[ns] * ya[ms];
				if ((lx[ns] * xc[nc] + ly[ns] * yc_box - C1) * (lx[ns] * xc[nc] + ly[ns] * yc_box - C2) < 0.0){ // potential disk-line contact
					dnm_line = (ly[ns] * xc[nc] - lx[ns] * yc_box - xa[ns] * ya[ms] + xa[ms] * ya[ns]) / lk[ns];
					if ((dnm_line < Dnm_a) && (dnm_line > 0.0)){
						if (dnm_line <= Dnm)
							F = Kac * (dnm_line - Dnm) / Dnm_sq / lk[ns];
						else if ((dnm_line > Dnm) && (dnm_line < Dnm_h))
							F = Kac_attr * (dnm_line - Dnm) / Dnm_sq / lk[ns];
						else
							F = Kac_attr * (Dnm_a - dnm_line) / Dnm_sq / lk[ns];

						dx_line = lx[ns] * dnm_line / lk[ns];
						dy_line = ly[ns] * dnm_line / lk[ns];

						Fx_a[ns] -= F * (yc_box - ya[ms] + dx_line);
						Fy_a[ns] -= F * (xa[ms] - xc[nc] + dy_line);
						Fx_a[ms] += F * (dy13 + dx_line);
						Fy_a[ms] -= F * (dx13 - dy_line);
						Fx_c[nc] -= F * ly[ns]; // 3rd law
						Fy_c[nc] += F * lx[ns];
					}
				}
				else { // check for disk-disk contact
					dnm = std::sqrt(dx13 * dx13 + dy13 * dy13);
					if (dnm < Dnm_a){
						// check for possible contact with line (ns-1)-ns
						d_alpha = std::atan2(dy13, dx13) - std::atan2(-lx[jft[ns]], ly[jft[ns]]);
						d_alpha -= std::round(d_alpha / 2.0 / PI) * 2.0 * PI;
						if (d_alpha <= 0.0)
							continue;

						if (dnm <= Dnm)
							F = Kac * (Dnm / dnm - 1.0) / Dnm_sq;
						else if ((dnm > Dnm) && (dnm < Dnm_h))
							F = Kac_attr * (Dnm / dnm - 1.0) / Dnm_sq;
						else
							F = Kac_attr * (1.0 - Dnm_a / dnm) / Dnm_sq;

						dFx = F * dx13;
						dFy = F * dy13;
						Fx_a[ns] -= dFx;
						Fy_a[ns] -= dFy;
						Fx_c[nc] += dFx; // 3rd law
						Fy_c[nc] += dFy;
					}
				}
			}
		}
	}
}

void BreastCancer_2D::InterForce_DPM_Disk_Smooth_VL(vector<vector <int>> &VL_ac,
						                        	int VL_count)
{
    int     ns, ms, na, nc;
	double 	C1, C2;
    double  Dnm, Dnm_sq, Dnm_a, Dnm_h;
    double  dx13, dy13, dnm, dx_line, dy_line, dnm_line;
	double 	d_alpha;
    double  F, dFx, dFy;
	double 	yc_box;
	vector<double> 	ya_cen(Na, 0.0);

	for (na = 0; na < Na; na++){
		for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
			ya_cen[na] += ya[ns];
		ya_cen[na] /= Ns[na];
	}

	for (int vl_idx = 0; vl_idx < VL_count; vl_idx++){
		ns = VL_ac[vl_idx][0];
		na = VL_ac[vl_idx][1];
		nc = VL_ac[vl_idx][2];
		ms = ift[ns];
		Dnm = R0[na] + Rc[nc];
		Dnm_a = a_cut_ac * Dnm;
		Dnm_h = a_cut_ac_h * Dnm;
		Dnm_sq = Dnm * Dnm;
		yc_box = yc[nc] - std::round((yc[nc] - ya_cen[na]) / Ly) * Ly;

		dx13 = xc[nc] - xa[ns];
		dy13 = yc_box - ya[ns];

		C1 = lx[ns] * xa[ns] + ly[ns] * ya[ns];
		C2 = lx[ns] * xa[ms] + ly[ns] * ya[ms];
		if ((lx[ns] * xc[nc] + ly[ns] * yc_box - C1) * (lx[ns] * xc[nc] + ly[ns] * yc_box - C2) < 0.0){ // potential disk-line contact
			dnm_line = (ly[ns] * xc[nc] - lx[ns] * yc_box - xa[ns] * ya[ms] + xa[ms] * ya[ns]) / lk[ns];
			if ((dnm_line < Dnm_a) && (dnm_line > 0.0)){
				if (dnm_line <= Dnm)
					F = Kac * (dnm_line - Dnm) / Dnm_sq / lk[ns];
				else if ((dnm_line > Dnm) && (dnm_line < Dnm_h))
					F = Kac_attr * (dnm_line - Dnm) / Dnm_sq / lk[ns];
				else
					F = Kac_attr * (Dnm_a - dnm_line) / Dnm_sq / lk[ns];

				dx_line = lx[ns] * dnm_line / lk[ns];
				dy_line = ly[ns] * dnm_line / lk[ns];

				Fx_a[ns] -= F * (yc_box - ya[ms] + dx_line);
				Fy_a[ns] -= F * (xa[ms] - xc[nc] + dy_line);
				Fx_a[ms] += F * (dy13 + dx_line);
				Fy_a[ms] -= F * (dx13 - dy_line);
				Fx_c[nc] -= F * ly[ns]; // 3rd law
				Fy_c[nc] += F * lx[ns];
			}
		}
		else { // check for disk-disk contact
			dnm = std::sqrt(dx13 * dx13 + dy13 * dy13);
			if (dnm < Dnm_a){
				// check for possible contact with line (ns-1)-ns
				d_alpha = std::atan2(dy13, dx13) - std::atan2(-lx[jft[ns]], ly[jft[ns]]);
				d_alpha -= std::round(d_alpha / 2.0 / PI) * 2.0 * PI;
				if (d_alpha <= 0.0)
					continue;

				if (dnm <= Dnm)
					F = Kac * (Dnm / dnm - 1.0) / Dnm_sq;
				else if ((dnm > Dnm) && (dnm < Dnm_h))
					F = Kac_attr * (Dnm / dnm - 1.0) / Dnm_sq;
				else
					F = Kac_attr * (1.0 - Dnm_a / dnm) / Dnm_sq;

				dFx = F * dx13;
				dFy = F * dy13;
				Fx_a[ns] -= dFx;
				Fy_a[ns] -= dFy;
				Fx_c[nc] += dFx; // 3rd law
				Fy_c[nc] += dFy;
			}
		}
	}
}

void BreastCancer_2D::Force_DPM_Disk_VL(vector<vector <int>> &VL_aa,
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

void BreastCancer_2D::Force_DPM_Disk_Smooth_VL(vector<vector <int>> &VL_aa,
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

	InterForce_DPM_Disk_Smooth_VL(VL_ac, VL_count[2]);

	Fw += WallForce_DPM();
	Fw += WallForce_Disk();
}

void BreastCancer_2D::PinningForce_DPM_Center()
{
    double  xa_cen, ya_cen;
    double  dFx, dFy;
    double  Ns_a;
    int     ns;
    if (Kpin > 0.000001){
        for (int na = 0; na < Na; na++){
            Ns_a = (double)Ns[na];
            xa_cen = 0.0;
            ya_cen = 0.0;
            for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                xa_cen += xa[ns];
                ya_cen += ya[ns];
            }
            xa_cen /= Ns_a;
            ya_cen /= Ns_a;
            dFx = Kpin * (xa_cen_pin[na] - xa_cen) / Ns_a;
            dFy = Kpin * (ya_cen_pin[na] - ya_cen) / Ns_a;
            for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                Fx_a[ns] += dFx;
                Fy_a[ns] += dFy;
            }
        }
    }
}

void BreastCancer_2D::PinningForce_DPM_Center_LipidLoss(vector<int> &reach_min,
														vector<double> &rsc_all,
														double rsc_min)
{
    double  xa_cen, ya_cen;
    double  dFx, dFy, Kpin_eff;
    double  Ns_a;
    int     ns;
    if (Kpin > 0.000001){
        for (int na = 0; na < Na; na++){
			if (reach_min[na] == 0){
				Ns_a = (double)Ns[na];
				xa_cen = 0.0;
				ya_cen = 0.0;
				for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
					xa_cen += xa[ns];
					ya_cen += ya[ns];
				}
				xa_cen /= Ns_a;
				ya_cen /= Ns_a;
				Kpin_eff = Kpin * (rsc_all[na] - rsc_min) / (1.0 - rsc_min) / Ns_a;
				dFx = Kpin_eff * (xa_cen_pin[na] - xa_cen);
				dFy = Kpin_eff * (ya_cen_pin[na] - ya_cen);
				for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
					Fx_a[ns] += dFx;
					Fy_a[ns] += dFy;
				}
			}
        }
    }
}

void BreastCancer_2D::PinningForce_DPM_Vertex()
{
    int     ns;
    if (Kpin > 0.000001){
        for (int na = 0; na < Na; na++){
            for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
                Fx_a[ns] += Kpin * (xa_pin[ns] - xa[ns]);
                Fy_a[ns] += Kpin * (ya_pin[ns] - ya[ns]);
            }
        }
    }
}

void BreastCancer_2D::Force_DPM_Disk_Pin_VL(vector<vector <int>> &VL_aa,
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

void BreastCancer_2D::Force_DPM_Disk_PinCen_VL(vector<vector <int>> &VL_aa,
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

void BreastCancer_2D::Force_DPM_Disk_Smooth_PinCen_VL(vector<vector <int>> &VL_aa,
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

	InterForce_DPM_Disk_Smooth_VL(VL_ac, VL_count[2]);

	Fw += WallForce_DPM();
	Fw += WallForce_Disk();
}

void BreastCancer_2D::Force_DPM_Disk_Smooth_PinCen_LipidLoss_VL(vector<vector <int>> &VL_aa,
																vector<vector <int>> &VL_cc,
																vector<vector <int>> &VL_ac,
																vector<int> &VL_count,
																vector<int> &reach_min,
																vector<double> &rsc_all,
																double rsc_min)
{
    clearForces_DPM();
    clearForces_Cancer();
    Fw = 0.0;

	ShapeForce_DPM();
    IntraForce_DPM();
    PinningForce_DPM_Center_LipidLoss(reach_min, rsc_all, rsc_min);
    InterForce_DPM_VL(VL_aa, VL_count[0]);

	InterForce_Disk_VL(VL_cc, VL_count[1]);

	InterForce_DPM_Disk_Smooth_VL(VL_ac, VL_count[2]);

	Fw += WallForce_DPM();
	Fw += WallForce_Disk();
}

void BreastCancer_2D::clearForces_DPM()
{
    for (int ns = 0; ns < Ns_tot; ns++){
        Fx_a[ns] = 0.0;
        Fy_a[ns] = 0.0;
    }
}

void BreastCancer_2D::clearForces_Cancer()
{
    for (int nc = 0; nc < Nc; nc++){
        Fx_c[nc] = 0.0;
        Fy_c[nc] = 0.0;
    }
}

double BreastCancer_2D::Pressure_DPM()
{
	double 		Pwall = 0.0;
	double 		Dnm, Dnm_sq, dFx;

	// adipocyte-wall force
	for (int na = 0; na < Na; na++){
	    Dnm = R0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++)
	    	if (xa[ns] > Lx - Dnm)
				Pwall += Kw_a * (Dnm - Lx + xa[ns]) / Dnm_sq;
	}

	return Pwall / Ly;
}

double BreastCancer_2D::Pressure_Disk()
{
	// disk-wall force
    double  Fw = 0.0;
	for (int n = 0; n < Nc; n++)
		if (xc[n] > Lx - Rc[n])
            Fw += Kw_c * (Rc[n] - Lx + xc[n]) / Rc[n] / Rc[n];

    return Fw / Ly;
}

double BreastCancer_2D::Pressure_DPM_Disk()
{
	return Pressure_Disk() + Pressure_DPM();
}