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

void BreastCancer_2D::VerletList_Disk(double r_cut,
			                            vector<double> &xc_save,
                                        vector<double> &yc_save,
			                            vector<vector <int>> &VL_cc,
                                        int &VL_count,
                                        int first_call)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list = r_factor * r_cut;
	double 		r_list_sq = r_list * r_list;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dr;
	int 		n, ny, m;

	if (first_call == 0){
		double dr_max = 0.0;
		for (n = 0; n < Nc; n++){
			dx = xc[n] - xc_save[n];
			dy = yc[n] - yc_save[n];
			dy -= std::round(dy / Ly) * Ly;
			dr = dx * dx + dy * dy;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq)
	        return;
	}

	VL_count = 0;

	for (n = 0; n < Nc - 1; n++){
	    for (m = n + 1; m < Nc; m++){
			dx = xc[m] - xc[n];
			if (std::abs(dx) < r_list){
				dy = yc[m] - yc[n];
				dy -= std::round(dy / Ly) * Ly;
				if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
					VL_cc[VL_count][0] = n;
					VL_cc[VL_count][1] = m;
					VL_count++;
				}
			}
		}
	}	

	for (n = 0; n < Nc; n++){
		xc_save[n] = xc[n];
        yc_save[n] = yc[n];
    }
}

void BreastCancer_2D::VerletList_DPM(double r_cut,
										vector<double> &xa_save,
										vector<double> &ya_save,
										vector<vector <int>> &VL_aa,
										int &VL_count,
										int first_call)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list = r_factor * r_cut;
	double 		r_list_sq = r_list * r_list;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dr;
	int 		i, j;
	int 		na, ma, ns, ms, nc, mc;

	if (first_call == 0){
		double dr_max = 0.0;
		for (i = 0; i < Ns_tot; i++){
			dx = xa[i] - xa_save[i];
			dy = ya[i] - ya_save[i];
			dy -= std::round(dy / Ly) * Ly;
			dr = dx * dx + dy * dy;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return;
	    }
	}

	vector<double>	x_cen(Na, 0.0);
	vector<double> 	y_cen(Na, 0.0);
	for (na = 0; na < Na; na++){
		for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
		    x_cen[na] += xa[ns];
		    y_cen[na] += ya[ns];
		}
		x_cen[na] /= Ns[na];
		y_cen[na] /= Ns[na];
	}

	double 		x_cen_m, y_cen_m, x_cen_n, y_cen_n;
	double 		dx_cen, dy_cen;
	double 		dx_cen_shift, dy_cen_shift;
	double 		x1, y1, x3, y3;
	double 		lc1, lc2, dl_n, dl_m;

	VL_count = 0;

	for (na = 0; na < Na - 1; na++){
	    x_cen_n = x_cen[na];
	    y_cen_n = y_cen[na];
	    for (ma = na + 1; ma < Na; ma++){
	        x_cen_m = x_cen[ma];
	        y_cen_m = y_cen[ma];
	        dx_cen = x_cen_m - x_cen_n;
	        dy_cen = y_cen_m - y_cen_n;
			dy_cen_shift = std::round(dy_cen / Ly) * Ly;
			dy_cen -= dy_cen_shift;
	        if (std::sqrt(dx_cen * dx_cen + dy_cen * dy_cen) > Da[na] + Da[ma] + r_list)
	            continue;

		    for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
	            x1 = xa[ns]; y1 = ya[ns];
		    	for (ms = idx_list[ma]; ms < idx_list[ma + 1]; ms++){
	                dx = xa[ms] - x1;
					if (std::abs(dx) < r_list){
						dy = ya[ms] - y1 - dy_cen_shift;
						if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
							VL_aa[VL_count][0] = ns;
							VL_aa[VL_count][1] = ms;
							VL_aa[VL_count][2] = na;
							VL_aa[VL_count][3] = ma;
							VL_count++;
						}
					}
		    	}
		    }
	    }
	}

	// intra-cell interaction
	/*
	for (na = 0; na < Na; na++){
	    for (ns = idx_start[na]; ns < idx_end[na] - 1; ns++){
	        for (ms = ns + 2; ms <= idx_end[na]; ms++){
	        	if ((ns == idx_start[na]) && (ms == idx_end[na]))
		    		continue;
	            dx = pos_a[ms] - pos_a[ns];
	            dy = pos_a[ms + Nhalf] - pos_a[ns + Nhalf];
	            dr = dx * dx + dy * dy;
	            if (dr < r_list_sq){
					VL_aa[VL_counter[0]][0] = ns;
					VL_aa[VL_counter[0]][1] = ms;
					VL_aa[VL_counter[0]][2] = na;
					VL_aa[VL_counter[0]][3] = na;
					VL_counter[0]++;
				}
	        }
	    }
	}
	*/
	
	for (i = 0; i < Ns_tot; i++){
		xa_save[i] = xa[i];
        ya_save[i] = ya[i];
    }
}

void BreastCancer_2D::VerletList_DPM_Disk(double r_cut,
											vector<double> &xa_save,
											vector<double> &ya_save,
											vector<double> &xc_save,
											vector<double> &yc_save,
											vector<vector <int>> &VL_aa,
											vector<vector <int>> &VL_cc,
											vector<vector <int>> &VL_ac,
											vector<int> &VL_count,
											int first_call)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list = r_factor * r_cut;
	double 		r_list_sq = r_list * r_list;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dr;
	int 		i, iy, j, jy;
	int 		na, ma, ns, ms, nc, mc;

	if (first_call == 0){
		double 	dr_max = 0.0;
		for (i = 0; i < Ns_tot; i++){
			dx = xa[i] - xa_save[i];
			dy = ya[i] - ya_save[i];
			dy -= std::round(dy / Ly) * Ly;
			dr = dx * dx + dy * dy;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
		for (i = 0; i < Nc; i++){
			dx = xc[i] - xc_save[i];
			dy = yc[i] - yc_save[i];
			dy -= std::round(dy / Ly) * Ly;
			dr = dx * dx + dy * dy;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return;
	    }
	}

	for (i = 0; i < VL_count.size(); i++)
		VL_count[i] = 0;

	double 		x1, y1, x3, y3;

	for (na = 0; na < Na - 1; na++){
		for (ma = na + 1; ma < Na; ma++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				x1 = xa[ns];
				y1 = ya[ns];
				for (ms = idx_list[ma]; ms < idx_list[ma + 1]; ms++){
					dx = xa[ms] - x1;
					if (std::abs(dx) < r_list){
						dy = ya[ms] - y1;
						dy -= std::round(dy / Ly) * Ly;
						if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
							VL_aa[VL_count[0]][0] = ns;
							VL_aa[VL_count[0]][1] = ms;
							VL_aa[VL_count[0]][2] = na;
							VL_aa[VL_count[0]][3] = ma;
							VL_count[0]++;
						}
					}
				}
			}
		}
	}

	// cancer-cancer force
	for (nc = 0; nc < Nc - 1; nc++){
		for (mc = nc + 1; mc < Nc; mc++){
			dx = xc[mc] - xc[nc];
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc[nc];
				dy -= std::round(dy / Ly) * Ly;
				if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
					VL_cc[VL_count[1]][0] = mc;
					VL_cc[VL_count[1]][1] = nc;
					VL_count[1]++;
				}
			}
		}
	}

	// adipocyte-cancer force
	for (nc = 0; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xa[ns] - xc[nc];
				if (std::abs(dx) < r_list){
					dy = ya[ns] - yc[nc];
					dy -= std::round(dy / Ly) * Ly;
					if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
						VL_ac[VL_count[2]][0] = ns;
						VL_ac[VL_count[2]][1] = na;
						VL_ac[VL_count[2]][2] = nc;
						VL_count[2]++;
					}
				}
			}
		}
	}

	for (i = 0; i < Ns_tot; i++){
		xa_save[i] = xa[i];
		ya_save[i] = ya[i];
	}
	for (nc = 0; nc < Nc; nc++){
		xc_save[nc] = xc[nc];
		yc_save[nc] = yc[nc];
	}
}

void BreastCancer_2D::VerletList_DPM_Disk_Append(int Nc_old,
													double r_cut,
													vector<double> &xc_save,
													vector<double> &yc_save,
													vector<vector <int>> &VL_cc,
													vector<vector <int>> &VL_ac,
													vector<int> &VL_count)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list = r_factor * r_cut;
	double 		r_list_sq = r_list * r_list;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dr;
	int 		i, iy, j, jy;
	int 		na, ma, ns, ms, nc, mc;

	// cancer-cancer force
	for (nc = 0; nc < Nc_old; nc++){
		for (mc = Nc_old; mc < Nc; mc++){
			dx = xc[mc] - xc_save[nc];
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc_save[nc];
				dy -= std::round(dy / Ly) * Ly;
				if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
					VL_cc[VL_count[1]][0] = mc;
					VL_cc[VL_count[1]][1] = nc;
					VL_count[1]++;
				}
			}
		}
	}
	for (nc = Nc_old; nc < Nc - 1; nc++){
		for (mc = nc + 1; mc < Nc; mc++){
			dx = xc[mc] - xc[nc];
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc[nc];
				dy -= std::round(dy / Ly) * Ly;
				if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
					VL_cc[VL_count[1]][0] = mc;
					VL_cc[VL_count[1]][1] = nc;
					VL_count[1]++;
				}
			}
		}
	}

	// adipocyte-cancer force
	for (nc = Nc_old; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xa[ns] - xc[nc];
				if (std::abs(dx) < r_list){
					dy = ya[ns] - yc[nc];
					dy -= std::round(dy / Ly) * Ly;
					if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
						VL_ac[VL_count[2]][0] = ns;
						VL_ac[VL_count[2]][1] = na;
						VL_ac[VL_count[2]][2] = nc;
						VL_count[2]++;
					}
				}
			}
		}
	}

    for (nc = Nc_old; nc < Nc; nc++){
        xc_save[nc] = xc[nc];
        yc_save[nc] = yc[nc];
    }
}

void BreastCancer_2D::VerletList_Disk_Append(int Nc_old,
												double r_cut,
												vector<double> &xc_save,
												vector<double> &yc_save,
												vector<vector <int>> &VL_cc,
												int &VL_count)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list = r_factor * r_cut;
	double 		r_list_sq = r_list * r_list;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dr;
	int 		i, iy, j, jy;
	int 		nc, mc;

	// cancer-cancer force
	for (nc = 0; nc < Nc_old; nc++){
		for (mc = Nc_old; mc < Nc; mc++){
			dx = xc[mc] - xc_save[nc];
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc_save[nc];
				dy -= std::round(dy / Ly) * Ly;
				if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
					VL_cc[VL_count][0] = mc;
					VL_cc[VL_count][1] = nc;
					VL_count++;
				}
			}
		}
	}
	for (nc = Nc_old; nc < Nc - 1; nc++){
		for (mc = nc + 1; mc < Nc; mc++){
			dx = xc[mc] - xc[nc];
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc[nc];
				dy -= std::round(dy / Ly) * Ly;
				if ((std::abs(dy) < r_list) && (dx * dx + dy * dy < r_list_sq)){
					VL_cc[VL_count][0] = mc;
					VL_cc[VL_count][1] = nc;
					VL_count++;
				}
			}
		}
	}

    for (nc = Nc_old; nc < Nc; nc++){
        xc_save[nc] = xc[nc];
        yc_save[nc] = yc[nc];
    }
}

void BreastCancer_2D::Periphery_List_Initial(double r_cut,
												vector<int> &is_in_periphery)
{
	double 		r_cut_sq = r_cut * r_cut;
	double 		dx, dy, dr;
	int 		na, ns, nc;
	int			is_in;

	for (nc = 0; nc < Nc; nc++){
		is_in = 0;
		for (na = 0; na < Na; na++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xa[ns] - xc[nc];
				if (std::abs(dx) < r_cut){
					dy = ya[ns] - yc[nc];
					dy -= std::round(dy / Ly) * Ly;
					if ((std::abs(dy) < r_cut) && (dx * dx + dy * dy < r_cut_sq)){
						is_in = 1;
						is_in_periphery[nc] = 1;
					}
				}
				if (is_in == 1)
					break;
			}
			if (is_in == 1)
				break;
		}
	}
}

void BreastCancer_2D::Periphery_List_Check(double r_cut,
											vector<int> &is_in_periphery,
											int nc)
{
	double 		r_cut_sq = r_cut * r_cut;
	double 		dx, dy, dr;
	int 		na, ns;
	int			is_in = 0;

	for (na = 0; na < Na; na++){
		for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
			dx = xa[ns] - xc[nc];
			if (std::abs(dx) < r_cut){
				dy = ya[ns] - yc[nc];
				dy -= std::round(dy / Ly) * Ly;
				if ((std::abs(dy) < r_cut) && (dx * dx + dy * dy < r_cut_sq)){
					is_in = 1;
					is_in_periphery[nc] = 1;
				}
			}
			if (is_in == 1)
				break;
		}
		if (is_in == 1)
			break;
	}
	if (is_in == 0)
		is_in_periphery[nc] = 0;
}

void BreastCancer_2D::Periphery_List_PF_Initial(double phi_c,
												vector<double> &phi,
												vector<int> &is_in_periphery)
{
	for (int nc = 0; nc < Nc; nc++){
		if (phi[nc] < phi_c)
			is_in_periphery[nc] = 1;
	}
}

void BreastCancer_2D::Periphery_List_PF_Update(double phi_c,
												vector<double> &phi,
												vector<int> &is_in_periphery,
												vector<int> is_dividing,
												vector<int> is_not_dividing_ori)
{
	for (int nc = 0; nc < Nc; nc++){
		if ((is_not_dividing_ori[nc]) || (is_dividing[nc]))
			continue;
		if (phi[nc] < phi_c)
			is_in_periphery[nc] = 1;
	}
}

void BreastCancer_2D::Contact_AC(vector<vector <int>> &VL_ac,
									int VL_count,
									vector<int> &CC,
									vector<int> &reach_min)
{
    int 			ns, nc, na;
    double 			Dnm;
    double 			dx, dy, dnm;
	int 			vl_idx;

	for (na = 0; na < Na; na++)
		CC[na] = 0;

	for (int vl_idx = 0; vl_idx < VL_count; vl_idx++){
		na = VL_ac[vl_idx][1];
		if (reach_min[na] == 0){
			ns = VL_ac[vl_idx][0];
			nc = VL_ac[vl_idx][2];
			Dnm = 1.2 * R0[na] + Rc[nc];
			
			dx = xc[nc] - xa[ns];
			if (std::abs(dx) < Dnm){
				dy = yc[nc] - ya[ns];
				dy -= std::round(dy / Ly) * Ly;
				if ((std::abs(dy) < Dnm) && (std::sqrt(dx * dx + dy * dy) < Dnm))
					CC[na]++;
			}
		}
	}
}

void BreastCancer_2D::Contact_AC(vector<vector <int>> &VL_ac,
									int VL_count,
									vector<int> &CC,
									vector<int> &reach_min,
                                    vector<int> &is_ll)
{
    int 			ns, nc, na;
    double 			Dnm;
    double 			dx, dy, dnm;
	int 			vl_idx;

	for (na = 0; na < Na; na++)
		CC[na] = 0;

	for (int vl_idx = 0; vl_idx < VL_count; vl_idx++){
		na = VL_ac[vl_idx][1];
		if ((is_ll[na] == 0) || (reach_min[na] == 0)){
			ns = VL_ac[vl_idx][0];
			nc = VL_ac[vl_idx][2];
			Dnm = 1.2 * R0[na] + Rc[nc];
			
			dx = xc[nc] - xa[ns];
			if (std::abs(dx) < Dnm){
				dy = yc[nc] - ya[ns];
				dy -= std::round(dy / Ly) * Ly;
				if ((std::abs(dy) < Dnm) && (std::sqrt(dx * dx + dy * dy) < Dnm))
					CC[na]++;
			}
		}
	}
}
