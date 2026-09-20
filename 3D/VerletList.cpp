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

void BreastCancer_3D::VerletList_Disk(double r_cut,
			                            vector<double> &xc_save,
                                        vector<double> &yc_save,
										vector<double> &zc_save,
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

	double 		dx, dy, dz, dr;
	int 		n, m;

	if (first_call == 0){
		double dr_max = 0.0;
		for (n = 0; n < Nc; n++){
			dx = xc[n] - xc_save[n];
			dy = yc[n] - yc_save[n];
			dz = zc[n] - zc_save[n];
			dx -= std::round(dx / Lx) * Lx;
			dz -= std::round(dz / Lz) * Lz;
			dr = dx * dx + dy * dy + dz * dz;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq)
	        return;
	}

	VL_count = 0;

	for (n = 0; n < Nc - 1; n++){
	    for (m = n + 1; m < Nc; m++){
			dx = xc[m] - xc[n];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < r_list){
				dy = yc[m] - yc[n];
				if (std::abs(dy) < r_list){
					dz = zc[m] - zc[n];
					dz -= std::round(dz / Lz) * Lz;
					if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
						VL_cc[VL_count][0] = n;
						VL_cc[VL_count][1] = m;
						VL_count++;
					}
				}
			}
		}
	}	

	for (n = 0; n < Nc; n++){
		xc_save[n] = xc[n];
        yc_save[n] = yc[n];
		zc_save[n] = zc[n];
    }
}

void BreastCancer_3D::VerletList_DPM(double r_cut,
										vector<double> &xa_save,
										vector<double> &ya_save,
										vector<double> &za_save,
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

	double 		dx, dy, dz, dr;
	int 		na, ma, ns, ms, nc, mc;

	if (first_call == 0){
		double dr_max = 0.0;
		for (ns = 0; ns < Ns_tot; ns++){
			dx = xa[ns] - xa_save[ns];
			dy = ya[ns] - ya_save[ns];
			dz = za[ns] - za_save[ns];
			dx -= std::round(dx / Lx) * Lx;
			dz -= std::round(dz / Lz) * Lz;
			dr = dx * dx + dy * dy + dz * dz;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return;
	    }
	}

	VL_count = 0;

	for (na = 0; na < Na - 1; na++){
	    for (ma = na + 1; ma < Na; ma++){
		    for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
		    	for (ms = idx_list[ma]; ms < idx_list[ma + 1]; ms++){
	                dx = xa[ms] - xa[ns];
					dx -= std::round(dx / Lx) * Lx;
					if (std::abs(dx) < r_list){
						dy = ya[ms] - ya[ns];
						if (std::abs(dy) < r_list){
							dz = za[ms] - za[ns];
							dz -= std::round(dz / Lz) * Lz;
							if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
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
	
	for (ns = 0; ns < Ns_tot; ns++){
		xa_save[ns] = xa[ns];
        ya_save[ns] = ya[ns];
		za_save[ns] = za[ns];
    }
}

void BreastCancer_3D::VerletList_DPM_Disk(double r_cut,
											vector<double> &xa_save,
											vector<double> &ya_save,
											vector<double> &za_save,
											vector<double> &xc_save,
											vector<double> &yc_save,
											vector<double> &zc_save,
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

	double 		dx, dy, dz, dr;
	int 		i, na, ma, ns, ms, nc, mc;

	if (first_call == 0){
		double 	dr_max = 0.0;
		for (i = 0; i < Ns_tot; i++){
			dx = xa[i] - xa_save[i];
			dy = ya[i] - ya_save[i];
			dz = za[i] - za_save[i];
			dx -= std::round(dx / Lx) * Lx;
			dz -= std::round(dz / Lz) * Lz;
			dr = dx * dx + dy * dy + dz * dz;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
		for (i = 0; i < Nc; i++){
			dx = xc[i] - xc_save[i];
			dy = yc[i] - yc_save[i];
			dz = zc[i] - zc_save[i];
			dx -= std::round(dx / Lx) * Lx;
			dz -= std::round(dz / Lz) * Lz;
			dr = dx * dx + dy * dy + dz * dz;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return;
	    }
	}

	for (i = 0; i < VL_count.size(); i++)
		VL_count[i] = 0;

	double 		x1, y1, z1;

	for (na = 0; na < Na - 1; na++){
		for (ma = na + 1; ma < Na; ma++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				x1 = xa[ns];
				y1 = ya[ns];
				z1 = za[ns];
				for (ms = idx_list[ma]; ms < idx_list[ma + 1]; ms++){
					dx = xa[ms] - x1;
					dx -= std::round(dx / Lx) * Lx;
					if (std::abs(dx) < r_list){
						dy = ya[ms] - y1;
						// dy -= std::round(dy / Ly) * Ly;
						if (std::abs(dy) < r_list){
							dz = za[ms] - z1;
							dz -= std::round(dz / Lz) * Lz;
							if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
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
	}

	// cancer-cancer force
	for (nc = 0; nc < Nc - 1; nc++){
		for (mc = nc + 1; mc < Nc; mc++){
			dx = xc[mc] - xc[nc];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc[nc];
				// dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < r_list){
					dz = zc[mc] - zc[nc];
					dz -= std::round(dz / Lz) * Lz;
					if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
						VL_cc[VL_count[1]][0] = mc;
						VL_cc[VL_count[1]][1] = nc;
						VL_count[1]++;
					}
				}
			}
		}
	}

	// adipocyte-cancer force
	for (nc = 0; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xa[ns] - xc[nc];
				dx -= std::round(dx / Lx) * Lx;
				if (std::abs(dx) < r_list){
					dy = ya[ns] - yc[nc];
					// dy -= std::round(dy / Ly) * Ly;
					if (std::abs(dy) < r_list){
						dz = za[ns] - zc[nc];
						dz -= std::round(dz / Lz) * Lz;
						if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
							VL_ac[VL_count[2]][0] = ns;
							VL_ac[VL_count[2]][1] = na;
							VL_ac[VL_count[2]][2] = nc;
							VL_count[2]++;
						}
					}
				}
			}
		}
	}

	for (i = 0; i < Ns_tot; i++){
		xa_save[i] = xa[i];
		ya_save[i] = ya[i];
		za_save[i] = za[i];
	}
	for (nc = 0; nc < Nc; nc++){
		xc_save[nc] = xc[nc];
		yc_save[nc] = yc[nc];
		zc_save[nc] = zc[nc];
	}
}

void BreastCancer_3D::VerletList_DPM_FixDisk(double r_cut,
												vector<double> &xa_save,
												vector<double> &ya_save,
												vector<double> &za_save,
												vector<vector <int>> &VL_aa,
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

	double 		dx, dy, dz, dr;
	int 		i, na, ma, ns, ms, nc, mc;

	if (first_call == 0){
		double 	dr_max = 0.0;
		for (i = 0; i < Ns_tot; i++){
			dx = xa[i] - xa_save[i];
			dy = ya[i] - ya_save[i];
			dz = za[i] - za_save[i];
			dx -= std::round(dx / Lx) * Lx;
			dz -= std::round(dz / Lz) * Lz;
			dr = dx * dx + dy * dy + dz * dz;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return;
	    }
	}

	for (i = 0; i < VL_count.size(); i++)
		VL_count[i] = 0;

	double 		x1, y1, z1;

	for (na = 0; na < Na - 1; na++){
		for (ma = na + 1; ma < Na; ma++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				x1 = xa[ns];
				y1 = ya[ns];
				z1 = za[ns];
				for (ms = idx_list[ma]; ms < idx_list[ma + 1]; ms++){
					dx = xa[ms] - x1;
					dx -= std::round(dx / Lx) * Lx;
					if (std::abs(dx) < r_list){
						dy = ya[ms] - y1;
						// dy -= std::round(dy / Ly) * Ly;
						if (std::abs(dy) < r_list){
							dz = za[ms] - z1;
							dz -= std::round(dz / Lz) * Lz;
							if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
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
	}

	// adipocyte-cancer force
	for (nc = 0; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xa[ns] - xc[nc];
				dx -= std::round(dx / Lx) * Lx;
				if (std::abs(dx) < r_list){
					dy = ya[ns] - yc[nc];
					// dy -= std::round(dy / Ly) * Ly;
					if (std::abs(dy) < r_list){
						dz = za[ns] - zc[nc];
						dz -= std::round(dz / Lz) * Lz;
						if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
							VL_ac[VL_count[1]][0] = ns;
							VL_ac[VL_count[1]][1] = na;
							VL_ac[VL_count[1]][2] = nc;
							VL_count[1]++;
						}
					}
				}
			}
		}
	}

	for (i = 0; i < Ns_tot; i++){
		xa_save[i] = xa[i];
		ya_save[i] = ya[i];
		za_save[i] = za[i];
	}
}

void BreastCancer_3D::VerletList_Disk_FixDPM(double r_cut,
												vector<double> &xc_save,
												vector<double> &yc_save,
												vector<double> &zc_save,
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

	double 		dx, dy, dz, dr;
	int 		i, na, ma, ns, ms, nc, mc;

	if (first_call == 0){
		double 	dr_max = 0.0;
		for (i = 0; i < Nc; i++){
			dx = xc[i] - xc_save[i];
			dy = yc[i] - yc_save[i];
			dz = zc[i] - zc_save[i];
			dx -= std::round(dx / Lx) * Lx;
			dz -= std::round(dz / Lz) * Lz;
			dr = dx * dx + dy * dy + dz * dz;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return;
	    }
	}

	for (i = 0; i < VL_count.size(); i++)
		VL_count[i] = 0;

	double 		x1, y1, z1;

	// cancer-cancer force
	for (nc = 0; nc < Nc - 1; nc++){
		for (mc = nc + 1; mc < Nc; mc++){
			dx = xc[mc] - xc[nc];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc[nc];
				// dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < r_list){
					dz = zc[mc] - zc[nc];
					dz -= std::round(dz / Lz) * Lz;
					if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
						VL_cc[VL_count[0]][0] = mc;
						VL_cc[VL_count[0]][1] = nc;
						VL_count[0]++;
					}
				}
			}
		}
	}

	// adipocyte-cancer force
	for (nc = 0; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xa[ns] - xc[nc];
				dx -= std::round(dx / Lx) * Lx;
				if (std::abs(dx) < r_list){
					dy = ya[ns] - yc[nc];
					// dy -= std::round(dy / Ly) * Ly;
					if (std::abs(dy) < r_list){
						dz = za[ns] - zc[nc];
						dz -= std::round(dz / Lz) * Lz;
						if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
							VL_ac[VL_count[1]][0] = ns;
							VL_ac[VL_count[1]][1] = na;
							VL_ac[VL_count[1]][2] = nc;
							VL_count[1]++;
						}
					}
				}
			}
		}
	}

	for (nc = 0; nc < Nc; nc++){
		xc_save[nc] = xc[nc];
		yc_save[nc] = yc[nc];
		zc_save[nc] = zc[nc];
	}
}

void BreastCancer_3D::VerletList_DPM_Disk(double r_cut,
											vector<double> &xa_save,
											vector<double> &ya_save,
											vector<double> &za_save,
											vector<double> &xc_save,
											vector<double> &yc_save,
											vector<double> &zc_save,
											vector<vector <int>> &VL_ac,
											int &VL_count,
											int first_call)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list = r_factor * r_cut;
	double 		r_list_sq = r_list * r_list;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dz, dr;
	int 		i, na, ma, ns, ms, nc, mc;

	if (first_call == 0){
		double 	dr_max = 0.0;
		for (i = 0; i < Ns_tot; i++){
			dx = xa[i] - xa_save[i];
			dy = ya[i] - ya_save[i];
			dz = za[i] - za_save[i];
			dx -= std::round(dx / Lx) * Lx;
			dz -= std::round(dz / Lz) * Lz;
			dr = dx * dx + dy * dy + dz * dz;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
		for (i = 0; i < Nc; i++){
			dx = xc[i] - xc_save[i];
			dy = yc[i] - yc_save[i];
			dz = zc[i] - zc_save[i];
			dx -= std::round(dx / Lx) * Lx;
			dz -= std::round(dz / Lz) * Lz;
			dr = dx * dx + dy * dy + dz * dz;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return;
	    }
	}

	VL_count = 0;

	// adipocyte-cancer force
	for (nc = 0; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xa[ns] - xc[nc];
				dx -= std::round(dx / Lx) * Lx;
				if (std::abs(dx) < r_list){
					dy = ya[ns] - yc[nc];
					// dy -= std::round(dy / Ly) * Ly;
					if (std::abs(dy) < r_list){
						dz = za[ns] - zc[nc];
						dz -= std::round(dz / Lz) * Lz;
						if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
							VL_ac[VL_count][0] = ns;
							VL_ac[VL_count][1] = na;
							VL_ac[VL_count][2] = nc;
							VL_count++;
						}
					}
				}
			}
		}
	}

	for (i = 0; i < Ns_tot; i++){
		xa_save[i] = xa[i];
		ya_save[i] = ya[i];
		za_save[i] = za[i];
	}
	for (nc = 0; nc < Nc; nc++){
		xc_save[nc] = xc[nc];
		yc_save[nc] = yc[nc];
		zc_save[nc] = zc[nc];
	}
}

void BreastCancer_3D::VerletList_DPM_Disk_Append(int Nc_old,
													double r_cut,
													vector<double> &xc_save,
													vector<double> &yc_save,
													vector<double> &zc_save,
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

	double 		dx, dy, dz, dr;
	int 		na, ma, ns, ms, nc, mc;

	// cancer-cancer force
	for (nc = 0; nc < Nc_old; nc++){
		for (mc = Nc_old; mc < Nc; mc++){
			dx = xc[mc] - xc_save[nc];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc_save[nc];
				//dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < r_list){
					dz = zc[mc] - zc_save[nc];
					dz -= std::round(dz / Lz) * Lz;
					if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
						VL_cc[VL_count[1]][0] = mc;
						VL_cc[VL_count[1]][1] = nc;
						VL_count[1]++;
					}
				}
			}
		}
	}
	for (nc = Nc_old; nc < Nc - 1; nc++){
		for (mc = nc + 1; mc < Nc; mc++){
			dx = xc[mc] - xc[nc];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc[nc];
				// dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < r_list){
					dz = zc[mc] - zc[nc];
					dz -= std::round(dz / Lz) * Lz;
					if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
						VL_cc[VL_count[1]][0] = mc;
						VL_cc[VL_count[1]][1] = nc;
						VL_count[1]++;
					}
				}
			}
		}
	}

	// adipocyte-cancer force
	for (nc = Nc_old; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xa[ns] - xc[nc];
				dx -= std::round(dx / Lx) * Lx;
				if (std::abs(dx) < r_list){
					dy = ya[ns] - yc[nc];
					// dy -= std::round(dy / Ly) * Ly;
					if (std::abs(dy) < r_list){
						dz = za[ns] - zc[nc];
						dz -= std::round(dz / Lz) * Lz;
						if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
							VL_ac[VL_count[2]][0] = ns;
							VL_ac[VL_count[2]][1] = na;
							VL_ac[VL_count[2]][2] = nc;
							VL_count[2]++;
						}
					}
				}
			}
		}
	}

    for (nc = Nc_old; nc < Nc; nc++){
        xc_save[nc] = xc[nc];
        yc_save[nc] = yc[nc];
		zc_save[nc] = zc[nc];
    }
}

void BreastCancer_3D::VerletList_DPM_Disk_Append(int Nc_old,
													double r_cut,
													vector<double> &xc_save,
													vector<double> &yc_save,
													vector<double> &zc_save,
													vector<vector <int>> &VL_ac,
													int &VL_count)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list = r_factor * r_cut;
	double 		r_list_sq = r_list * r_list;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dz, dr;
	int 		na, ma, ns, ms, nc, mc;

	// adipocyte-cancer force
	for (nc = Nc_old; nc < Nc; nc++){
		for (na = 0; na < Na; na++){
			for (ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
				dx = xa[ns] - xc[nc];
				dx -= std::round(dx / Lx) * Lx;
				if (std::abs(dx) < r_list){
					dy = ya[ns] - yc[nc];
					// dy -= std::round(dy / Ly) * Ly;
					if (std::abs(dy) < r_list){
						dz = za[ns] - zc[nc];
						dz -= std::round(dz / Lz) * Lz;
						if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
							VL_ac[VL_count][0] = ns;
							VL_ac[VL_count][1] = na;
							VL_ac[VL_count][2] = nc;
							VL_count++;
						}
					}
				}
			}
		}
	}
}

void BreastCancer_3D::VerletList_Disk_Append(int Nc_old,
												double r_cut,
												vector<double> &xc_save,
												vector<double> &yc_save,
												vector<double> &zc_save,
												vector<vector <int>> &VL_cc,
												int &VL_count)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list = r_factor * r_cut;
	double 		r_list_sq = r_list * r_list;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dz, dr;
	int 		nc, mc;

	// cancer-cancer force
	for (nc = 0; nc < Nc_old; nc++){
		for (mc = Nc_old; mc < Nc; mc++){
			dx = xc[mc] - xc_save[nc];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc_save[nc];
				// dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < r_list){
					dz = zc[mc] - zc_save[nc];
					dz -= std::round(dz / Lz) * Lz;
					if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
						VL_cc[VL_count][0] = mc;
						VL_cc[VL_count][1] = nc;
						VL_count++;
					}
				}
			}
		}
	}
	for (nc = Nc_old; nc < Nc - 1; nc++){
		for (mc = nc + 1; mc < Nc; mc++){
			dx = xc[mc] - xc[nc];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < r_list){
				dy = yc[mc] - yc[nc];
				// dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < r_list){
					dz = zc[mc] - zc[nc];
					dz -= std::round(dz / Lz) * Lz;
					if ((std::abs(dz) < r_list) && (dx * dx + dy * dy + dz * dz < r_list_sq)){
						VL_cc[VL_count][0] = mc;
						VL_cc[VL_count][1] = nc;
						VL_count++;
					}
				}
			}
		}
	}

    for (nc = Nc_old; nc < Nc; nc++){
        xc_save[nc] = xc[nc];
        yc_save[nc] = yc[nc];
		zc_save[nc] = zc[nc];
    }
}

void BreastCancer_3D::Periphery_List_PF_Initial(double phi_c,
												vector<double> &phi,
												vector<int> &is_in_periphery)
{
	for (int nc = 0; nc < Nc; nc++){
		if (phi[nc] < phi_c)
			is_in_periphery[nc] = 1;
	}
}

void BreastCancer_3D::Periphery_List_PF_Update(double phi_c,
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

void BreastCancer_3D::Contact_AC(vector<vector <int>> &VL_ac,
									int VL_count,
									vector<int> &CC,
									vector<int> &reach_min)
{
    int 			ns, nc, na;
    double 			Dnm;
    double 			dx, dy, dz, dnm;
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
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < Dnm){
				dy = yc[nc] - ya[ns];
				// dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < Dnm){
					dz = zc[nc] - za[ns];
					dz -= std::round(dz / Lz) * Lz;
					if ((std::abs(dz) < Dnm) && (std::sqrt(dx * dx + dy * dy + dz * dz) < Dnm))
						CC[na]++;
				}
			}
		}
	}
}

void BreastCancer_3D::Contact_AC(vector<vector <int>> &VL_ac,
									int VL_count,
									vector<int> &CC,
									vector<int> &reach_min,
                                    vector<int> &is_ll)
{
    int 			ns, nc, na;
    double 			Dnm;
    double 			dx, dy, dz, dnm;
	int 			vl_idx;

	for (na = 0; na < Na; na++)
		CC[na] = 0;

	for (int vl_idx = 0; vl_idx < VL_count; vl_idx++){
		na = VL_ac[vl_idx][1];
		if ((reach_min[na] == 0) || (is_ll[na] == 0)){
			ns = VL_ac[vl_idx][0];
			nc = VL_ac[vl_idx][2];
			Dnm = 1.2 * R0[na] + Rc[nc];
			
			dx = xc[nc] - xa[ns];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < Dnm){
				dy = yc[nc] - ya[ns];
				// dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < Dnm){
					dz = zc[nc] - za[ns];
					dz -= std::round(dz / Lz) * Lz;
					if ((std::abs(dz) < Dnm) && (std::sqrt(dx * dx + dy * dy + dz * dz) < Dnm))
						CC[na]++;
				}
			}
		}
	}
}
