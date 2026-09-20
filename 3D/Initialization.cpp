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
#include "Bumpy_3D.h"

using std::vector;

void BreastCancer_3D::SetParticleParameters()
{
	double 	D0_unit, Alf_min, A_sum_unit;
	vector<double> 	A0_unit(Nf, 0.0);
    LoadUnitParameter(D0_unit, Alf_min, A0_unit, A_sum_unit);
	
	Alf *= Alf_min;

	double 	V0_single = 1.0;
    double 	Ra_single = std::pow(V0_single * Alf / (std::pow(A_sum_unit, 1.5) / std::sqrt(PI) / 6.0), 1.0 / 3.0);
	double 	D0_intra_single = D0_unit * Ra_single;
	double 	D0_single = delta * D0_intra_single;
	double 	D0_sq_single = D0_single * D0_single;
	double 	D0_intra_sq_single = D0_intra_single * D0_intra_single;
    double  Ra_sq = Ra_single * Ra_single;
	for (int na = 0; na < Na; na++){
		V0[na] = V0_single;
		Ra[na] = Ra_single;
		D0_intra[na] = D0_intra_single;
		R0[na] = 0.5 * D0_single;
		D0_intra_sq[na] = D0_intra_sq_single;
		for (int nf = 0; nf < Nf; nf++)
			A0[na][nf] = A0_unit[nf] * Ra_sq;
	}
}

bool BreastCancer_3D::overlap(vector<double> &xyz)
{
	double dx, dy, dz, dr, Dnm;
	
	for (int na = 0; na < Na - 1; na ++){
		for (int ma = na + 1; ma < Na; ma++){
			Dnm = Ra[na] + Ra[ma];
			dx = xyz[ma] - xyz[na];
			dx -= std::round(dx / Lx) * Lx;
			if (std::abs(dx) < Dnm){
				dy = xyz[ma + Na] - xyz[na + Na];
				if (std::abs(dy) < Dnm){
					dz = xyz[ma + 2 * Na] - xyz[na + 2 * Na];
					dz -= std::round(dz / Lz) * Lz;
					if (std::abs(dx) < Dnm){
						dr = std::sqrt(dx * dx + dy * dy + dz * dz);
						if (dr < Dnm)
							return true;
					}
				}
			}
		}
	}

	for (int na = 0; na < Na; na++)
		if ((xyz[na + Na] > Ly - Ra[na]) || (xyz[na + Na] < Ra[na]))
			return true;

	return false;

}

void BreastCancer_3D::AdipocyteInitialization(double ipf)
{
    Lx = std::pow(PI * Na / ipf, 1.0 / 3.0);
	Ly = Lx;
	Lz = Lx;
	
    std::mt19937 gen(pid);
    std::uniform_real_distribution<> dis_uniform(0.0, 1.0);

    vector<double>  Ccen(3 * Na, 0.0);
    while (true){
		for (int na = 0; na < Na; na++){
			Ccen[na] = Lx * dis_uniform(gen);
			Ccen[na + Na] = (Ly - 2.0 * Ra[0]) * dis_uniform(gen) + Ra[0];
			Ccen[na + 2 * Na] = Lz * dis_uniform(gen);
		}
		if (!overlap(Ccen))
			break;
	}

	vector<double> 	xyz_unit(3 * Ns, 0.0);
	LoadUnitPosition(xyz_unit);

	int i;
	double Ra_perturb;
    for (int na = 0; na < Na; na++){
    	for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
			Ra_perturb = (1.0 + (dis_uniform(gen) - 0.5) * 0.01) * Ra[na];
			i = ns - idx_list[na];
			xa[ns] = Ccen[na] + Ra_perturb * xyz_unit[3 * i];
			ya[ns] = Ccen[na + Na] + Ra_perturb * xyz_unit[3 * i + 1];
			za[ns] = Ccen[na + 2 * Na] + Ra_perturb * xyz_unit[3 * i + 2];
		}
    }
}

void BreastCancer_3D::AdipocyteInitialization_FromRigid(string dir_prefix)
{
	Bumpy_3D BumpyPacking(Na, Ns, delta, 1.0, 1.0, Pt_c, Pt_e, pid, dir_prefix);
	BumpyPacking.Packing_AboveJamming();
	Lx = BumpyPacking.ReturnBumpPos(xa, ya, za);
	Ly = Lx;
	Lz = Lx;
}

void BreastCancer_3D::CancerInitialization(double ipf)
{
    std::mt19937 gen_r(pid); //Standard mersenne_twister_engine seeded with pid
	std::mt19937 gen(pid + 1000);
    std::uniform_real_distribution<> dis_uniform(0.0, 1.0);
	std::normal_distribution<> dis_norm(0.0, 1.0);
	double 	Gn_test, Gn_min = 10.0;
	double 	V_sum = 0.0;
    for (int i = 0; i < Nc; i++) {
		while (true){
			Gn_test = 1 + 0.15 * dis_norm(gen_r);
			if ((Gn_test >= Rs_ratio_min) && (Gn_test <= Rs_ratio_max)){
				Rc[i] = Gn_test * Rs; // Rs set in the header
				Rc_sq[i] = Rc[i] * Rc[i];
				Dc[i] = 2.0 * Rc[i];
				Dc_cube[i] = Dc[i] * Dc[i] * Dc[i];
				break;
			}
		}
		V_sum += Rc_sq[i] * Rc[i];
		if (Gn_min > Gn_test)
			Gn_min = Gn_test;
    }	
	V_sum *= 4.0 / 3.0 * PI;
	Ly = V_sum / ipf / Lx / Lz; // square box length
	for (int i = 0; i < Nc; i++){
		xc[i] = dis_uniform(gen) * Lx;
		yc[i] = dis_uniform(gen) * Ly;
		zc[i] = dis_uniform(gen) * Lz;
	}
}