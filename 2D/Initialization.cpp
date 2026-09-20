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

bool BreastCancer_2D::overlap(vector<double> &xy)
{
	double dx, dy, dr, Dnm;
	
	for (int na = 0; na < Na - 1; na ++){
		for (int ma = na + 1; ma < Na; ma++){
			Dnm = Da[na] + Da[ma];
			dx = xy[ma] - xy[na];
			if (std::abs(dx) < Dnm){
				dy = xy[ma + Na] - xy[na + Na];
                dy -= std::round(dy / Ly) * Ly;
				if (std::abs(dy) < Dnm){
					dr = std::sqrt(dx * dx + dy * dy);
					if (dr < Dnm)
						return true;
				}
			}
		}
	}

	for (int na = 0; na < Na; na++)
		if ((xy[na] > Lx - Da[na]) || (xy[na] < Da[na]))
			return true;

	return false;

}

void BreastCancer_2D::AdipocyteInitialization(double ipf)
{
    Lx = std::sqrt(0.37 * PI * Na * std::sin(PI / (double)Nss) * std::sin(PI / (double)Nss) / D0_a / D0_a / ipf);
	Ly = Lx;
	
    std::mt19937 gen(pid);
    std::uniform_real_distribution<> dis_uniform(0.0, 1.0);

    vector<double>  Ccen(2 * Na, 0.0);
    while (true){
		for (int na = 0; na < Na; na++){
			Ccen[na] = (Lx - 2.0 * Da[0]) * dis_uniform(gen) + Da[0];
			Ccen[na + Na] = Ly * dis_uniform(gen);
		}
		if (!overlap(Ccen))
			break;
	}

    for (int na = 0; na < Na; na++){
		double d_theta = 2.0 * PI / (double)Ns[na];
		double theta;
		double theta_0 = dis_uniform(gen) * d_theta;
    	for (int ns = idx_list[na]; ns < idx_list[na + 1]; ns++){
			theta = theta_0 + (ns - idx_list[na]) * d_theta;
			xa[ns] = Ccen[na] + 0.5 * Da[na] * std::cos(theta);
			ya[ns] = Ccen[na + Na] + 0.5 * Da[na] * std::sin(theta);
		}
    }
}

void BreastCancer_2D::CancerInitialization(double ipf)
{
    std::mt19937 gen_r(pid); //Standard mersenne_twister_engine seeded with pid
	std::mt19937 gen(pid + 1000);
    std::uniform_real_distribution<> dis_uniform(0.0, 1.0);
	std::normal_distribution<> dis_norm(0.0, 1.0);
	double 	Gn_test, Gn_min = 10.0;
	double 	A_sum = 0.0;
    for (int i = 0; i < Nc; i++) {
		while (true){
			Gn_test = 1 + 0.15 * dis_norm(gen_r);
			if ((Gn_test >= Rs_ratio_min) && (Gn_test <= Rs_ratio_max)){
				Rc[i] = Gn_test * Rs; // Rs set in the header
				Rc_sq[i] = Rc[i] * Rc[i];
				Dc[i] = 2.0 * Rc[i];
				break;
			}
		}
		A_sum += Rc_sq[i];
		if (Gn_min > Gn_test)
			Gn_min = Gn_test;
    }	
	A_sum *= PI;
	Lx = A_sum / ipf / Ly; // square box length
	for (int i = 0; i < Nc; i++){
		xc[i] = dis_uniform(gen) * Lx;
		yc[i] = dis_uniform(gen) * Ly;
	}
}
