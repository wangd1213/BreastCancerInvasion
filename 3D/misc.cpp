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

template <typename T>
vector<size_t> BreastCancer_3D::sort_indices(const vector<T> &v)
{
  	// initialize original index locations
  	vector<size_t> idx(v.size());
  	iota(idx.begin(), idx.end(), 0);

  	// sort indexes based on comparing values in v
  	// using std::stable_sort instead of std::sort
  	// to avoid unnecessary index re-orderings
  	// when v contains elements of equal values 
  	std::stable_sort(idx.begin(), idx.end(),
       				 [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  	return idx;
}

double BreastCancer_3D::MaxForce_DPM()
{
    double Acc_abs;
    double Acc_max = std::abs(Fx_a[0]);
	for (int n = 1; n < Ns_tot; n++){
		Acc_abs = std::abs(Fx_a[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
    for (int n = 0; n < Ns_tot; n++){
		Acc_abs = std::abs(Fy_a[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	for (int n = 0; n < Ns_tot; n++){
		Acc_abs = std::abs(Fz_a[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}

    return Acc_max;
}

double BreastCancer_3D::MaxForce_Disk()
{
    double Acc_abs;
    double Acc_max = std::abs(Fx_c[0]);
	for (int n = 1; n < Nc; n++){
		Acc_abs = std::abs(Fx_c[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
    for (int n = 0; n < Nc; n++){
		Acc_abs = std::abs(Fy_c[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	for (int n = 0; n < Nc; n++){
		Acc_abs = std::abs(Fz_c[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}

    return Acc_max;
}

double BreastCancer_3D::MaxForce_All()
{
	double Acc_max = MaxForce_DPM();
    double Acc_max2 = MaxForce_Disk();
	return (Acc_max > Acc_max2) ? Acc_max : Acc_max2;
}

bool BreastCancer_3D::CheckInvasionBoundary()
{
	int 	N_contact = 0;
	for (int nc = 0; nc < Nc; nc++)
		if (yc[nc] > Ly - 1.5 * Rc[nc]) // expand cancer cell diameters to accommodate 3D
			N_contact++;
	if (N_contact * Rs * Rs * PI / Lx / Lz > 0.2)
		return true;
	else
		return false;
}
