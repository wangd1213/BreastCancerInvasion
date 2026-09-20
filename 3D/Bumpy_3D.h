#ifndef BUMPY_3D_H
#define BUMPY_3D_H

#include <vector>
#include <random>
#include <string>

using std::vector;
using std::string;

class Bumpy_3D
{
private:
    int N;
    int Ns;
    int Ns_tot;
    int pid;
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> Q0;
    vector<double> Q1;
    vector<double> Q2;
    vector<double> Q3;
    vector<double> x_bump;
    vector<double> y_bump;
    vector<double> z_bump;
    vector<double> X_unit;
    vector<double> Y_unit;
    vector<double> Z_unit;
    double R;
    double R_bump;
    double D_bump;
    double D_max;
    vector<int> idx_list;
    vector<double> Fx;
    vector<double> Fy;
    vector<double> Fz;
    vector<double> Tx;
    vector<double> Ty;
    vector<double> Tz;

    double delta;
    double K;
    double Kw;
    int Pt_c;
    int Pt_e;
    double Pt;
    
    double L;

    string          dir_run;
    string          P_str;
    string          dva_name;
    string          rigid_unit_name;
    string          rigid_name;
    
    const double PI = 3.14159265358979323846264338328;

public:
    // Constructor
    Bumpy_3D(int N_set, int Ns_set, double delta,
                double K_set, double Kw_set,
                int Pt_c_set, int Pt_e_set, int seed, string dir_pre);
    
    // Destructor
    ~Bumpy_3D() = default;
    
    // Initialization methods
    void Initialization(double ipf);
    bool overlap(vector<double> &x, vector<double> &y, vector<double> &z);

    // packing generation
    void Packing_AboveJamming(double ipf = 0.01);
    
    // energy minimization
    void FIRE_VL();

    // verlet list
    void VerletList(double r_cut, vector<double> &x_bump_save, vector<double> &y_bump_save,
                    vector<double> &z_bump_save, vector<vector <int>> &VL,
                    int &VL_count, int first_call);
    
    // forces
    void InterForce_VL(vector<vector <int>> &VL, int VL_count);
    void WallForce();
    void Force_VL(vector<vector <int>> &VL, int VL_count);
    void clearForces();
    double Energy();

    // miscellaneous
    void GetBumpPos();
    double MaxForce();
    void SetParticleParameters();

    // read and write configurations
    void SaveConfig();
    void LoadUnitParameter(double &D0_unit, double &Alf_min, double &A_sum_unit);
    void LoadUnitBumpy();
    void LoadConfig();
    double ReturnBumpPos(vector<double> &x_bump_return, vector<double> &y_bump_return, vector<double> &z_bump_return);
};

#endif