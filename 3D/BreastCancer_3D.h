#ifndef BREASTCANCER_3D_H
#define BREASTCANCER_3D_H

#include <vector>
#include <random>
#include <string>

using std::vector;
using std::string;

class BreastCancer_3D
{
private:
    int Na;
    int Ns;
    int Nf;
    int Ne;
    int Nnn;
    int Ns_tot;
    int Nc;
    int pid;
    vector<double> xa;
    vector<double> ya;
    vector<double> za;
    vector<double> xa_pin;
    vector<double> ya_pin;
    vector<double> za_pin;
    vector<double> xa_cen_pin;
    vector<double> ya_cen_pin;
    vector<double> za_cen_pin;
    vector<double> xc;
    vector<double> yc;
    vector<double> zc;
    vector<double> Ra;
    vector<double> R0;
    vector<double> D0_intra;
    vector<double> D0_intra_sq;
    vector<double> V0;
    vector<vector <double>> A0;
    vector<double> theta0;
    vector<int> idx_list;
    vector<vector <int>> edgelist;
    vector<vector <int>> f_unit;
    vector<vector <int>> nnlist;
    vector<double> Rc;
    vector<double> Dc;
    vector<double> Rc_sq;
    vector<double> Dc_cube;
    vector<double> Fx_a;
    vector<double> Fy_a;
    vector<double> Fz_a;
    vector<double> Fx_c;
    vector<double> Fy_c;
    vector<double> Fz_c;
    double Fw;
    double Vw;

    double delta;
    double Alf;
    double KV;
    double KA;
    double Kaa;
    double Kb;
    double Kb_std;
    vector<double> Kb_all;
    double Kw_a;
    double Kpin;
    double Kcc;
    double Kcc_attr;
    double Kac_attr;
    double Kac;
    double Kw_c;
    double a_cut_c;
    double a_cut_h;
    double a_cut_ac;
    double a_cut_ac_h;
    int Pt_c;
    int Pt_e;
    double Pt;
    
    double Lx;
    double Ly;
    double Lz;

    double Rs = 0.15;
    double Rs_ratio_max = 1.25;
    double Rs_ratio_min = 0.75;

    string          dir_run;
    string          P_str;
    string          Kr_str_c;
    string          Kr_str_ac;
    string          K_adip_str;
    string          dva_name;
    string          f_unit_name;
    string          edgelist_name;
    string          nnn_list_name;
    string          theta_name;
    string          unit_pos_name;
    string          rigid_unit_name;
    string          rigid_name;
    string          adip_name;
    string          cancer_name;
    string          all_name;
    string          pos_inv_name;
    string          vel_inv_name;
    
    const double PI = 3.14159265358979323846264338328;

public:
    // Constructor
    BreastCancer_3D(int Na_set, int Ns_set, int Nc_set, double delta,
					int dAlf_c_set, int dAlf_e_set, double KV_set, double KA_set,
                    double Kaa_set, int Kb_c_set, int Kb_e_set, double Kb_std_set,
                    double Kw_a_set, double Kpin_set, double Kcc_set,
                    double a_cut_c_set, double Kcc_attr_ratio, double a_cut_ac_set,
                    double Kac_attr_ratio, double Kac_set, double Kw_c_set,
                    int Pt_c_set, int Pt_e_set, int seed, string dir_pre);
    
    // Destructor
    ~BreastCancer_3D() = default;
    
    // Initialization methods
    void CancerInitialization(double ipf = 0.01);
    void AdipocyteInitialization(double ipf = 0.01);
    void AdipocyteInitialization_FromRigid(string dir_pre);
    bool overlap(vector<double> &xyz);

    // packing generation
    void AdipocytePacking_FromRigid(string dir_pre);
    void AdipocytePacking_Compress(string dir_pre, int P_start_c, int P_start_e,
                                    int P_end_c, int P_end_e, int P_num);
    void AdipocytePacking();
    void CancerPacking();
    void AllPacking();
    void AllPacking_Immediate();
    void AllPacking_ThreeSteps();
    void AllPacking_Compress();
    void AllPacking_DeCompress();

    // molecular dynamics
    void Invasion_Periphery_PF_Growth(double phi_c, int R_divide_c, int R_divide_e,
                                        double gamma_a_ratio, double gamma_c_ratio, double gamma_ac_ratio,
                                        double dt, int save_every_step, int Nt);
    void Invasion_Periphery_PF_Growth_LipidLoss(double phi_c, int T_divide_c, int T_divide_e,
                                                int R_lip_c, int R_lip_e, double gamma_a_ratio,
                                                double gamma_c_ratio, double gamma_ac_ratio,
                                                double dt, int save_every_step, int Nt);
    void Invasion_Periphery_PF_Growth_OD(double phi_c, int R_divide_c, int R_divide_e,
                                            double dt, int save_every_step, int Nt);
    void Invasion_Periphery_PF_Growth_LipidLoss_OD(double phi_c, int T_divide_c, int T_divide_e,
                                                    int R_lip_c, int R_lip_e, double dt,
                                                    int save_every_step, int Nt);
    void Invasion_Periphery_PF_Growth_LipidLossLinear_OD(double phi_c, int T_divide_c, int T_divide_e,
                                                            int R_lip_c, int R_lip_e, double dt,
                                                            int save_every_step, int Nt);
    
    // energy minimization
    void FIRE_DPM_VL();
    void FIRE_DPM_Enthalpy_VL(double Pt);
    void FIRE_Disk_VL();
    void FIRE_Disk_Enthalpy_VL(double Pt);
    void FIRE_DPM_Disk_VL();
    void FIRE_DPM_Disk_Enthalpy_VL(double Pt);
    void FIRE_DPM_FixDisk_Enthalpy_VL(double Pt);
    void FIRE_Disk_FixDPM_Enthalpy_VL(double Pt);

    // verlet list
    void VerletList_DPM(double r_cut, vector<double> &xa_save, vector<double> &ya_save,
                        vector<double> &za_save, vector<vector <int>> &VL_aa,
                        int &VL_count, int first_call);
    void VerletList_Disk(double r_cut, vector<double> &xc_save, vector<double> &yc_save,
			                vector<double> &zc_save, vector<vector <int>> &VL_cc,
                            int &VL_count, int first_call);
    void VerletList_DPM_Disk(double r_cut, vector<double> &xa_save, vector<double> &ya_save,
                                vector<double> &za_save, vector<double> &xc_save,
                                vector<double> &yc_save, vector<double> &zc_save,
                                vector<vector <int>> &VL_aa, vector<vector <int>> &VL_cc,
                                vector<vector <int>> &VL_ac, vector<int> &VL_count,
                                int first_call);
    void VerletList_DPM_FixDisk(double r_cut, vector<double> &xa_save, vector<double> &ya_save,
                                vector<double> &za_save, vector<vector <int>> &VL_aa,
                                vector<vector <int>> &VL_ac, vector<int> &VL_count,
                                int first_call);
    void VerletList_Disk_FixDPM(double r_cut, vector<double> &xc_save, vector<double> &yc_save,
                                vector<double> &zc_save, vector<vector <int>> &VL_cc,
                                vector<vector <int>> &VL_ac, vector<int> &VL_count,
                                int first_call);
    void VerletList_DPM_Disk(double r_cut, vector<double> &xa_save, vector<double> &ya_save,
                                vector<double> &za_save, vector<double> &xc_save,
                                vector<double> &yc_save, vector<double> &zc_save,
                                vector<vector <int>> &VL_ac, int &VL_count,
                                int first_call);
    void VerletList_DPM_Disk_Append(int Nc_old, double r_cut, vector<double> &xc_save,
                                    vector<double> &yc_save, vector<double> &zc_save,
                                    vector<vector <int>> &VL_cc, vector<vector <int>> &VL_ac,
                                    vector<int> &VL_count);
    void VerletList_DPM_Disk_Append(int Nc_old, double r_cut, vector<double> &xc_save,
                                    vector<double> &yc_save, vector<double> &zc_save,
                                    vector<vector <int>> &VL_ac, int &VL_count);
    void VerletList_Disk_Append(int Nc_old, double r_cut, vector<double> &xc_save,
                                vector<double> &yc_save, vector<double> &zc_save,
                                vector<vector <int>> &VL_cc, int &VL_count);
    void Contact_AC(vector<vector <int>> &VL_ac, int VL_count, vector<int> &CC, vector<int> &reach_min);
    void Contact_AC(vector<vector <int>> &VL_ac, int VL_count, vector<int> &CC, vector<int> &reach_min, vector<int> &is_ll);
    void Periphery_List_Initial(double r_cut, vector<int> &is_in_periphery);
    void Periphery_List_Check(double r_cut, vector<int> &is_in_periphery, int nc);
    void Periphery_List_PF_Initial(double phi_c, vector<double> &phi, vector<int> &is_in_periphery);
    void Periphery_List_PF_Update(double phi_c, vector<double> &phi, vector<int> &is_in_periphery,
								    vector<int> is_dividing, vector<int> is_dividing_ori);

    // local packing fraction
    double LocalPF(double R_cut, double R_cut_sq, double R, double D, double R_sq,
                    double D_cube, double dnm, double dnm_sq);
    void LocalPhi(vector<double> &xc, vector<double> &yc, vector<double> &zc,
                    double r_shell, vector<vector <int>> &VL_phi_cc,
                    int VL_phi_count, vector<double> &phi);
    void LocalPhi_Update(vector<double> &xc, vector<double> &yc, vector<double> &zc,
                            double r_shell, vector<vector <int>> &VL_phi_cc,
                            int VL_phi_count, vector<int> &is_dividing,
                            vector<int> &is_dividing_ori, vector<double> &phi);
    
    // forces
    void ShapeForce_DPM();
    void IntraForce_DPM();
    void PinningForce_DPM_Center();
    void PinningForce_DPM_Vertex();
    void InterForce_DPM();
    void Force_DPM();
    void InterForce_DPM_VL(vector<vector <int>> &VL_aa, int VL_count);
    double WallForce_DPM();
    void Force_DPM_VL(vector<vector <int>> &VL_aa, int VL_count);
    void InterForce_Disk();
    void InterForce_Disk_Test();
    void InterForce_Disk_VL(vector<vector <int>> &VL_cc, int VL_count);
    double WallForce_Disk();
    void Force_Disk_VL(vector<vector <int>> &VL_cc, int VL_count);
    void InterForce_DPM_Disk_VL(vector<vector <int>> &VL_ac, int VL_count);
    void InterForce_DPM_Disk();
    void Force_DPM_Disk_VL(vector<vector <int>> &VL_aa, vector<vector <int>> &VL_cc,
                            vector<vector <int>> &VL_ac, vector<int> &VL_count);
    void Force_DPM_FixDisk_VL(vector<vector <int>> &VL_aa, vector<vector <int>> &VL_ac, vector<int> &VL_count);
    void Force_Disk_FixDPM_VL(vector<vector <int>> &VL_cc, vector<vector <int>> &VL_ac, vector<int> &VL_count);
    void Force_DPM_Disk_Pin_VL(vector<vector <int>> &VL_aa, vector<vector <int>> &VL_cc,
                                vector<vector <int>> &VL_ac, vector<int> &VL_count);
    void Force_DPM_Disk_PinCen_VL(vector<vector <int>> &VL_aa, vector<vector <int>> &VL_cc,
                                    vector<vector <int>> &VL_ac, vector<int> &VL_count);
    void PinningForce_DPM_Center_LipidLoss(vector<int> &reach_min, vector<double> &rsc_all, double rsc_min);
    void clearForces_DPM();
    void clearForces_Cancer();
    
    // Utility methods
    double Pressure_DPM();
    double Pressure_Disk();
    double Pressure_DPM_Disk();

    // miscellaneous
    template <typename T> vector<size_t> sort_indices(const vector<T> &v);
    double MaxForce_DPM();
    double MaxForce_Disk();
    double MaxForce_All();
    void SetDamping(double gamma_a_ratio, double gamma_c_ratio, double gamma_ac_ratio,
                    double &gamma_a, double &gamma_c, double &gamma_ac);
    void SetOverDamping(double &gamma_a, double &gamma_c, double &gamma_w);
    bool CheckInvasionBoundary();
    void SetParticleParameters();
    void Adipocyte_Check();
    void Bumpy_Test(string dir_pre);

    // read and write configurations
    void SaveConfig_Adipocyte(string filename);
    void SaveConfig_Adipocyte_Comp(string filename);
    void SaveConfig_Cancer(string filename);
    void SaveConfig_All(string filename);
    void LoadUnitParameter(double &D0_unit, double &Alf_min,
							vector<double> &A0_unit,
							double &A_sum_unit);
    void LoadUnitPosition(vector<double> &xyz_unit);
    void Load_BumpyToDPM();
    void LoadConfig_Adipocyte(string filename, bool full_load = 1);
    void LoadConfig_Cancer(string filename);
    void LoadConfig_All(string filename);
    void SaveInv_Growth(string posname, string velname,
                        vector<double> &Vx_a, vector<double> &Vy_a,
                        vector<double> &Vz_a, vector<double> &Vx_c,
                        vector<double> &Vy_c, vector<double> &Vz_c);
    void SaveInv_Growth_LipidLoss(string posname, string velname, vector<double> &rsc_all,
                                    vector<double> &Vx_a, vector<double> &Vy_a,
                                    vector<double> &Vz_a, vector<double> &Vx_c,
                                    vector<double> &Vy_c, vector<double> &Vz_c);
    void SaveInv_Growth_OD(string posname);
    void SaveInv_Growth_LipidLoss_OD(string posname, vector<double> &rsc_all);
};

#endif
