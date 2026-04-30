#ifndef PIPE_PHASE_CHANGE_LBM_HPP
#define PIPE_PHASE_CHANGE_LBM_HPP

#include "main.hpp"

using std::vector;

class PipePhaseChangeLBM
{
public:
    int nx, ny, nstep;
    double dt, Twall, liquid_mass, vapor_mass;

    double tau_f;
    double tau_t;
    double tau_phi;

    double u_in;
    double T_in;

    double rho_l, rho_g;
    double Tsat;
    double Kphase;
    double latent;
    double cp;

    vector<vector<double>> rho;
    vector<vector<double>> ux;
    vector<vector<double>> uy;
    vector<vector<double>> T;
    vector<vector<double>> phi;
    vector<vector<double>> mdot;

    vector<vector<vector<double>>> f, f_new;
    vector<vector<vector<double>>> g, g_new;
    vector<vector<vector<double>>> h, h_new;

    int ex[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int ey[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    double w[9] = {4.0 / 9.0,
                   1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
                   1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};

public:
    PipePhaseChangeLBM(int nx_, int ny_, int nstep_, double dt_);

    void Initialize();
    void ComputeMdot();
    void UpdateRhoFromPhi();

    void CollisionF();
    void CollisionT();
    void CollisionPhi();

    void Streaming();
    void ApplyBoundary();

    void Macroscopic();

    void Run();

    void export_results(const string &filename);

    double feq(int i, double rho_, double ux_, double uy_);
    double geq(int i, double T_, double ux_, double uy_);
    double heq(int i, double phi_, double ux_, double uy_);
};

#endif