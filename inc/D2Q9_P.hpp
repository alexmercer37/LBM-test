#ifndef _D2Q9_P_H
#include "main.hpp"
class D2Q9_P
{
protected:
    int m, n, nstep;

    double dt, sigma, tau_e, tau_v, u_in, rho_out;

    vector<double> s_diag, f_post, f_local, m_local, meq, S_local;

    vector<int> ex = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    vector<int> ey = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    vector<double> w = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    vector<double> w_f = {0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0};

    vector<vector<int>> e = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};

    vector<vector<double>> Fx, Fy, rho, psi, u, v;
    vector<vector<vector<double>>> f, f_new;

public:
    D2Q9_P(int m_, int n_, int nstep_, double dt_) : m(m_), n(n_), nstep(nstep_), dt(dt_)
    {

        this->sigma = 1.0;
        this->tau_v = 1.0;
        this->tau_e = 1.1;

        this->u_in = 0.02;
        this->rho_out = 0.1;

        s_diag = {0.0, 1.0 / tau_e, 1.0 / tau_e, 0.0, 1.0 / tau_e, 0.0, 1.0 / tau_e, 1.0 / tau_v, 1.0 / tau_v};

        meq.assign(9, 0.0);
        f_post.assign(9, 0.0);

        f_local.assign(9, 0.0);
        m_local.assign(9, 0.0);
        S_local.assign(9, 0.0);

        v.assign(m, vector<double>(n, 0.0));
        u.assign(m, vector<double>(n, 0.0));

        Fx.assign(m, vector<double>(n, 0.0));
        Fy.assign(m, vector<double>(n, 0.0));
    }

    void D2Q9_p();
    void export_results(const string &filename);

    ~D2Q9_P();
};

class D2Q9_T : D2Q9_P
{
private:
    double alpha, Cv, omega, Tc, a, b, R, psi_wall, G;

    vector<vector<double>> k1, k2, k3, k4, T_tmp, current_T;

public:
    D2Q9_T(int m_, int n_, double dt_, int nstep_, double psi_wall_) : D2Q9_P(m_, n_, nstep_, dt_), psi_wall(psi_wall_)
    {

        this->R = 1.0;
        this->G = -1.0;
        this->Cv = 2.5;

        this->alpha = 0.02;
        this->omega = 0.344;

        this->a = 2.0 / 49.0;
        this->b = 2.0 / 21.0;
        this->Tc = 0.0778 * (a / (R * b));

        k1.assign(m, vector<double>(n, 0.0));
        current_T = k2 = k3 = k4 = T_tmp = k1;

        rho.assign(m, vector<double>(n, 0.0));
        psi.assign(m, vector<double>(n, 0.0));

        f.assign(9, vector<vector<double>>(m, vector<double>(n, 0.0)));
        f_new = f;

        init();
    }

    double get_dTdt(int j, int k, const std::vector<std::vector<double>> &current_T);
    double calculate_dpdT(double rho_val, double T_val);

    void compute_rho_u();
    void compute_i_foeces();
    void rk4_step();

    void init()
    {

        double rho_l = 1.8;           // 饱和液体密度
        double rho_g = 0.2;           // 饱和气体密度
        double T_init = 0.9 * Tc;     // 初始饱和温度
        double interface_y = n / 3.0; // 界面位置（底部1/3处）
        double width = 4.0;           // 界面平滑宽度（建议3.0-5.0）

        for (int k = 0; k < n; ++k)
        {

            double rho_val = 0.5 * (rho_l + rho_g) - 0.5 * (rho_l - rho_g) * tanh(2.0 * (k - interface_y) / width);

            for (int j = 0; j < m; ++j)
            {

                rho[j][k] = rho_val;
                psi[j][k] = rho_val; // P-R模型中 psi = rho

                current_T[j][k] = T_init;

                u[j][k] = 0.0;
                v[j][k] = 0.0;
            }
        }

        for (int i = 0; i < 9; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                for (int k = 0; k < n; ++k)
                {

                    double r = rho[j][k];
                    double ux = u[j][k];
                    double uy = v[j][k];

                    double u2 = ux * ux + uy * uy;
                    double eu = ex[i] * ux + ey[i] * uy;

                    f[i][j][k] = w[i] * r * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);

                    f_new[i][j][k] = f[i][j][k];
                }
            }
        }
    }
};

#endif