#ifndef _D1Q2_H

#include "main.hpp"

// #define D1Q2_F
// #define D1Q2_T

// #define use_matlab ;

class D1Q2
{
private:
    int m, nstep;
    double L, dx, dt, alpha, omega, Tw_left, Tw_right;

    vector<double> T;
    vector<double> f1;
    vector<double> f2;
    vector<double> x;

    void reset()
    {

        T.assign(m, 0.5);
        f1.assign(m, 0);
        f2.assign(m, 0);

        for (int i = 0; i < m; ++i)
        {
            x[i] = i * dx;
            f1[i] = T[i] / 2;
            f2[i] = T[i] / 2;
        }
    }

public:
    D1Q2(double L_, int m_, int nstep_, double dt_, double alpha_, double tl, double tr) : L(L_), m(m_), nstep(nstep_), dt(dt_), alpha(alpha_), Tw_left(tl), Tw_right(tr)
    {

        dx = L / (m - 1);
        omega = 1.0 / (alpha * dt / (dx * dx) + 0.5);

        x.resize(m);
        T.resize(m);
        f1.resize(m);
        f2.resize(m);

        reset();
    }

    void D1Q2_t();
    void D1Q2_f();
    void export_results(const std::string &filename);

    ~D1Q2();
};

#endif