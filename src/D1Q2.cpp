#include "../inc/D1Q2.hpp"

void D1Q2::export_results(const string &filename)
{
    ofstream outFile(filename);
    outFile << "x,T,Flux\n";

    for (int i = 0; i < m; i++)
    {
        double x_pos = i * dx;
        double flux = omega * (f1[i] - f2[i]) / dx;

        outFile << x_pos << "," << T[i] << "," << flux << "\n";
    }

    outFile.close();
    std::cout << "Successfully saved to: " << filename << std::endl;
}

#ifdef use_matlab

#include "mex.h"
void lbm_compute(int m, double alpha, double dt, double dx, int nstep, double Tw_left, double Tw_right, double *T_out)
{
    double omega = 1.0 / (alpha * dt / (dx * dx) + 0.5);

    std::vector<double> T(m, 0.5);
    std::vector<double> f1(m, 0.25);
    std::vector<double> f2(m, 0.25);

    for (int k = 0; k < nstep; ++k)
    {

        for (int i = 0; i < m; ++i)
        {
            double feq = 0.5 * T[i];
            f1[i] = (1.0 - omega) * f1[i] + omega * feq;
            f2[i] = (1.0 - omega) * f2[i] + omega * feq;
        }

        for (int i = m - 1; i > 0; --i)
            f1[i] = f1[i - 1];

        for (int i = 0; i < m - 1; ++i)
            f2[i] = f2[i + 1];

        f1[0] = Tw_left - f2[0];
        f2[m - 1] = Tw_right - f1[m - 1];

        for (int i = 0; i < m; ++i)
            T[i] = f1[i] + f2[i];
    }

    for (int i = 0; i < m; ++i)
        T_out[i] = T[i];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int m = (int)mxGetScalar(prhs[0]);
    double alpha = mxGetScalar(prhs[1]);
    double dt = mxGetScalar(prhs[2]);
    double dx = mxGetScalar(prhs[3]);
    int nstep = (int)mxGetScalar(prhs[4]);
    double Tw_left = mxGetScalar(prhs[5]);
    double Tw_right = mxGetScalar(prhs[6]);

    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    double *T_out = mxGetPr(plhs[0]);

    lbm_compute(m, alpha, dt, dx, nstep, Tw_left, Tw_right, T_out);
}
#endif

#ifdef D1Q2_T

void D1Q2::D1Q2_t()
{

    reset();
    cout << "Omega: " << omega << " dx: " << dx << endl;

    for (int k1 = 1; k1 <= nstep; ++k1)
    {

        for (int i = 0; i < m; ++i)
        {
            double feq = T[i] / 2;
            f1[i] = (1 - omega) * f1[i] + omega * feq;
            f2[i] = (1 - omega) * f2[i] + omega * feq;
        } // 碰撞

        for (int i = m - 1; i > 0; --i)
        {
            f1[i] = f1[i - 1];
        } // 迁移

        for (int i = 0; i < m - 1; ++i)
        {
            f2[i] = f2[i + 1];
        } // 迁移

        f1[0] = Tw_left - f2[0];
        f2[m - 1] = Tw_right - f1[m - 1];

        for (int i = 0; i < m; ++i)
        {
            T[i] = f1[i] + f2[i];
        }
    }

    D1Q2::export_results("T_D1Q2.csv");
}

#endif

#ifdef D1Q2_F

void D1Q2::D1Q2_f()
{

    reset();

    T.assign(m, 0.0);
    f1.assign(m, 0.0);
    f2.assign(m, 0.0);

    for (int k1 = 0; k1 < nstep; ++k1)
    {

        for (int i = 0; i < m; ++i)
        {
            double feq = 0.5 * T[i];
            f1[i] = (1 - omega) * f1[i] + omega * feq;
            f2[i] = (1 - omega) * f2[i] + omega * feq;
        } // 碰撞

        for (int i = m - 1; i > 0; --i)
        {
            f1[i] = f1[i - 1];
        } // 迁徙

        for (int i = 0; i < m - 1; ++i)
        {
            f2[i] = f2[i + 1];
        } // 迁徙

        f1[0] = f2[0] + Tw_left * dx / (omega * 2);
        f2[m - 1] = Tw_right - f1[m - 1];

        for (int j = 0; j < m; ++j)
        {
            T[j] = f1[j] + f2[j];
        }
    }

    D1Q2::export_results("F_D1Q2.csv");
}

#endif

D1Q2::~D1Q2()
{

    std::cout << "Object destroyed: D1Q3 simulation resources have been released." << std::endl;
}
